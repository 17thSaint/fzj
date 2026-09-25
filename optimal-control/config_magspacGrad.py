"""QuOCS dCRAB optimization of the magnetic spacing ramp, with the gradient B(t) as the control.

Same transport problem as config_magspacRamp.py -- dd groundstate manifold at a_start -> the one
at a_end over ramptime T -- but instead of shaping a(t) and reading the gradient back, the
gradient itself is optimized and a(t) is its driven response in the trap,

    a(t) = a_start + (1/w) int_0^t B(t') sin(w(t - t')) dt'

(magspac-ramp-control-functions.jl: magspac_spacing_from_gradient). Shaping a(t) directly let
dCRAB buy fidelity by shaking the lattice: read-back |B| of 3300 unconstrained, and still 8231
with a low-passed basis, because the Nelder-Mead coefficients are unbounded and the a(t) bounds
clip into slope kinks (delta kicks in B). With B as the control the gradient cap is a hard
QuOCS amplitude limit, |B| <= B_MAX, honoured exactly by every pulse the search evaluates.

What that gives up is the end condition, which pinning a(t) used to provide for free. Here:
  - the update is scaled by (1 - t/T), freezing B(T) at the guess's final value, the holding
    gradient w^2 (a_end - a_start) that keeps a displaced mode parked at a_end;
  - FoM = fidelity - R, with R = sqrt((a(T) - a_end)^2 + (adot(T)/w)^2) the amplitude the mode
    rings with after T under that holding gradient: zero only if the lattice arrives at a_end
    at rest.
B(0) is left free: a step in the gradient at t = 0 is physical (it only sets addot(0)).
Pulses driving a below a_floor are rejected without being evolved.

The initial guess is the ZVD staircase B = Bf [1/4, 3/4, 1] of the ZVD spacing block in
tevo-daily-things.jl, sampled on the dt/2 grid. The ramp lasts exactly as long as that pulse,
T = 2 pi/w, the time ZVD needs to settle, and dt is chosen so that T is a whole number of
half-steps: then both ZVD switches, pi/w and 2 pi/w, sit on samples. The interior one is exact
(half-weight sample), but the last sample is the pinned holding gradient, so the trapezoid reads
the final jump 3/4 Bf -> Bf as a ramp over the last half-step rather than a step at T: the guess
rings with R = 4.6e-3 (7.5e-5 if that sample were 3/4 Bf). That is a physical pulse all the same,
and the optimizer removes the residual through the R term.
"""

import math

import numpy as np

from config_magspacRamp import magspacRamp
from quocs_common import (best_pulse, dcrab_algorithm_settings, finish_run, fourier_pulse,
                          halfstep_bins, jl, run_optimization)

PULSE_NAME = "magspacGrad"
TIME_NAME = "time_magspacGrad"


class magspacGrad(magspacRamp):
    """magspacRamp's lattice, endpoints and fast evolution, with the gradient as the control."""

    def __init__(self, args_dict: dict = None):
        if args_dict is None:
            args_dict = {}

        # the trapped mode the gradient drives, and the pulse grid: the ZVD pulse's own duration
        # T = 2 pi/w on 256 half-steps (see the module constants below); dt also sets the
        # trapezoid of the response integral, the evolution itself sub-steps adaptively
        self.trap_frequency = args_dict.setdefault("trap_frequency", 10.0)
        self.a_floor = args_dict.setdefault("a_floor", 0.1)
        args_dict.setdefault("ramptime", 2 * math.pi / self.trap_frequency)
        args_dict.setdefault("dt", 2 * (args_dict["ramptime"] / 256))

        super().__init__(args_dict)

    def evaluate(self, gradient_pulse) -> tuple:
        """(FoM, fidelity, ringing, min a) of one gradient pulse."""
        return tuple(jl.compute_fom_magspac_gradient([gradient_pulse], self.to_julia_dict(), self._setup))

    def get_FoM(self, pulses: list = [], parameters: list = [], timegrids: list = []) -> dict:
        fom = jl.compute_fom_magspac_gradient(pulses, self.to_julia_dict(), self._setup)[0]
        return {"FoM": float(fom)}


# must match magspacGrad's defaults above
a_start = 0.5
a_end = 2.0
trap_frequency = 10.0
ramptime = 2 * math.pi / trap_frequency     # the ZVD pulse's own duration
# 256 half-steps: a power of two, so ramptime / (dt/2) is exactly 256.0 in floating point and
# the ceil in halfstep_bins / pulse_ramp cannot add a spurious sample past T (252 half-steps
# gives 252.00000000000003). h = 0.00245, dt = 0.00491 < dt_crit(a = 0.5) = 0.00667
dt = 2 * (ramptime / 256)
B_MAX = 300.0           # hard gradient cap, 2x the ZVD pulse's peak

hold_gradient = trap_frequency ** 2 * (a_end - a_start)
assert hold_gradient <= B_MAX, f"holding a_end needs B = {hold_gradient}, above B_MAX = {B_MAX}"

zvd_guess = np.asarray(jl.zvd_gradient_pulse(a_start, a_end, trap_frequency, ramptime, dt))
assert zvd_guess.size == halfstep_bins(ramptime, dt), "ZVD guess is not on the pulse's half-step grid"
assert zvd_guess[-1] == hold_gradient, "ZVD guess does not end on the holding gradient"

pulse = fourier_pulse(pulse_name=PULSE_NAME,
                      time_name=TIME_NAME,
                      ramptime=ramptime,
                      dt=dt,
                      lower_limit=-B_MAX,
                      upper_limit=B_MAX,
                      amplitude_variation=50.0,
                      initial_guess_lambda="",
                      # pins B(T) at the guess's holding gradient, leaves B(0) free
                      scaling_lambda="lambda t: 1.0 - t / t[-1]")
pulse["initial_guess"] = {"function_type": "list_function", "list_function": zvd_guess}

optimization_dictionary = {
    "optimization_client_name": "magspacGrad_dCRAB",
    # ~10 ms per evaluation, so the budget is set by convergence, not cost
    "algorithm_settings": dcrab_algorithm_settings(super_iteration_number=10,
                                                   max_eval_total=5000),
    "pulses": [pulse],
    "parameters": [],
    "times": [{"time_name": TIME_NAME, "initial_value": ramptime}],
}


def describe(fom: magspacGrad, gradient_pulse, label: str) -> None:
    value, fidelity, ringing, min_a = fom.evaluate(gradient_pulse)
    print(f"{label}: FoM {value:.6f} = fidelity {fidelity:.6f} - ringing {ringing:.2e}; "
          f"min a {min_a:.3f}, max |B| {np.max(np.abs(gradient_pulse)):.1f} (cap {B_MAX:.0f}), "
          f"B(0) {gradient_pulse[0]:.1f}, B(T) {gradient_pulse[-1]:.1f}")


def main():
    fom = magspacGrad({"a_start": a_start, "a_end": a_end, "ramptime": ramptime, "dt": dt,
                       "trap_frequency": trap_frequency})
    describe(fom, zvd_guess, "ZVD guess")

    optimization_obj = run_optimization(optimization_dictionary, fom)

    best_controls, _ = finish_run(optimization_obj,
                                  [(PULSE_NAME, TIME_NAME, "Magnetic gradient B(t)")],
                                  title=f"Optimized gradient, |B| <= {B_MAX:.0f}, FoM {{fom:.4f}}",
                                  prefix="magspacGrad")
    describe(fom, best_pulse(best_controls, PULSE_NAME, TIME_NAME)[1], "Optimized")


if __name__ == "__main__":
    main()

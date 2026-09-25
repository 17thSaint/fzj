"""QuOCS dCRAB optimization of the magnetic spacing ramp.

Optimizes the spacing a(t) of the dipole-dipole interaction profile
U_r = intstren / (r a)^3, ramped a_start -> a_end, to transport the dd groundstate
manifold at a_start into the one at a_end, using compute_fidelity_magspac_ramp /
setup_magspac_ramp (exact-diag/magspac-ramp-control-functions.jl).

The initial guess is the ZVD spacing ramp: the closed-form response of the trapped mode
xddot + w^2 x = B(t) to the three-step ZVD gradient staircase, the same pulse as the ZVD
spacing block of tevo-daily-things.jl (fidelity 0.9111 there, vs 0.8963 for a sudden quench
and 0.9253 for a linear ramp over the same T = 1).

The update is scaled by 16 (t/T)^2 (1 - t/T)^2, which pins both the value AND the slope of
a(t) at the two ends to the guess's. The ZVD guess leaves a_start and arrives at a_end at
rest, so every pulse in the search does too: a(t) stays one that a finite gradient can
produce (a nonzero slope at an end would need a delta kick in B), and B(t) can be read back
with get_magnetic_gradient_from_spacing. Hard limits: a >= 0.1 (the evolution is
unconditionally stable there, see the Julia file) and effectively no upper bound.
"""

import numpy as np

from quocs_common import (JuliaFoM, dcrab_algorithm_settings, finish_run, fourier_pulse,
                          halfstep_bins, include_julia, jl, run_optimization)

include_julia("magspac-ramp-control-functions.jl")

PULSE_NAME = "magspacRamp"
TIME_NAME = "time_magspacRamp"


class magspacRamp(JuliaFoM):

    def __init__(self, args_dict: dict = None):
        if args_dict is None:
            args_dict = {}

        # lattice / interaction parameters
        self.Lx = args_dict.setdefault("Lx", 4)
        self.Ly = args_dict.setdefault("Ly", 4)
        self.N = args_dict.setdefault("N", 2)
        self.intstren = args_dict.setdefault("intstren", 10.0)

        # ramp endpoints; the duration is fixed (only the shape is optimized) and must match
        # the ramptime/dt the pulse dictionary's bins are derived from
        self.a_start = args_dict.setdefault("a_start", 0.5)
        self.a_end = args_dict.setdefault("a_end", 2.0)
        self.ramptime = args_dict.setdefault("ramptime", 1.0)

        # figure of merit / time evolution settings; dt only sets the half-step pulse grid
        # (dt/2), the evolution itself sub-steps adaptively
        self.speccount = args_dict.setdefault("speccount", 2)
        self.dt = args_dict.setdefault("dt", 0.005)

        # endpoint manifolds and the H = H0 + D(a) decomposition, computed once
        self._setup = jl.setup_magspac_ramp(self.to_julia_dict())

    def get_FoM(self, pulses: list = [], parameters: list = [], timegrids: list = []) -> dict:
        fidelity = jl.compute_fidelity_magspac_ramp(pulses, self.to_julia_dict(), self._setup)
        return {"FoM": fidelity}


# must match magspacRamp's defaults above
a_start = 0.5
a_end = 2.0
ramptime = 1.0
dt = 0.005
trap_frequency = 10.0   # of the ZVD guess only; the FoM never sees the gradient
a_min = 0.1
a_max = 1.0e3           # "no upper bound": far above anything the search reaches

zvd_guess = np.asarray(jl.zvd_spacing_pulse(a_start, a_end, trap_frequency, ramptime, dt))
assert zvd_guess.size == halfstep_bins(ramptime, dt), "ZVD guess is not on the pulse's half-step grid"

pulse = fourier_pulse(pulse_name=PULSE_NAME,
                      time_name=TIME_NAME,
                      ramptime=ramptime,
                      dt=dt,
                      lower_limit=a_min,
                      upper_limit=a_max,
                      amplitude_variation=0.5,
                      initial_guess_lambda="",
                      scaling_lambda="lambda t: 16.0 * (t / t[-1])**2 * (1.0 - t / t[-1])**2")
pulse["initial_guess"] = {"function_type": "list_function", "list_function": zvd_guess}

optimization_dictionary = {
    "optimization_client_name": "magspacRamp_dCRAB",
    # ~10 ms per evaluation, so the budget is set by convergence, not cost
    "algorithm_settings": dcrab_algorithm_settings(super_iteration_number=10,
                                                   max_eval_total=5000),
    "pulses": [pulse],
    "parameters": [],
    "times": [{"time_name": TIME_NAME, "initial_value": ramptime}],
}


def main():
    fom = magspacRamp({"a_start": a_start, "a_end": a_end, "ramptime": ramptime, "dt": dt})
    guess_fidelity = fom.get_FoM([zvd_guess])["FoM"]
    print(f"ZVD initial guess fidelity: {guess_fidelity:.6f}")

    optimization_obj = run_optimization(optimization_dictionary, fom)

    _, fidelity = finish_run(optimization_obj,
                             [(PULSE_NAME, TIME_NAME, "Magnetic spacing a(t)")],
                             title=f"Optimized spacing ramp, fidelity {{fom:.4f}} (ZVD guess {guess_fidelity:.4f})",
                             prefix="magspacRamp")
    print(f"Optimized fidelity: {fidelity:.6f} (ZVD guess {guess_fidelity:.6f})")


if __name__ == "__main__":
    main()

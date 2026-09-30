"""QuOCS regular-CRAB optimization of the magnetic gradient B(t), Nelder-Mead inner search.

Same problem, FoM and hard cap as config_magspacGrad.py (4x4 N=2 dd, a 0.5 -> 2.0 over
T = 2 pi/w, FoM = fidelity - ringing R, |B| <= 300); only the optimizer setup differs:

  - regular CRAB rather than dCRAB: a single super-iteration, i.e. one random draw of the
    frequencies followed by one Nelder-Mead search over their coefficients;
  - 4 basis frequencies, drawn as randomized principal harmonics: quocslib's Uniform
    distribution is stratified, so limits [0.5, 4.5] put exactly one frequency in each
    [k - 1/2, k + 1/2] cycles per T, k = 1..4;
  - both ends pinned at the gradients that hold the lattice at rest: B(0) = 0 (a_start is the
    trap's equilibrium) and B(T) = w^2 (a_end - a_start), the holding gradient. The guess
    interpolates smoothly between them, B = B_hold (1 - cos(pi t/T))/2 (zero slope at both
    ends), and the default parabolic envelope (t/T)(1 - t/T) freezes the update at both ends.
"""

import sys

import numpy as np

from config_magspacGrad import (B_MAX, PULSE_NAME, TIME_NAME, a_end, a_start, describe, dt,
                                hold_gradient, magspacGrad, ramptime, trap_frequency)
from quocs_common import (best_pulse, dcrab_algorithm_settings, finish_run, fourier_pulse,
                          halfstep_bins, run_optimization)

BASIS_FREQUENCIES = 4

times = np.linspace(0.0, ramptime, halfstep_bins(ramptime, dt))
smooth_guess = hold_gradient * (1.0 - np.cos(np.pi * times / ramptime)) / 2
assert smooth_guess[0] == 0.0 and smooth_guess[-1] == hold_gradient, "guess does not hit the holding gradients"

pulse = fourier_pulse(pulse_name=PULSE_NAME,
                      time_name=TIME_NAME,
                      ramptime=ramptime,
                      dt=dt,
                      lower_limit=-B_MAX,
                      upper_limit=B_MAX,
                      amplitude_variation=50.0,
                      initial_guess_lambda="",
                      basis_vector_number=BASIS_FREQUENCIES)
pulse["initial_guess"] = {"function_type": "list_function", "list_function": smooth_guess}
# principal harmonics k = 1..4 cycles per T, each randomized by up to +-1/2
pulse["basis"]["random_super_parameter_distribution"].update(
    {"lower_limit": 0.5, "upper_limit": BASIS_FREQUENCIES + 0.5})

optimization_dictionary = {
    "optimization_client_name": "magspacGrad_CRAB",
    # one super-iteration = regular CRAB; ~10 ms per evaluation
    "algorithm_settings": dcrab_algorithm_settings(super_iteration_number=1,
                                                   max_eval_total=5000),
    "pulses": [pulse],
    "parameters": [],
    "times": [{"time_name": TIME_NAME, "initial_value": ramptime}],
}


def main(seed: int = None, super_iterations: int = 1):
    # super_iterations > 1 turns this into dCRAB from the same guess and pins: each super-iteration
    # draws 4 new frequencies (same band) and gets its own Nelder-Mead search, capped at 400 evals
    if super_iterations > 1:
        settings = optimization_dictionary["algorithm_settings"]
        settings["super_iteration_number"] = super_iterations
        settings["dsm_settings"]["stopping_criteria"]["max_eval"] = 400
        optimization_dictionary["optimization_client_name"] = "magspacGradSmooth_dCRAB"
    # plain CRAB draws the frequencies once, so the result depends on the draw: a seed fixes it
    # and is recorded in the run's folder name
    if seed is not None:
        optimization_dictionary["algorithm_settings"]["random_number_generator"] = {"seed_number": seed}
        optimization_dictionary["optimization_client_name"] += f"_seed{seed}"

    fom = magspacGrad({"a_start": a_start, "a_end": a_end, "ramptime": ramptime, "dt": dt,
                       "trap_frequency": trap_frequency})
    describe(fom, smooth_guess, "Smooth guess")

    optimization_obj = run_optimization(optimization_dictionary, fom)

    best_controls, _ = finish_run(optimization_obj,
                                  [(PULSE_NAME, TIME_NAME, "Magnetic gradient B(t)")],
                                  title=f"CRAB gradient, {super_iterations} SI x {BASIS_FREQUENCIES} freqs, seed {seed}, |B| <= {B_MAX:.0f}, FoM {{fom:.4f}}",
                                  prefix="magspacGradCRAB")
    describe(fom, best_pulse(best_controls, PULSE_NAME, TIME_NAME)[1], "Optimized")


if __name__ == "__main__":
    # python config_magspacGradCRAB.py [seed] [super_iterations]
    main(int(sys.argv[1]) if len(sys.argv) > 1 else None,
         int(sys.argv[2]) if len(sys.argv) > 2 else 1)

"""QuOCS dCRAB optimization of the interaction strength ramp.

Optimizes the pulse shape of an interaction_strength ramp connecting the
strongly interacting ULR ground-state manifold to the FCI ground-state
manifold, using compute_fidelity_intstren_ramp / setup_intstren_ramp
(exact-diag/intstren-ramp-control-functions.jl).

config_intstrenRamp_AD.py is the gradient-based (automatic differentiation)
counterpart of this optimization.
"""

from quocs_common import (JuliaFoM, dcrab_algorithm_settings, fourier_pulse,
                          include_julia, jl, linear_ramp_lambda, load_best_controls,
                          plot_best_pulses, run_optimization)

include_julia("intstren-ramp-control-functions.jl")


class intstrenRamp(JuliaFoM):

    def __init__(self, args_dict: dict = None):
        if args_dict is None:
            args_dict = {}

        # hamiltonian / lattice parameters
        self.Lx = args_dict.setdefault("Lx", 4)
        self.Ly = args_dict.setdefault("Ly", 4)
        self.N = args_dict.setdefault("N", 2)
        self.lr = args_dict.setdefault("lr", "all")
        self.if_periodic_x = args_dict.setdefault("if_periodic_x", True)
        self.if_periodic_y = args_dict.setdefault("if_periodic_y", True)

        # interaction profile over distance: "flat" (default), "exp", "gaussian",
        # "rydberg" or "dd", with their respective length/width parameters (defaults
        # mirror get_normal_model_params_ed); forwarded to Julia by
        # setup_intstren_ramp's scaling loop
        self.scaling_type = args_dict.setdefault("scaling_type", "flat")
        self.corr_length = args_dict.setdefault("corr_length", float(self.Ly))
        self.sigma = args_dict.setdefault("sigma", 1.0)
        self.blockade_radius = args_dict.setdefault("blockade_radius", 1.0)
        self.magnetic_spacing = args_dict.setdefault("magnetic_spacing", 1.0)

        # ramp targets. The ramp duration is fixed (only the pulse shape is optimized)
        # and must match the ramptime value the pulse dictionary's bins are derived from
        self.intstren_start = args_dict.setdefault("intstren_start", 10.0)
        self.intstren_end = args_dict.setdefault("intstren_end", 0.0)
        self.ramptime = args_dict.setdefault("ramptime", 1.0)

        # figure of merit / time evolution settings
        self.speccount = args_dict.setdefault("speccount", 2)
        self.dt = args_dict.setdefault("dt", 0.005)

        # endpoint manifolds from ED, computed once and reused every evaluation
        self._setup = jl.setup_intstren_ramp(self.to_julia_dict())

    def get_FoM(self, pulses: list = [], parameters: list = [], timegrids: list = []) -> dict:
        fidelity = jl.compute_fidelity_intstren_ramp(pulses, self.to_julia_dict(), self._setup)
        return {"FoM": fidelity}


# must match intstrenRamp's intstren_start / intstren_end / ramptime / dt defaults above
intstren_start = 10.0
intstren_end = 0.0
intstren_max = 12.0
ramptime = 1.0
dt = 0.005

optimization_dictionary = {
    "optimization_client_name": "intstrenRamp_dCRAB",
    "algorithm_settings": dcrab_algorithm_settings(super_iteration_number=5,
                                                   max_eval_total=200),
    "pulses": [
        fourier_pulse(pulse_name="intstrenRamp",
                      time_name="time_intstrenRamp",
                      ramptime=ramptime,
                      dt=dt,
                      lower_limit=0.0,
                      upper_limit=intstren_max,
                      amplitude_variation=3.0,
                      initial_guess_lambda=linear_ramp_lambda(intstren_start, intstren_end)),
    ],
    "parameters": [],
    "times": [{"time_name": "time_intstrenRamp", "initial_value": ramptime}],
}


def main():
    optimization_obj = run_optimization(optimization_dictionary, intstrenRamp({
        "intstren_start": intstren_start,
        "intstren_end": intstren_end,
        "ramptime": ramptime,
        "dt": dt,
    }))

    best_controls = load_best_controls(optimization_obj)
    fidelity = abs(optimization_obj.opt_alg_obj.best_FoM)
    plot_best_pulses(best_controls,
                     [("intstrenRamp", "time_intstrenRamp", "Interaction strength")],
                     title=f"Optimized Intstren-Ramp Fidelity: {fidelity:.4f}",
                     filename="intstrenRamp_optimized.png")


if __name__ == "__main__":
    main()

"""QuOCS dCRAB optimization of the two-stage pinned-state preparation ramp.

Optimizes the pulse shapes of two sequential hopping ramps (ty first, then tx)
that connect a real-space corner-pinned product state to the isotropic-hopping
FCI ground-state manifold, using compute_fidelity_pinned_ramp /
setup_pinned_ramp (exact-diag/pinned-ramp-control-functions.jl).
"""

from quocs_common import (JuliaFoM, dcrab_algorithm_settings, fourier_pulse,
                          include_julia, jl, linear_ramp_lambda, load_best_controls,
                          plot_best_pulses, run_optimization)

include_julia("pinned-ramp-control-functions.jl")


class pinnedRamp(JuliaFoM):

    def __init__(self, args_dict: dict = None):
        if args_dict is None:
            args_dict = {}

        # hamiltonian / lattice parameters
        self.Lx = args_dict.setdefault("Lx", 4)
        self.Ly = args_dict.setdefault("Ly", 4)
        self.N = args_dict.setdefault("N", 2)
        self.interaction_strength = args_dict.setdefault("interaction_strength", 0.0)
        self.lr = args_dict.setdefault("lr", "all")
        self.if_periodic_x = args_dict.setdefault("if_periodic_x", True)
        self.if_periodic_y = args_dict.setdefault("if_periodic_y", True)

        # pinned starting configuration: flat (col,row) pairs, one per particle
        self.starting_config = args_dict.setdefault("starting_config", [1, 1, 1, 2])

        # ramp targets. Ramp durations are fixed (only pulse shape is optimized) and
        # must match the ramptime values the pulse dictionaries' bins are derived from
        self.end_tx = args_dict.setdefault("end_tx", 1.0)
        self.end_ty = args_dict.setdefault("end_ty", 1.0)
        self.ramptime_firstramp = args_dict.setdefault("ramptime_firstramp", 0.5)
        self.ramptime_secondramp = args_dict.setdefault("ramptime_secondramp", 0.5)

        # figure of merit / time evolution settings
        self.speccount = args_dict.setdefault("speccount", 2)
        self.dt = args_dict.setdefault("dt", 0.005)

        # pinned starting state and target manifold, computed once and reused
        self._setup = jl.setup_pinned_ramp(self.to_julia_dict())

    def get_FoM(self, pulses: list = [], parameters: list = [], timegrids: list = []) -> dict:
        fidelity = jl.compute_fidelity_pinned_ramp(pulses, self.to_julia_dict(), self._setup)
        return {"FoM": fidelity}


# must match pinnedRamp's ramptime_firstramp / ramptime_secondramp / dt defaults above
ramptime_firstramp = 0.5
ramptime_secondramp = 0.5
dt = 0.005

optimization_dictionary = {
    "optimization_client_name": "pinnedRamp_dCRAB",
    "algorithm_settings": dcrab_algorithm_settings(super_iteration_number=5,
                                                   max_eval_total=200),
    "pulses": [
        fourier_pulse(pulse_name="tyRamp",
                      time_name="time_tyRamp",
                      ramptime=ramptime_firstramp,
                      dt=dt,
                      lower_limit=0.0,
                      upper_limit=1.2,
                      amplitude_variation=0.3,
                      initial_guess_lambda=linear_ramp_lambda(0.0, 1.0)),
        fourier_pulse(pulse_name="txRamp",
                      time_name="time_txRamp",
                      ramptime=ramptime_secondramp,
                      dt=dt,
                      lower_limit=0.0,
                      upper_limit=1.2,
                      amplitude_variation=0.3,
                      initial_guess_lambda=linear_ramp_lambda(0.0, 1.0)),
    ],
    "parameters": [],
    "times": [
        {"time_name": "time_tyRamp", "initial_value": ramptime_firstramp},
        {"time_name": "time_txRamp", "initial_value": ramptime_secondramp},
    ],
}


def main():
    optimization_obj = run_optimization(optimization_dictionary, pinnedRamp({
        "ramptime_firstramp": ramptime_firstramp,
        "ramptime_secondramp": ramptime_secondramp,
        "dt": dt,
    }))

    best_controls = load_best_controls(optimization_obj)
    fidelity = abs(optimization_obj.opt_alg_obj.best_FoM)
    plot_best_pulses(best_controls,
                     [("tyRamp", "time_tyRamp", "ty"),
                      ("txRamp", "time_txRamp", "tx")],
                     title=f"Optimized Pinned-Ramp Fidelity: {fidelity:.4f}",
                     filename="pinnedRamp_optimized.png")


if __name__ == "__main__":
    main()

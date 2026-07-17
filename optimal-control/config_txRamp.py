"""QuOCS dCRAB optimization of the tx hopping ramp.

Optimizes the pulse shape of a tx ramp (anisotropic -> isotropic hopping) that
connects the tx ~ 0 ground state to the isotropic ground state, using
compute_fidelity / setup_tx_ramp (exact-diag/control-functions.jl).
"""

from quocs_common import (JuliaFoM, dcrab_algorithm_settings, fourier_pulse, jl,
                          linear_ramp_lambda, load_best_controls, plot_best_pulses,
                          run_optimization)


class txRamp(JuliaFoM):

    def __init__(self, args_dict: dict = None):
        if args_dict is None:
            args_dict = {}

        # control parameter
        self.tx_initial = args_dict.setdefault("tx_initial", 0.01)
        self.tx_final = args_dict.setdefault("tx_final", 1.0)

        # hamiltonian parameters
        self.ty = args_dict.setdefault("ty", 1.0)
        self.if_periodic_x = args_dict.setdefault("if_periodic_x", True)
        self.if_periodic_y = args_dict.setdefault("if_periodic_y", True)
        self.interaction_strength = args_dict.setdefault("interaction_strength", 0.0)
        self.lr = args_dict.setdefault("lr", "all")

        # lattice parameters
        self.Lx = args_dict.setdefault("Lx", 4)
        self.Ly = args_dict.setdefault("Ly", 4)
        self.N = args_dict.setdefault("N", 2)
        self.filling = args_dict.setdefault("filling", 0.5)

        # ramp duration is fixed (only the pulse shape is optimized) and must match
        # the ramptime value the pulse dictionary's bins are derived from
        self.ramptime = args_dict.setdefault("ramptime", 2.0)
        self.dt = args_dict.setdefault("dt", 0.05)

        # other parameters
        self.if_reading = args_dict.setdefault("if_reading", False)
        self.nev = args_dict.setdefault("nev", 10)

        # running parameters
        self.if_find_data = args_dict.setdefault("if_find_data", False)
        self.if_save_data = args_dict.setdefault("if_save_data", False)

        # endpoint ground states from ED, computed once and reused every evaluation
        self._setup = jl.setup_tx_ramp(self.to_julia_dict())

    def get_FoM(self, pulses: list = [], parameters: list = [], timegrids: list = []) -> dict:
        fidelity = jl.compute_fidelity(pulses, self.to_julia_dict(), self._setup)
        return {"FoM": fidelity}


# must match txRamp's ramptime / dt defaults above
ramptime = 2.0
dt = 0.05

optimization_dictionary = {
    "optimization_client_name": "txRamp_dCRAB",
    "algorithm_settings": dcrab_algorithm_settings(super_iteration_number=5,
                                                   max_eval_total=200),
    "pulses": [
        fourier_pulse(pulse_name="txRamp",
                      time_name="time_txRamp",
                      ramptime=ramptime,
                      dt=dt,
                      lower_limit=0.0,
                      upper_limit=2.0,
                      amplitude_variation=0.3,
                      initial_guess_lambda=linear_ramp_lambda(0.001, 1.0)),
    ],
    "parameters": [],
    "times": [{"time_name": "time_txRamp", "initial_value": ramptime}],
}


def main():
    optimization_obj = run_optimization(optimization_dictionary, txRamp({
        "ramptime": ramptime,
        "dt": dt,
    }))

    best_controls = load_best_controls(optimization_obj)
    fidelity = abs(optimization_obj.opt_alg_obj.best_FoM)
    plot_best_pulses(best_controls,
                     [("txRamp", "time_txRamp", "tx")],
                     title=f"Optimized Ramp Fidelity: {fidelity:.4f}",
                     filename="txRamp_optimized.png")


if __name__ == "__main__":
    main()

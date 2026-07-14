
import math
import matplotlib.pyplot as plt
import os
from quocslib.utils.AbstractFoM import AbstractFoM
import functools

# Set up the Julia environment from the exact-diag project
# Julia version is limited to 1.11
juliaup_bin = "/home/patrick/.juliaup/bin"
os.environ["PATH"] = f"{juliaup_bin}:{os.environ.get('PATH', '')}"
os.environ["JULIAUP_CHANNEL"] = "1.11"
# juliacall does not read JULIA_PROJECT -- it needs PYTHON_JULIACALL_PROJECT (paired
# with PYTHON_JULIACALL_EXE) or it silently falls back to its own private juliapkg-managed
# environment, which doesn't have exact-diag's dependencies (JLD2, ITensors, ...)
os.environ["PYTHON_JULIACALL_PROJECT"] = os.path.abspath("../exact-diag")
os.environ["PYTHON_JULIACALL_EXE"] = os.path.join(juliaup_bin, "julia")
# CONFIG['opt_handle_signals'] is None triggers a juliacall init-time NameError
# (references an undefined 'Base') on this Julia/juliacall version combination -- set
# explicitly to skip that code path (see https://juliapy.github.io/PythonCall.jl/stable/faq)
os.environ["PYTHON_JULIACALL_HANDLE_SIGNALS"] = "yes"

from juliacall import Main as jl

jl.include("../exact-diag/control-functions.jl")
jl.include("../exact-diag/time-evolution.jl")
jl.include("../exact-diag/intstren-ramp-control-functions.jl")


def python_dict_to_julia_dict(data: dict):
    pairs = [jl.Pair(str(k), v) for k, v in data.items()]
    return jl.Dict(pairs)


class intstrenRamp(AbstractFoM):

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

        # ramp targets. The ramp duration is fixed (only the pulse shape is optimized)
        # and must match the ramptime value used below to build the pulse dictionary
        self.intstren_start = args_dict.setdefault("intstren_start", 10.0)
        self.intstren_end = args_dict.setdefault("intstren_end", 0.0)
        self.ramptime = args_dict.setdefault("ramptime", 1.0)

        # figure of merit / time evolution settings
        self.speccount = args_dict.setdefault("speccount", 2)
        self.dt = args_dict.setdefault("dt", 0.005)

    def to_julia_dict(self):
        payload = {
            k: v
            for k, v in vars(self).items()
            if not k.startswith("_") and not callable(v)
        }
        return python_dict_to_julia_dict(payload)

    def get_FoM(self, pulses: list = [], parameters: list = [], timegrids: list = []) -> dict:

        fidelity = jl.compute_fidelity_intstren_ramp(pulses, self.to_julia_dict())

        return {"FoM": fidelity}


optimization_dictionary = {"optimization_client_name": "intstrenRamp_dCRAB"}

optimization_dictionary["algorithm_settings"] = {"algorithm_name": "dCRAB"}

optimization_dictionary["algorithm_settings"]["optimization_direction"] = "maximization"
optimization_dictionary["algorithm_settings"]["super_iteration_number"] = 5
optimization_dictionary["algorithm_settings"]["max_eval_total"] = 200

# will need expert advice on choices made here
dsm_settings = {
        "general_settings": {
            "dsm_algorithm_name": "NelderMead",
            "is_adaptive": False
        },
        "stopping_criteria": {
            "xatol": 1e-4,
            "fatol": 1e-6,
            "change_based_stop": {
                "cbs_funct_evals": 200,
                "cbs_change": 0.01
            }
        }
    }

optimization_dictionary["algorithm_settings"]["dsm_settings"] = dsm_settings

# these three must match intstrenRamp's intstren_start / ramptime / dt defaults above
intstren_start = 10.0
intstren_end = 0.0
ramptime = 1.0
dt = 0.005

# pulse_ramp (exact-diag/time-evolution.jl) samples the pulse on the RK4
# half-step grid (spacing dt/2) up to ending_time, so the pulse array QuOCS optimizes
# must have exactly ceil(2 * ramptime / dt) + 1 bins
bins_intstren = math.ceil(2 * ramptime / dt) + 1

pulse_intstren = {"pulse_name": "intstrenRamp",
           "upper_limit": 12.0,
           "lower_limit": 0.0,
           "bins_number": bins_intstren,
           "amplitude_variation": 3.0,
           "time_name": "time_intstrenRamp",
           "shaping_options": [
               "add_base_pulse",
               "add_new_update_pulse",
               "scale_pulse",
               "add_initial_guess",
               "limit_pulse"
           ]
           }

pulse_intstren["initial_guess"] = {
    "function_type": "lambda_function",
    "lambda_function": "lambda t: 10.0 + (0.0 - 10.0) * (t / t[-1])"
    }

pulse_intstren["scaling_function"] = {
    "function_type": "lambda_function",
    "lambda_function": "lambda t: (t / t[-1]) * (1.0 - t / t[-1])"
    }

pulse_intstren["basis"] = {
                "basis_name": "Fourier",
                "basis_vector_number": 5,
                "random_super_parameter_distribution": {
                    "distribution_name": "Uniform",
                    "lower_limit": 0.01,
                    "upper_limit": 10.0
                }
            }


time_intstrenRamp = {"time_name": "time_intstrenRamp",
                "initial_value": ramptime}


optimization_dictionary["pulses"] = [pulse_intstren]
optimization_dictionary["parameters"] = []
optimization_dictionary["times"] = [time_intstrenRamp]


from quocslib.Optimizer import Optimizer
import time

optimization_obj = Optimizer(optimization_dictionary, intstrenRamp({
    "is_maximization": True,
    "intstren_start": intstren_start,
    "intstren_end": intstren_end,
    "ramptime": ramptime,
    "dt": dt,
}))


time1 = time.time()
optimization_obj.execute()
time2 = time.time()
print("The optimization took {seconds} seconds".format(seconds=time2 - time1))


opt_alg_obj = optimization_obj.opt_alg_obj
controls = opt_alg_obj.get_best_controls()

pulse_intstren_opt, timegrid_intstren_opt = controls["pulses"][0], controls["timegrids"][0]

fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.plot(timegrid_intstren_opt, pulse_intstren_opt)
ax.set_xlabel("Time")
ax.set_ylabel("Interaction strength")
ax.grid()
fig.suptitle("Optimized Intstren-Ramp Fidelity: {fidelity:.4f}".format(fidelity=opt_alg_obj.best_FoM))
plt.savefig("local-figs/intstrenRamp_optimized.png")




"fin"


import matplotlib.pyplot as plt
import os
from quocslib.utils.AbstractFoM import AbstractFoM
from quocslib.timeevolution.piecewise_integrator import pw_evolution
import functools

# Set up the Julia environment from the exact-diag project
# Julia version is limited to 1.11 
juliaup_bin = "/home/patrick/.juliaup/bin"
os.environ["PATH"] = f"{juliaup_bin}:{os.environ.get('PATH', '')}"
os.environ["JULIAUP_CHANNEL"] = "1.11"
os.environ["JULIA_PROJECT"] = os.path.abspath("../exact-diag")

from juliacall import Main as jl

jl.include("../exact-diag/control-functions.jl")
jl.include("../exact-diag/time-evolution.jl")


def python_dict_to_julia_dict(data: dict):
    pairs = [jl.Pair(str(k), v) for k, v in data.items()]
    return jl.Dict(pairs)


class txRamp(AbstractFoM):

    def __init__(self, args_dict: dict = None):
        if args_dict is None:
            args_dict = {}    

        #Dict([("output_level",1),("periodic_potential_strength",ppstren),("tx",anis),("ty",1.0),("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10)])
        
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

        # other parameters
        self.if_reading = args_dict.setdefault("if_reading", False)
        self.nev = args_dict.setdefault("nev", 10)

        # running parameters
        self.if_find_data = args_dict.setdefault("if_find_data", False)
        self.if_save_data = args_dict.setdefault("if_save_data", False)


    def to_julia_dict(self):
        payload = {
            k: v
            for k, v in vars(self).items()
            if not k.startswith("_") and not callable(v)
        }
        return python_dict_to_julia_dict(payload)


    def get_FoM(self, pulses: list = [], parameters: list = [], timegrids: list = []) -> dict:
        
        fidelity = jl.compute_fidelity(pulses, self.to_julia_dict())

        return {"FoM": fidelity}




optimization_dictionary = {"optimization_client_name": "txRamp_dCRAB"}

optimization_dictionary["algorithm_settings"] = { "algorithm_name": "dCRAB"}

optimization_dictionary["algorithm_settings"]["optimization_direction"] = "maximization"
optimization_dictionary["algorithm_settings"]["super_iteration_number"] = 3
optimization_dictionary["algorithm_settings"]["max_eval_total"] = 20

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

pulse_tx = {"pulse_name": "txRamp",
           "upper_limit": 2.0,
           "lower_limit": 0.0,
           "bins_number": 79,
           "amplitude_variation": 0.3,
           "time_name": "time_txRamp",
           "shaping_options": [
               "add_base_pulse",
               "add_new_update_pulse",
               "scale_pulse",
               "add_initial_guess",
               "limit_pulse"
           ]
           }



pulse_tx["initial_guess"] = {
    "function_type": "lambda_function",
    "lambda_function": "lambda t: 0.01 + (1.0 - 0.01) * (t / t[-1])"
    }


pulse_tx["scaling_function"] = {
    "function_type": "lambda_function",
    "lambda_function": "lambda t: (t / t[-1]) * (1.0 - t / t[-1])"
    }



pulse_tx["basis"] = {
                "basis_name": "Fourier",
                "basis_vector_number": 5,
                "random_super_parameter_distribution": {
                    "distribution_name": "Uniform",
                    "lower_limit": 0.01,
                    "upper_limit": 10.0
                }
            }


time_txRamp = {"time_name": "time_txRamp",
                "initial_value": 1.0 }


optimization_dictionary["pulses"] = [pulse_tx]
optimization_dictionary["parameters"] = []
optimization_dictionary["times"] = [time_txRamp]


from quocslib.Optimizer import Optimizer
import time

optimization_obj = Optimizer(optimization_dictionary, txRamp({"is_maximization": True}))


time1 = time.time()
optimization_obj.execute()
time2 = time.time()
print("The optimization took {seconds} seconds".format(seconds=time2 - time1))


opt_alg_obj = optimization_obj.opt_alg_obj
controls = opt_alg_obj.get_best_controls()

pulse,timegrid = controls["pulses"][0], controls["timegrids"][0]

fig = plt.figure()
plt.plot(timegrid, pulse)
plt.xlabel("Time")
plt.ylabel("tx")
plt.title("Optimized Ramp Fidelity: {fidelity:.4f}".format(fidelity=opt_alg_obj.best_FoM))
plt.grid()
plt.savefig("local-figs/txRamp_optimized.png")










"fin"
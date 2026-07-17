"""Shared infrastructure for the QuOCS optimal-control configs.

Importing this module boots the Julia environment (juliacall) from the exact-diag
project and includes the base control/time-evolution files; the configured `jl`
handle is exported for the configs to call their figure-of-merit functions with.

Also collects everything the config files would otherwise copy-paste: the Julia
dict conversion, the pulse/bins conventions of pulse_ramp (time-evolution.jl),
the standard dCRAB settings block, the timed optimizer runner, and reading /
plotting of the best controls from the BestDump npz.
"""

import glob
import math
import os
import time

_HERE = os.path.dirname(os.path.abspath(__file__))
EXACT_DIAG_DIR = os.path.normpath(os.path.join(_HERE, "..", "exact-diag"))
FIGURES_DIR = os.path.join(_HERE, "local-figs")

# Set up the Julia environment from the exact-diag project before importing juliacall.
# Julia version is limited to 1.11
_JULIAUP_BIN = "/home/patrick/.juliaup/bin"
os.environ["PATH"] = f"{_JULIAUP_BIN}:{os.environ.get('PATH', '')}"
os.environ["JULIAUP_CHANNEL"] = "1.11"
# juliacall does not read JULIA_PROJECT -- it needs PYTHON_JULIACALL_PROJECT (paired
# with PYTHON_JULIACALL_EXE) or it silently falls back to its own private juliapkg-managed
# environment, which doesn't have exact-diag's dependencies (JLD2, ITensors, ...)
os.environ["PYTHON_JULIACALL_PROJECT"] = EXACT_DIAG_DIR
os.environ["PYTHON_JULIACALL_EXE"] = os.path.join(_JULIAUP_BIN, "julia")
# CONFIG['opt_handle_signals'] is None triggers a juliacall init-time NameError
# (references an undefined 'Base') on this Julia/juliacall version combination -- set
# explicitly to skip that code path (see https://juliapy.github.io/PythonCall.jl/stable/faq)
os.environ["PYTHON_JULIACALL_HANDLE_SIGNALS"] = "yes"

from juliacall import Main as jl  # noqa: E402  (env vars above must precede this)

from quocslib.utils.AbstractFoM import AbstractFoM  # noqa: E402


def include_julia(filename: str) -> None:
    """Include an exact-diag Julia file by name, e.g. "control-functions.jl"."""
    jl.include(os.path.join(EXACT_DIAG_DIR, filename))


# base files every figure of merit depends on; the exact-diag include chain (pwd-based
# include_other_files) still requires running from inside the fzj tree
include_julia("control-functions.jl")
include_julia("time-evolution.jl")


def python_dict_to_julia_dict(data: dict):
    pairs = [jl.Pair(str(k), v) for k, v in data.items()]
    return jl.Dict(pairs)


class JuliaFoM(AbstractFoM):
    """Base class for figures of merit backed by exact-diag Julia functions.

    Subclasses set their parameters as public attributes in __init__;
    to_julia_dict() serializes exactly those (private "_" attributes such as
    cached Julia setup objects are excluded).
    """

    def to_julia_dict(self):
        payload = {
            k: v
            for k, v in vars(self).items()
            if not k.startswith("_") and not callable(v)
        }
        return python_dict_to_julia_dict(payload)


def halfstep_bins(ramptime: float, dt: float) -> int:
    """Number of pulse samples pulse_ramp (time-evolution.jl) expects.

    The pulse is sampled on the RK4 half-step grid (spacing dt/2) up to
    ending_time, so the array QuOCS optimizes must have exactly
    ceil(2 * ramptime / dt) + 1 bins.
    """
    return math.ceil(2 * ramptime / dt) + 1


def fourier_pulse(pulse_name: str,
                  time_name: str,
                  ramptime: float,
                  dt: float,
                  lower_limit: float,
                  upper_limit: float,
                  amplitude_variation: float,
                  initial_guess_lambda: str,
                  basis_vector_number: int = 5) -> dict:
    """Standard dCRAB/AD pulse dictionary: Fourier basis, parabolic scaling that
    pins both endpoints of the update, and hard amplitude limits."""
    return {
        "pulse_name": pulse_name,
        "upper_limit": upper_limit,
        "lower_limit": lower_limit,
        "bins_number": halfstep_bins(ramptime, dt),
        "amplitude_variation": amplitude_variation,
        "time_name": time_name,
        "shaping_options": [
            "add_base_pulse",
            "add_new_update_pulse",
            "scale_pulse",
            "add_initial_guess",
            "limit_pulse",
        ],
        "initial_guess": {
            "function_type": "lambda_function",
            "lambda_function": initial_guess_lambda,
        },
        "scaling_function": {
            "function_type": "lambda_function",
            "lambda_function": "lambda t: (t / t[-1]) * (1.0 - t / t[-1])",
        },
        "basis": {
            "basis_name": "Fourier",
            "basis_vector_number": basis_vector_number,
            "random_super_parameter_distribution": {
                "distribution_name": "Uniform",
                "lower_limit": 0.01,
                "upper_limit": 10.0,
            },
        },
    }


def linear_ramp_lambda(start: float, end: float) -> str:
    """Initial-guess lambda string for a linear ramp from start to end."""
    return f"lambda t: {start} + ({end} - {start}) * (t / t[-1])"


def dcrab_algorithm_settings(super_iteration_number: int = 5,
                             max_eval_total: int = 200) -> dict:
    """The standard dCRAB algorithm_settings block (NelderMead inner search)."""
    return {
        "algorithm_name": "dCRAB",
        "optimization_direction": "maximization",
        "super_iteration_number": super_iteration_number,
        "max_eval_total": max_eval_total,
        "dsm_settings": {
            "general_settings": {
                "dsm_algorithm_name": "NelderMead",
                "is_adaptive": False,
            },
            "stopping_criteria": {
                "xatol": 1e-4,
                "fatol": 1e-6,
                "change_based_stop": {
                    "cbs_funct_evals": max_eval_total,
                    "cbs_change": 0.01,
                },
            },
        },
    }


def run_optimization(optimization_dictionary: dict, fom_object):
    """Build the QuOCS Optimizer, execute it, and report the wall time."""
    from quocslib.Optimizer import Optimizer

    optimization_obj = Optimizer(optimization_dictionary, fom_object)
    start = time.time()
    optimization_obj.execute()
    print(f"The optimization took {time.time() - start:.1f} seconds")
    return optimization_obj


def load_best_controls(optimization_obj):
    """Load the best controls from the run's BestDump npz.

    Use this instead of opt_alg_obj.get_best_controls(): the latter reconstructs
    the pulse from the *current* (last super iteration's) random basis, which no
    longer matches the recorded best coefficients when a super iteration is cut
    off by the evaluation budget.
    """
    import numpy as np

    pattern = os.path.join(optimization_obj.results_path, "*_best_controls.npz")
    return np.load(glob.glob(pattern)[0])


def best_pulse(best_controls, pulse_name: str, time_name: str):
    """Extract (timegrid, pulse) for one control as real 1-D numpy arrays."""
    import numpy as np

    pulse = np.real(np.asarray(best_controls[pulse_name])).ravel()
    timegrid = np.real(np.asarray(best_controls[f"time_grid_for_{pulse_name}"])).ravel()
    return timegrid, pulse


def plot_best_pulses(best_controls, pulse_specs, title: str, filename: str) -> str:
    """Plot the optimized pulses side by side and save under local-figs.

    pulse_specs: list of (pulse_name, time_name, ylabel) tuples.
    """
    import matplotlib.pyplot as plt

    n = len(pulse_specs)
    fig, axes = plt.subplots(1, n, figsize=(6 if n == 1 else 5 * n, 4), squeeze=False)
    for ax, (pulse_name, time_name, ylabel) in zip(axes[0], pulse_specs):
        timegrid, pulse = best_pulse(best_controls, pulse_name, time_name)
        ax.plot(timegrid, pulse)
        ax.set_xlabel("Time")
        ax.set_ylabel(ylabel)
        ax.grid()
    fig.suptitle(title)
    os.makedirs(FIGURES_DIR, exist_ok=True)
    figure_path = os.path.join(FIGURES_DIR, filename)
    fig.savefig(figure_path)
    print(f"Saved {figure_path}")
    return figure_path

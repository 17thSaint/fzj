
import matplotlib.pyplot as plt
import os
import numpy as np
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

jl.include("../exact-diag/execute-ed.jl")
jl.include("../exact-diag/time-evolution.jl")


result = jl.sqrt(2.0)
print(result)























"fin"
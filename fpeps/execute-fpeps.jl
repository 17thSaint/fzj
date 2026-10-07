#####################################################
#=

This file contains set-up for fPEPS simulations

Depends on:


=#
######################################################

include("../other-funcs/include-other-files.jl")

using QuantumNaturalGradient, QuantumNaturalfPEPS, ITensors, ITensorMPS

lx,ly,n = 8,4,4
intstren = 50.0
bond_dim = 2
xi = 0.5

params_dict = Dict([("outputlevel",1),("scaling","dd"),("magnetic_spacing",xi),("lr","all"),("hopping_anisotropy",1.0),("Lx",lx),("Ly",ly),("particles",n),("if_save_data",false),("filling",0.5),("if_find_data",false),("onsite_strength",intstren),("if_periodic_phys",false),("if_periodic_synth",false)])

hilbert = siteinds("Boson",lx,ly)
peps = PEPS(hilbert; bond_dim=bond_dim)
QuantumNaturalfPEPS.multiply_algebraic_spectrum!(peps, 3.)

include_other_files(["other-funcs/basic-2d-observables.jl","synth-dims/long-range-ttn.jl"])

model_paras = get_normal_model_params(params_dict)
net = build_HH_net(model_paras)

ham_opsum = long_range_HH_ham(net,model_paras[:ts],model_paras[:alpha]; model_paras...)

# Generate Operators for QNG
Oks_and_Eks = QuantumNaturalfPEPS.generate_Oks_and_Eks(peps, ham_opsum)

# Setup the Integrator and Solver
integrator = QuantumNaturalGradient.Euler(lr=0.05)
solver = QuantumNaturalGradient.EigenSolver()

# Define a Parameters object to be evolved
θ = QuantumNaturalGradient.Parameters(peps)

# Evolve for a fixed (small) number of iterations as a demo
@time loss_value, trained_θ, misc = QuantumNaturalGradient.evolve(Oks_and_Eks, θ; 
        integrator, 
        verbosity=2,
        solver,
        sample_nr=1000,
        maxiter=10,)






















"fin"
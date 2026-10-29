#####################################################
#=

This file runs a new fPEPS groundstate search (cluster script, started from fpeps/), saving every save_every
iterations (default 10) unless if_save_data = false. Resume a saved run with continue_fpeps (savings.jl) instead

Depends on:
        other-funcs/include-other-files.jl (include_other_files)
        fpeps/load-fpeps.jl, of which
                build-ham.jl (get_fpeps_model_paras, long_range_HH_ham, particle_sector_counts)
                savings.jl (make_fpeps_filename, save_initial_fpeps, saving_fpeps_callback)
                other-funcs/data-storage-funcs.jl (get_folder_location, make_parameters_filename)
                QuantumNaturalfPEPS (PEPS, generate_Oks_and_Eks), QuantumNaturalGradient (evolve, Parameters, Euler, EigenSolver)
        Distributed (worker processes for generate_Oks_and_Eks)

=#
######################################################

using Pkg
Pkg.activate(".")

include("../other-funcs/include-other-files.jl")
include_other_files(["fpeps/load-fpeps.jl"])

using Distributed

nr_procs, nr_threads = 2, 2
addprocs(nr_procs, exeflags="-t $nr_threads")

@everywhere begin
    using QuantumNaturalfPEPS
    using LinearAlgebra
    BLAS.set_num_threads(1) # Each thread should only use one Blas thread
end
BLAS.set_num_threads(prod([nr_procs, nr_threads])) # For computing the double layers

# model
lx,ly,n = 8,7,4
intstren = 18.90527
dataloc = get_folder_location("cluster-data/synth-dims/fpeps/flat-ints")

# optimisation (bond_dim also goes in params_dict, it's part of the filename)
bond_dim = 5
eigencut, contract_cutoff, sample_cutoff = 1e-6, 1e-6, 5e-4
learning_rate, sample_nr, maxiter = 0.05, 3000, 5000

params_dict = Dict([("outputlevel",1),("dataloc",dataloc),("if_check_fluxes",false),("bond_dim",bond_dim),("lr","all"),("hopping_anisotropy",1.0),("Lx",lx),("Ly",ly),("particles",n),("if_save_data",true),("filling",0.5),("onsite_strength",intstren)])
model_paras = get_fpeps_model_paras(params_dict)

hilbert = ITensors.siteinds("S=1/2",lx,ly)
peps = PEPS(ComplexF64, hilbert; bond_dim=bond_dim, contract_cutoff=contract_cutoff, sample_cutoff=sample_cutoff)
QuantumNaturalfPEPS.multiply_algebraic_spectrum!(peps, 3.)

max_counts = particle_sector_counts(lx, ly, n)
println("Sampling in the particle number N = $n sector, max_counts = $max_counts")

ham_opsum = long_range_HH_ham(model_paras)

# Generate Operators for QNG
Oks_and_Eks = QuantumNaturalfPEPS.generate_Oks_and_Eks(peps, ham_opsum; threaded=true, multiproc=true, max_counts)

# Setup the Integrator and Solver
integrator = QuantumNaturalGradient.Euler(lr=learning_rate)
solver = QuantumNaturalGradient.EigenSolver(eigencut)

# Define a Parameters object to be evolved
θ = QuantumNaturalGradient.Parameters(peps)

# Saving: like the TTNs, the PEPS goes in wavefuncfpeps-<params>.jld2 and the light data + metadata in fpeps-<params>.jld2
filename = "fpeps-" * make_parameters_filename(make_fpeps_filename(model_paras)) * ".jld2"
metadata = merge(model_paras, Dict([("sample_nr",sample_nr),("maxiter",maxiter),("learning_rate",learning_rate),("solver",string(typeof(solver))),("integrator",string(typeof(integrator))),("niter",0),("contract_cutoff",contract_cutoff),("sample_cutoff",sample_cutoff),("eigencut",eigencut)]))
filepath = model_paras["if_save_data"] ? joinpath(dataloc, save_initial_fpeps(peps, metadata; location=dataloc, name=filename)) : nothing

# Callback: overwrite the saved PEPS and history every save_every iterations and at maxiter (filepath = nothing when not saving)
callback = saving_fpeps_callback(filepath; save_every=model_paras["save_every"])

loss_value, trained_θ, misc = QuantumNaturalGradient.evolve(Oks_and_Eks, θ;
        integrator,
        verbosity=2,
        solver,
        sample_nr,
        maxiter,
        callback,)

println("Final energy: $(real(loss_value))")


"fin"

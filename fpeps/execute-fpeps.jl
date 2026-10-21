#####################################################
#=

This file contains set-up for fPEPS simulations

Depends on:
        build-ham.jl (long_range_HH_ham, save_initial_fpeps, save_fpeps, continue_fpeps)
        other-funcs/basic-2d-stuff.jl


=#
######################################################

include("../other-funcs/include-other-files.jl")
include("build-ham.jl")
include("../other-funcs/data-storage-funcs.jl")

using QuantumNaturalGradient, QuantumNaturalfPEPS, ITensors, ITensorMPS, JLD2#, Distributed

#=nr_procs, nr_threads = 2, 8
addprocs(nr_procs, exeflags="-t $nr_threads")

@everywhere begin
    using QuantumNaturalfPEPS
    using LinearAlgebra
    BLAS.set_num_threads(1) # Each thread should only use one Blas thread
end
BLAS.set_num_threads(16) # For computing the double layers
=#

lx,ly,n = 4,4,2
intstren = 50.0
bond_dim = 6

params_dict = Dict([("outputlevel",1),("if_check_fluxes",false),("lr","all"),("hopping_anisotropy",1.0),("Lx",lx),("Ly",ly),("particles",n),("if_save_data",false),("filling",0.5),("if_find_data",false),("onsite_strength",intstren),("if_periodic_x",false),("if_periodic_y",false)])

hilbert = ITensors.siteinds("S=1/2",lx,ly)
peps = PEPS(ComplexF64, hilbert; bond_dim=bond_dim)
QuantumNaturalfPEPS.multiply_algebraic_spectrum!(peps, 3.)

# [#Up = occupied, #Dn = empty]; sz_sector_counts(peps, Q) fixes total Sz = Q, i.e. n - lx*ly/2 here
max_counts = [n, lx*ly - n]
println("Sampling in the particle number N = $(params_dict["particles"]) sector, max_counts = $max_counts")

ham_opsum = long_range_HH_ham(params_dict)

# Generate Operators for QNG
Oks_and_Eks = QuantumNaturalfPEPS.generate_Oks_and_Eks(peps, ham_opsum; threaded=true, multiproc=true, max_counts)

# Setup the Integrator and Solver
learning_rate, sample_nr, maxiter = 0.01, 2000, 20
integrator = QuantumNaturalGradient.Euler(lr=learning_rate)
solver = QuantumNaturalGradient.EigenSolver()

# Define a Parameters object to be evolved
θ = QuantumNaturalGradient.Parameters(peps)

# Saving: like the TTNs, the PEPS goes in wavefuncfpeps-<params>.jld2 and the light data + metadata in fpeps-<params>.jld2
if_save_data = params_dict["if_save_data"]
dataloc = get(params_dict, "dataloc", get_folder_location("cluster-data/synth-dims/fpeps"))
filename_keys = ["Lx","Ly","particles","onsite_strength","lr","hopping_anisotropy","if_periodic_x","if_periodic_y"]
filename_dict = merge(Dict([(k,params_dict[k]) for k in filename_keys]), Dict([("bond_dim",bond_dim)]))
metadata = merge(params_dict, Dict([("bond_dim",bond_dim),("sample_nr",sample_nr),("maxiter",maxiter),("learning_rate",learning_rate),("solver",string(typeof(solver))),("integrator",string(typeof(integrator))),("location",dataloc),("niter",0)]))
filepath = if_save_data ? joinpath(dataloc, save_initial_fpeps(peps, metadata; location=dataloc, name="fpeps-" * make_parameters_filename(filename_dict) * ".jld2")) : nothing

# Callback: overwrite the saved PEPS and history at each iteration (returning false would stop evolve)
function callback(; misc, kwargs...)
    if_save_data && save_fpeps(peps, misc, filepath)
    return nothing
end

@time loss_value, trained_θ, misc = QuantumNaturalGradient.evolve(Oks_and_Eks, θ;
        integrator,
        verbosity=2,
        solver,
        sample_nr,
        maxiter,
        callback,)


println("Final energy: $(real(loss_value))")





















"fin"
####################################################
#=

This file runs a new fPEPS groundstate search (cluster script, started from fpeps/),
sequentially for values of Gamma using the previous final groundstate as the initial state
for the next Gamma. Saves every save_every iterations (default 10) unless if_save_data = false.

Depends on:
        other-funcs/include-other-files.jl (include_other_files)
        fpeps/load-fpeps.jl, of which
                build-ham.jl (get_fpeps_model_paras, long_range_HH_ham, particle_sector_counts, column_exclusion_peps)
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

nr_procs, nr_threads = 12, 6
addprocs(nr_procs, exeflags="-t $nr_threads")

@everywhere begin
    using QuantumNaturalfPEPS
    using LinearAlgebra
    BLAS.set_num_threads(1) # Each thread should only use one Blas thread
end
BLAS.set_num_threads(prod([nr_procs, nr_threads])) # For computing the double layers

args_dict = make_args_dict(ARGS)

# model
lx,ly,n = args_dict["Lx"],args_dict["Ly"],args_dict["particles"]
intstren = 50.0
exponent = args_dict["exponent"]
dataloc = choose_dataloc(args_dict)
if exponent == 0.0
    scaling_version = "flat"
elseif exponent == 3.0
    scaling_version = "dd"
else
    scaling_version = "powerlaw"
end

# optimisation (bond_dim also goes in params_dict, it's part of the filename)
bond_dim = 5
eigencut, contract_cutoff, sample_cutoff = 1e-6, 1e-6, 5e-4
learning_rate, sample_nr, maxiter = 0.05, 3000, 5000

# particle number conservation
max_counts = particle_sector_counts(lx, ly, n)
println("Sampling in the particle number N = $n sector, max_counts = $max_counts")

# range of Gamma values
gamma_vals = [35.0231,55.1244,86.8161,136.789,216.006,339.679,537.868,842.9,1329.95,2093.11]

previous_filepath = nothing

for (i,gamma) in enumerate(gamma_vals)

    magspac = exponent == 0.0 ? 1.0 : (intstren / gamma) .^ (1/exponent)


    params_dict = Dict([("outputlevel",1),("nrgtol",1e-2),("dataloc",dataloc),("scaling",scaling_version),("magnetic_spacing",magspac),("bond_dim",bond_dim),("lr","all"),("hopping_anisotropy",1.0),("Lx",lx),("Ly",ly),("particles",n),("if_save_data",true),("filling",0.5),("onsite_strength",intstren)])
    model_paras = get_fpeps_model_paras(params_dict)

    # initial state: use the previous final state, otherwise a random starting state
    if i == 1
        hilbert = ITensors.siteinds("S=1/2",lx,ly)
        peps = PEPS(ComplexF64, hilbert; bond_dim=bond_dim, contract_cutoff=contract_cutoff, sample_cutoff=sample_cutoff)
        QuantumNaturalfPEPS.multiply_algebraic_spectrum!(peps, 3.)
    else
        wavefunc_data, _ = read_data(add_wavefunc_to_filepath(previous_filepath); output_level=0)
        peps = wavefunc_data["peps"]
    end

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

    global previous_filepath = filepath

    # Callback: overwrite the saved PEPS and history every save_every iterations and at maxiter (filepath = nothing when not saving)
    callback = saving_fpeps_callback(filepath; save_every=model_paras["save_every"], nrgtol=model_paras["nrgtol"])

    loss_value, trained_θ, misc = QuantumNaturalGradient.evolve(Oks_and_Eks, θ;
            integrator,
            verbosity=2,
            solver,
            sample_nr,
            maxiter,
            callback,)

end


"fin"
















"fin"
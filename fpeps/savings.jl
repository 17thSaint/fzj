#####################################################
#=

This file contains saving, saving callbacks and restarts for fPEPS groundstate searches

Depends on (all loaded by fpeps/load-fpeps.jl):
	fpeps/build-ham.jl (get_normal_model_paras, long_range_HH_ham, for continue_fpeps)
	other-funcs/data-storage-funcs.jl (write_data, read_data, modify_data, get_folder_location, add_wavefunc_to_filepath, remove_wavefunc_from_filename)
	QuantumNaturalfPEPS (PEPS, generate_Oks_and_Eks), QuantumNaturalGradient (evolve, Parameters, Euler, EigenSolver)

Like the TTNs, the PEPS goes in wavefuncfpeps-<params>.jld2 and the light data (history, energy, rng)
+ metadata in fpeps-<params>.jld2, both overwritten by saving_fpeps_callback every save_every
(params_dict key, default 10) iterations and at maxiter

=#
######################################################


# cached double-layer envs are big and rebuilt on demand, so they're not saved
function strip_peps(peps)
    return PEPS(peps.tensors, peps.bond_dim; sample_dim=peps.sample_dim, sample_cutoff=peps.sample_cutoff,
        contract_dim=peps.contract_dim, contract_cutoff=peps.contract_cutoff,
        double_contract_dim=peps.double_contract_dim, double_contract_cutoff=peps.double_contract_cutoff,
        show_warning=peps.show_warning, mask=peps.mask)
end

# make filename
function make_fpeps_filename(model_parameters::Dict)
	# Start with the usual stuff for every filename
	num_particles = model_parameters["particles"]
	alpha = model_parameters["alpha"]
	onsite_strength = model_parameters["onsite_strength"]
	anis = model_parameters["hopping_anisotropy"]
	bond_dim = model_parameters["bond_dim"]
	
	filename_dict::Dict{String,Any} = Dict([("lr",model_parameters["lr"]),("particles",num_particles),("alpha",round(alpha,digits=4)),("onsite_strength",onsite_strength),("hopping_anisotropy",anis),("bond_dim",bond_dim)])

	filename_dict["scaling"] = model_parameters["scaling"]
	if model_parameters["scaling"] == "gaussian"
		filename_dict["sigma"] = model_parameters["sigma"]
	elseif model_parameters["scaling"] == "exp"
		filename_dict["corr_length"] = model_parameters["corr_length"]
	elseif model_parameters["scaling"] == "rydberg"
		filename_dict["blockade_radius"] = model_parameters["blockade_radius"]
	elseif model_parameters["scaling"] == "dd"
		filename_dict["magnetic_spacing"] = model_parameters["magnetic_spacing"]
	elseif model_parameters["scaling"] == "powerlaw"
		filename_dict["magnetic_spacing"] = model_parameters["magnetic_spacing"]
		filename_dict["exponent"] = model_parameters["exponent"]
	elseif model_parameters["scaling"] != "flat"
		error("ULR Scaling Type Not Recognized: $(model_parameters["scaling"])")
	end

	if model_parameters["ts"] != 1.0
		filename_dict["ts"] = model_parameters["ts"]
	end

	# the interaction axis, default y (rydberg is always both)
	if model_parameters["which_dir"] != "y" && model_parameters["scaling"] != "rydberg"
		filename_dict["which_dir"] = model_parameters["which_dir"]
	end

	filename_dict["Lx"] = model_parameters["Lx"]
	filename_dict["Ly"] = model_parameters["Ly"]

	if model_parameters["if_pinning"]
		filename_dict["if_pinning"] = true
	end

	return filename_dict
end


# main file first, so the wavefunc file takes its (possibly mk-renamed) name
function save_initial_fpeps(peps, metadata::Dict; kwargs...)
    location, filename = kwargs[:location], kwargs[:name]
    actual_filename = write_data(filename, Dict([("history",nothing),("energy",nothing),("rng",nothing)]), location, metadata)
    write_data("wavefunc"*actual_filename, Dict([("peps",strip_peps(peps))]), location, metadata)
    return actual_filename
end

# misc from QuantumNaturalGradient.evolve; history/energy/rng/niter are what misc_restart needs
function save_fpeps(peps, misc::Dict, filepath::String)
    modify_data(Dict([("peps",strip_peps(peps))]), add_wavefunc_to_filepath(filepath), "all_data")
    modify_data(Dict([("history",misc["history"]),("energy",misc["energy"]),("rng",misc["rng"])]), filepath, "all_data")
    modify_data(Dict([("niter",misc["niter"])]), filepath, "metadata")
end

# Callback for QuantumNaturalGradient.evolve that overwrites the saved PEPS and history every save_every iterations
# and at maxiter (returning false would stop evolve). Like the TTN saving observers everything it needs is carried
# with it and the PEPS comes from the optimisation state, no globals. filepath = nothing gives a callback that does nothing
function saving_fpeps_callback(filepath::Union{String,Nothing}; save_every::Int=10)
    return function (; state, misc, kwargs...)
        if !isnothing(filepath) && (misc["niter"] % save_every == 0 || misc["niter"] >= state.maxiter)
            save_fpeps(state.θ.obj, misc, filepath)
        end
        return nothing
    end
end

# Continue a saved groundstate search for extra_iters more iterations. filepath is the fpeps-*.jld2 or its
# wavefuncfpeps-*.jld2 (a bare filename is looked for in cluster-data/synth-dims/fpeps). The PEPS, history,
# energy, rng and niter are read back and the Hamiltonian rebuilt from the metadata; sample_nr, the Euler
# learning rate and save_every default to the saved ones. With if_save_data the same two files keep being overwritten
function continue_fpeps(filepath::String, extra_iters::Int; if_save_data::Bool=true, save_every=nothing, sample_nr=nothing, learning_rate=nothing,
        integrator=nothing, solver=QuantumNaturalGradient.EigenSolver(), threaded=true, multiproc=true, verbosity=2)

    occursin("/",filepath) ? nothing : filepath = joinpath(get_folder_location("cluster-data/synth-dims/fpeps"), filepath)
    filepath = joinpath(dirname(filepath), remove_wavefunc_from_filename(basename(filepath)))
    data, metadata = read_data(filepath)
    wavefunc_data, _ = read_data(add_wavefunc_to_filepath(filepath))
    peps = wavefunc_data["peps"]
    params_dict = Dict{String,Any}(metadata)
    # runs saved before get_normal_model_paras stored only the raw params_dict, so derive t_x, long_range_strengths, ...
    # from it. Newer runs already have them (runs from 2026-10-29 also keep lr, ts, filling and the scaling inputs)
    haskey(params_dict, "long_range_strengths") || (params_dict = merge(params_dict, get_normal_model_paras(params_dict)))

    # niter = 0 means no iteration finished before the run stopped, so start the optimisation fresh
    niter = params_dict["niter"]
    misc_restart = niter == 0 ? nothing : Dict([("history",data["history"]),("energy",data["energy"]),("rng",data["rng"]),("niter",niter)])
    maxiter = niter + extra_iters

    lx,ly,n = params_dict["Lx"],params_dict["Ly"],params_dict["particles"]
    max_counts = [n, lx*ly - n]
    ham_opsum = long_range_HH_ham(params_dict)
    Oks_and_Eks = QuantumNaturalfPEPS.generate_Oks_and_Eks(peps, ham_opsum; threaded, multiproc, max_counts)

    sample_nr = isnothing(sample_nr) ? params_dict["sample_nr"] : sample_nr
    learning_rate = isnothing(learning_rate) ? params_dict["learning_rate"] : learning_rate
    save_every = isnothing(save_every) ? get(params_dict, "save_every", 10) : save_every
    if isnothing(integrator)
        occursin("Euler", params_dict["integrator"]) ? nothing : @warn "saved integrator was $(params_dict["integrator"]), continuing with Euler"
        integrator = QuantumNaturalGradient.Euler(lr=learning_rate)
    end
    println("Continuing $(basename(filepath)) from niter = $niter to maxiter = $maxiter")

    if if_save_data
        modify_data(Dict([("maxiter",maxiter),("sample_nr",sample_nr),("learning_rate",learning_rate),("save_every",save_every),
            ("solver",string(typeof(solver))),("integrator",string(typeof(integrator)))]), filepath, "metadata")
    end

    callback = saving_fpeps_callback(if_save_data ? filepath : nothing; save_every)

    θ = QuantumNaturalGradient.Parameters(peps)
    @time loss_value, trained_θ, misc = QuantumNaturalGradient.evolve(Oks_and_Eks, θ;
            integrator, verbosity, solver, sample_nr, maxiter, callback, misc_restart)

    println("Final energy: $(real(loss_value))")
    return loss_value, trained_θ.obj, misc
end

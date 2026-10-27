#####################################################
#=

This file contains set-up for fPEPS simulations

Depends on:
	other-funcs/basic-2d-stuff.jl (coordinate, find_dist, check_fluxes, long_range_scaling, get_hopping_strengths)
	ITensors (OpSum)

All parameters are read from a String-keyed dict, axes are x (1st index, 1:Lx) and
y (2nd index, 1:Ly), matching hilbert = siteinds("S=1/2", Lx, Ly)

Hardcore bosons on S=1/2 sites: Up = occupied, Dn = empty, so b† = S+, b = S-, n = ProjUp

=#
######################################################


include("../other-funcs/basic-2d-stuff.jl")


# Peierls hopping from site along axis (1 = x, 2 = y), Landau gauge with the phase on the
# flux_direction hoppings: x hops pick up exp(2πiα(y-1)), y hops exp(2πiα(x-1))
function get_hopping_coeff(site,axis::Int,t_x::Float64,t_y::Float64,alpha::Float64,flux_direction::String)
	if axis == 1
		coeff = -t_x
		flux_direction == "x" ? coeff *= exp(im*2*pi*alpha*(site[2]-1)) : nothing
	else
		coeff = -t_y
		flux_direction == "y" ? coeff *= exp(im*2*pi*alpha*(site[1]-1)) : nothing
	end
	return coeff
end

function long_range_HH_ham(model_paras::Dict{String,Any})

	Lx,Ly = model_paras["Lx"],model_paras["Ly"]
	L = (Lx,Ly)
	num_sites = Lx*Ly
	outputlevel = get(model_paras, "outputlevel", 1)
	outputlevel > 0 ? println("Lx = ",Lx,", Ly = ",Ly) : nothing

	# lattice
	if_periodic_x = get(model_paras, "if_periodic_x", false)
	if_periodic_y = get(model_paras, "if_periodic_y", false)
	if_per = (if_periodic_x,if_periodic_y)

	# hopping and flux
	t_strength::Float64 = get(model_paras, "ts", 1.0)
	t_x,t_y = get_hopping_strengths(t_strength, Float64(get(model_paras, "hopping_anisotropy", 1.0)))
	twist_angle = get(model_paras, "twist_angle", [0.0,0.0])
	if haskey(model_paras, "alpha")
		alpha = Float64(model_paras["alpha"])
	else
		particles = get(model_paras, "particles", Int(ceil(Lx/2)))
		filling = get(model_paras, "filling", 0.5)
		alpha = filling == 0.0 ? 0.0 : particles/(filling*(Lx - !if_periodic_x)*(Ly - !if_periodic_y))
	end
	get(model_paras, "no_magF", false) ? alpha = 0.0 : nothing
	flux_direction = get(model_paras, "flux_direction", "y")
	if if_periodic_y && !if_periodic_x
		flux_direction = "y"
	elseif if_periodic_x && !if_periodic_y
		flux_direction = "x"
	end
	if get(model_paras, "if_check_fluxes", true) && alpha != 0.0
		flux_direction = check_fluxes(alpha,Lx,Ly,if_periodic_x,if_periodic_y,flux_direction; output_level=outputlevel)
	end

	# interactions, which_dir is the axis the long-range interaction acts along ("x", "y" or "both")
	onsite_strength::Float64 = get(model_paras, "onsite_strength", 0.0)
	scaling = get(model_paras, "scaling", "flat")
	which_dir = scaling == "rydberg" ? "both" : get(model_paras, "which_dir", "y")
	interaction_axis_length = which_dir == "x" ? Lx : which_dir == "y" ? Ly : max(Lx,Ly)
	scaling_distance = get(model_paras, "lr", 0)
	scaling_distance == "all" ? scaling_distance = interaction_axis_length - 1 : nothing
	scaling_kwargs = Dict(Symbol(k) => v for (k,v) in model_paras)
	long_range_strengths = long_range_scaling(scaling_distance,interaction_axis_length,onsite_strength; scaling_kwargs...)
	outputlevel > 0 ? display(long_range_strengths) : nothing

	# on-site potentials
	chem_strength = get(model_paras, "chem_strength", 0.0)
	if_pinning = get(model_paras, "if_pinning", false)
	pinning_strength::Float64 = get(model_paras, "pinning_strength", 1e-4)

	ham = ITensors.OpSum()

	for i in 1:num_sites
		starting_site = coordinate(i,Lx,Ly)
		for axis in [1,2]
			ending_site = collect(starting_site)
			ending_site[axis] += 1
			twist = 1.0

			# enforce boundary conditions
			if ending_site[axis] > L[axis]
				if_per[axis] ? nothing : continue
				ending_site[axis] = 1
				twist = exp(im*twist_angle[axis]*2*pi)
			end

			coeff = get_hopping_coeff(starting_site,axis,t_x,t_y,alpha,flux_direction) * twist
			ham += (coeff,"S+",starting_site,"S-",Tuple(ending_site))
			ham += (conj(coeff),"S+",Tuple(ending_site),"S-",starting_site)
		end
	end

	# long_range_strengths[d+1] couples each pair a distance d apart along the interaction axis, open boundaries,
	# the on-site d = 0 term drops out for hardcore bosons
	for i in 1:num_sites, j in i+1:num_sites
		s1,s2 = coordinate(i,Lx,Ly),coordinate(j,Lx,Ly)
		_,(dx,dy),_ = find_dist(s1,s2,(Lx,Ly))
		on_axis = which_dir == "x" ? dy == 0 : which_dir == "y" ? dx == 0 : dx == 0 || dy == 0
		on_axis || continue
		stren = long_range_strengths[Int(dx + dy) + 1]
		stren == 0.0 && continue
		ham += (stren,"ProjUp",s1,"ProjUp",s2)
	end

	if chem_strength != 0.0
		for i in 1:num_sites
			ham -= (chem_strength,"ProjUp",coordinate(i,Lx,Ly))
		end
	end

	if if_pinning
		ham += (pinning_strength,"ProjUp",(1,1))
	end

	return ham
end

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
	longrange_dist = model_parameters["lr"]
	num_particles = model_parameters["particles"]
	alpha = model_parameters["alpha"]
	onsite_strength = model_parameters["onsite_strength"]
	anis = model_parameters["hopping_anisotropy"]
	bond_dim = model_parameters["bond_dim"]
	
	filename_dict::Dict{String,Any} = Dict([("lr",longrange_dist),("particles",num_particles),("alpha",round(alpha,digits=4)),("onsite_strength",onsite_strength),("hopping_anisotropy",anis),("bond_dim",bond_dim)])

	if model_parameters["scaling"] != "flat"
		filename_dict["scaling"] = model_parameters["scaling"]
		if model_parameters["scaling"] == "gaussian"
			filename_dict["sigma"] = model_parameters["sigma"]
		elseif model_parameters["scaling"] == "exp"
			filename_dict["corr_length"] = model_parameters["corr_length"]
		elseif model_parameters["scaling"] == "rydberg"
			filename_dict["blockade_radius"] = model_parameters["blockade_radius"]
		else
			error("ULR Scaling Type Not Recognized: $(model_parameters["scaling"])")
		end
	end

	if model_parameters["ts"] != 1.0
		filename_dict["ts"] = model_parameters["ts"]
	end

	filename_dict["Lx"] = model_parameters["lattice_size"][1]
	filename_dict["Ly"] = model_parameters["lattice_size"][2]

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

# Continue a saved groundstate search for extra_iters more iterations. filepath is the fpeps-*.jld2 or its
# wavefuncfpeps-*.jld2 (a bare filename is looked for in cluster-data/synth-dims/fpeps). The PEPS, history,
# energy, rng and niter are read back and the Hamiltonian rebuilt from the metadata; sample_nr and the Euler
# learning rate default to the saved ones. With if_save_data the same two files keep being overwritten
function continue_fpeps(filepath::String, extra_iters::Int; if_save_data::Bool=true, sample_nr=nothing, learning_rate=nothing,
        integrator=nothing, solver=QuantumNaturalGradient.EigenSolver(), threaded=true, multiproc=true, verbosity=2)

    occursin("/",filepath) ? nothing : filepath = joinpath(get_folder_location("cluster-data/synth-dims/fpeps"), filepath)
    filepath = joinpath(dirname(filepath), remove_wavefunc_from_filename(basename(filepath)))
    data, metadata = read_data(filepath)
    wavefunc_data, _ = read_data(add_wavefunc_to_filepath(filepath))
    peps = wavefunc_data["peps"]
    params_dict = Dict{String,Any}(metadata)

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
    if isnothing(integrator)
        occursin("Euler", params_dict["integrator"]) ? nothing : @warn "saved integrator was $(params_dict["integrator"]), continuing with Euler"
        integrator = QuantumNaturalGradient.Euler(lr=learning_rate)
    end
    println("Continuing $(basename(filepath)) from niter = $niter to maxiter = $maxiter")

    if if_save_data
        modify_data(Dict([("maxiter",maxiter),("sample_nr",sample_nr),("learning_rate",learning_rate),
            ("solver",string(typeof(solver))),("integrator",string(typeof(integrator)))]), filepath, "metadata")
    end

    function callback(; state, misc, kwargs...)
        if_save_data && save_fpeps(state.θ.obj, misc, filepath)
        return nothing
    end

    θ = QuantumNaturalGradient.Parameters(peps)
    @time loss_value, trained_θ, misc = QuantumNaturalGradient.evolve(Oks_and_Eks, θ;
            integrator, verbosity, solver, sample_nr, maxiter, callback, misc_restart)

    println("Final energy: $(real(loss_value))")
    return loss_value, trained_θ.obj, misc
end

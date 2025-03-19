#####################################################
#=

This file contains functions for pfaffian searches

=#
######################################################


include("../other-funcs/include-other-files.jl")
include_other_files(["synth-dims/long-range-ttn.jl"])
include_other_files(["other-funcs/basic-2d-plottings.jl","review-practice-codes/plottings.jl"])

function c3(wavefunc::TTN.TreeTensorNetwork; kwargs...)

    lx,ly = get_lattice_dims(wavefunc)

    c3_mat::Matrix{ComplexF64} = zeros(ComplexF64,lx,ly)
    for s_phys in 1:lx
        for s_synth in 1:ly
            c3_mat[s_phys,s_synth] += expect(wavefunc,"N * N * N",(s_phys,s_synth))
            c3_mat[s_phys,s_synth] += -3*expect(wavefunc,"N * N",(s_phys,s_synth))
            c3_mat[s_phys,s_synth] += 2*expect(wavefunc,"N",(s_phys,s_synth))
        end
    end

    return c3_mat,sum(c3_mat)

end

function pfaffian_ham(latsize::Tuple; kwargs...)

    if_periodic::Vector{Bool} = [kwargs[:if_periodic_phys],kwargs[:if_periodic_synth]]
    t_strength::Float64 = kwargs[:ts]
    phi::Float64 = kwargs[:phi]
    nbody_inter::Int = kwargs[:nbody_inter]
    interaction_strength::Float64 = kwargs[:interaction_strength]
    twist_angle::Vector{Float64} = kwargs[:twist_angle]

    hopping = TTN.OpSum()
    for s_phys in 1:latsize[1]
        for s_synth in 1:latsize[2]
            starting_site = [s_phys,s_synth]
            twist = 0
            for which_axis in [1,2]
                ending_site = starting_site .+ ((which_axis == 1,which_axis == 2))

                if ending_site[which_axis] > latsize[which_axis] 
                    if if_periodic[which_axis]
                        ending_site[which_axis] = 1
                        twist = 1
                    else
                        continue
                    end
                end

                if ending_site[which_axis] < 1
                    if if_periodic[which_axis]
                        ending_site[which_axis] = latsize[which_axis]
                        twist = 2
                    else
                        continue
                    end
                end

                if starting_site[1] == ending_site[1] # Synthetic Dimension Hopping
                    coeff = t_strength
                elseif starting_site[2] == ending_site[2] # Physical Dimension Hopping
                    coeff = t_strength * exp(im*2*pi*(phi*starting_site[2]))
                end

                twist == 1 && (coeff *= exp(im*twist_angle[which_axis]*2*pi))
				twist == 2 && (coeff *= exp(-im*twist_angle[which_axis]*2*pi))

                hopping += (coeff,"Adag",Tuple(starting_site),"A",Tuple(ending_site))
                hopping += (conj(coeff),"Adag",Tuple(ending_site),"A",Tuple(starting_site))
            end
        end
    end

    interaction = TTN.OpSum()
    if nbody_inter == 2
        for s_phys in 1:latsize[1]
            for s_synth in 1:latsize[2]
                interaction += (interaction_strength/2,"N * N",(s_phys,s_synth))
                interaction += (-interaction_strength/2,"N",(s_phys,s_synth))
            end
        end
    elseif nbody_inter == 3
        for s_phys in 1:latsize[1]
            for s_synth in 1:latsize[2]
                interaction += (interaction_strength/6,"N * N * N",(s_phys,s_synth))
                interaction += (-3*interaction_strength/6,"N * N",(s_phys,s_synth))
                interaction += (2*interaction_strength/6,"N",(s_phys,s_synth))
            end
        end
    end

    return hopping + interaction

end

function make_pfaffian_filename(model_parameters::Dict)
	layer_count = model_parameters["layers"]
	num_particles = model_parameters["particles"]
	alpha = model_parameters["alpha"]
	if_periodic_phys = model_parameters["if_periodic_phys"]
	if_periodic_synth = model_parameters["if_periodic_synth"]
	interaction_strength = model_parameters["interaction_strength"]
    nbody_inter = model_parameters["nbody_inter"]
	
	filename_dict = Dict([("layers",layer_count),("particles",num_particles),("alpha",round(alpha,digits=4)),("if_periodic_phys",if_periodic_phys),("if_periodic_synth",if_periodic_synth),("nbody_inter",nbody_inter),("interaction_strength",interaction_strength)])

	if model_parameters["ts"] != 1.0
		filename_dict["ts"] = model_parameters["ts"]
	end

	if model_parameters["twist_angle"] != [0.0,0.0]
		filename_dict["twist_angle1"] = model_parameters["twist_angle"][1]
		filename_dict["twist_angle2"] = model_parameters["twist_angle"][2]
	end

	return make_parameters_filename(filename_dict)
end

function get_pfaffian_model_params(params_dict::Dict)

	# DMRG parameters
	sweep_type::String = get(params_dict, "sweep_type", "dmrg")
	nrgtol = get(params_dict, "nrgtol", 5E-5)
	cutoff = get(params_dict, "cutoff", 1E-8)
	evolve = get(params_dict, "evolve", true)
	expander_fraction = get(params_dict, "expander_fraction", 0.0)
	expan = TTN.DefaultExpander(expander_fraction)
	noise = get(params_dict, "noise", [0.0])
	syms = get(params_dict, "syms", true)
	nswps = get(params_dict, "num_sweeps", 100)
	output_level = get(params_dict, "output_level", 1)
	seed_ttn = get(params_dict, "seed_ttn", nothing)


	# Lattice/TTN Parameters
	layer_count = Int(get(params_dict, "layers", 4))
	mdim = get(params_dict, "mdim", 300)
	if_periodic_phys = get(params_dict, "if_periodic_phys", false)
	if_periodic_synth = get(params_dict, "if_periodic_synth", false)
	max_occ = get(params_dict, "max_occ", 4)

	# Get Lattice parameters whose values depend on other parameters
	if layer_count % 2 == 0
		phys_edge_length,synth_edge_length = Int(sqrt(2^layer_count)),Int(sqrt(2^layer_count))
		num_particles = get(params_dict, "particles", Int(phys_edge_length/2))
	else
		phys_edge_length,synth_edge_length = Int(sqrt(2^(layer_count+1))),Int(sqrt(2^(layer_count+1))/2)
		num_particles = get(params_dict, "particles", Int(sqrt(2^(layer_count+1))/2))
	end

	make_smaller_lattice = get(params_dict, "make_smaller_lattice", [phys_edge_length,synth_edge_length])
	if make_smaller_lattice != [phys_edge_length,synth_edge_length]
		phys_edge_length,synth_edge_length = make_smaller_lattice
	end


	# Hamiltonian parameters
	alpha = get(params_dict, "alpha", nothing)

	hopping_amplitude = get(params_dict, "ts", 1.0)
	anis = get(params_dict, "hopping_anisotropy", 1.0)

	interaction_strength::Float64 = get(params_dict, "interaction_strength", 0.0)
    nbody_inter = get(params_dict, "nbody_inter", 2)

	twist_angle::Vector{Float64} = [get(params_dict, "tw1", 0.0),get(params_dict, "tw2", 0.0)]

	if isnothing(alpha)
		filling = get(params_dict, "filling", 1.0)
		phys_shift,synth_shift = !if_periodic_phys,!if_periodic_synth
		alpha = num_particles/(filling*(phys_edge_length - phys_shift)*(synth_edge_length - synth_shift))
		filling == 0.0 ? alpha = 0.0 : nothing
	else
		mag_off = alpha == 0.0
	end
	if_check_fluxes = get(params_dict, "if_check_fluxes", true)
	if_check_fluxes ? check_fluxes(alpha,phys_edge_length,synth_edge_length,if_periodic_phys,if_periodic_synth,"phys"; output_level=output_level) : nothing


	# What to calculate
	if_redo = get(params_dict, "if_redo", false)
	if_densmat = get(params_dict, :if_densmat, true)
	save_data = get(params_dict, "if_save_data", true)
	if_cluster = any([occursin("local",pwd()),occursin("Local",pwd()),occursin("geraghty",pwd())])
	if_continuous_saving = get(params_dict,"if_continuous_saving",if_cluster || layer_count >= 7)
	save_data ? nothing : if_continuous_saving = false
	es_count = get(params_dict, "es_count", 0)
	
	all_measurements = get(params_dict, "all_measurements", String[])
	measurement_functions::Vector{NamedTuple} = construct_measurement_info(all_measurements)
	measurements::Dict{String,Any} = Dict()
	for info_tuple in measurement_functions
		measurements[info_tuple[:name]] = nothing
	end

	dataloc = get_folder_location("cluster-data/pfaffian")
	loc = get(params_dict, "dataloc", dataloc)
	


	# hardware parameters
	if_gpu = get(params_dict, "if_gpu", false)

	
	model_paras_dict = Dict("hopping_anisotropy"=>anis,
						"layers"=>layer_count,
						"particles"=>num_particles,
						"ts"=>hopping_amplitude,
						"syms"=>syms,
						"cutoff"=>cutoff,
						"seed_ttn"=>seed_ttn,
						"twist_angle"=>twist_angle,
						"if_continuous_saving"=>if_continuous_saving,
						"output_level"=>output_level,
						"nrgtol"=>nrgtol,
						"if_densmat"=>if_densmat,
						"if_redo"=>if_redo,
						"restricted_size"=>make_smaller_lattice,
						"if_periodic_phys"=>if_periodic_phys,
						"if_periodic_synth"=>if_periodic_synth,
						"alpha"=>alpha,
						"interaction_strength"=>interaction_strength,
                        "nbody_inter"=>nbody_inter,
						"if_gpu"=>if_gpu,
						"noise"=>noise,
						"if_save_data"=>save_data,
						"if_sweep"=>evolve,
						"sweep_type"=>sweep_type,
						"expander"=>expan,
						"max_occ"=>max_occ,
						"mdim"=>mdim,
						"num_sweeps"=>nswps,
						"phi"=>alpha,
						"output_level"=>0,
						"location"=>loc)
		
	if length(all_measurements) > 0
		model_paras_dict["measurements"] = measurements
		model_paras_dict["measurement_functions"] = measurement_functions
	end
	filename = make_pfaffian_filename(model_paras_dict)
	model_paras_dict["name"] = "ttn-"*filename
	
	return dict_to_symbols(model_paras_dict)
end

function run_pfaffian_generic(params_dict::Dict)

	if_find_data = get(params_dict, "if_find_data", true)

	model_paras = get_pfaffian_model_params(params_dict)
	metadata_dict = named_tuple_to_dict(model_paras)
    lx,ly = get_lattice_dims_from_layers(model_paras[:layers])

	es_count = get(params_dict, "es_count", 0)
	if_redo = get(params_dict, "if_redo", false)

		#
	println(model_paras[:name])
	filename_dict = get_params_dict_from_filename(model_paras[:name])
	println("Location for data is $(model_paras[:location])")
	if_exists,found_data = if_find_data ? check_data_exists(filename_dict,"ttn"; location=model_paras[:location],output_level=false) : (false,nothing)

	if if_exists
		if_wavefunc = !isnothing(found_data[1]["ttn"])
		if es_count > 0 # when need excited states start by counting how many inside the data file
			count_found_states = length(findall(x -> occursin("ttn",x),collect(keys(found_data[1]))))

			if count_found_states < es_count + 1 # found states less than asked for count means run for higher states
				println("Not Enough States in Data File, Running for $(es_count - count_found_states + 1) more States")
				ortho_states = Vector{TTN.TreeTensorNetwork}(undef,count_found_states)
				ortho_states[1] = found_data[1]["ttn"]
				for i in 2:count_found_states
					ortho_states[i] = found_data[1]["ttn_$(i-1)"]
				end
				model_paras[:ham] = found_data[2]["ham"]
				metadata_dict["ham"] = model_paras[:ham]
				all_states, hamilt, all_obs, all_densmats, all_runtimes = find_excited_states(params_dict["layers"],es_count,model_paras[:particles],ortho_states; model_paras...,metadata=metadata_dict)

			else # found states is less than or equal to asked for count means use found states
				println("Found Data")
				ortho_states = Vector{TTN.TreeTensorNetwork}(undef,es_count+1)
				densmats = Vector{Matrix{ComplexF64}}(undef,es_count+1)
				obss = Vector{TTN.AbstractObserver}(undef,es_count+1)
				runtimes = zeros(es_count+1)
				if_wavefunc ? ortho_states[1] = found_data[1]["ttn"] : nothing
				densmats[1] = found_data[1]["densmat"]
				obss[1] = found_data[2]["observer"]
				runtimes[1] = found_data[2]["runtime"]
				for i in 2:es_count+1
					if_wavefunc ? ortho_states[i] = found_data[1]["ttn_$(i-1)"] : nothing
					densmats[i] = found_data[1]["densmat_$(i-1)"]
					obss[i] = found_data[2]["observer_$(i-1)"]
					runtimes[i] = found_data[2]["runtime_$(i-1)"]
				end
				
				all_states, hamilt, all_obs, all_densmats, all_runtimes = ortho_states, found_data[2]["ham"], obss, densmats, runtimes
			end

		else # if only ask for the ground state then if data if found then have all needed results
			println("Found Data")
			og_ttn = if_wavefunc ? found_data[1]["ttn"] : nothing
			gs_dens = found_data[1]["densmat"]
			gs_obs = found_data[2]["observer"]
			hamilt = found_data[2]["ham"]
			gs_runtime = found_data[2]["runtime"]
			gs_sp = nothing
		end
	else # if no data found then run from scratch starting from the ground state
		println("Starting Script using $(model_paras[:particles]) particles on $(2^model_paras[:layers]) sites $(model_paras[:nbody_inter])-body interactions at strength $(model_paras[:interaction_strength]) with Flux = $(round(model_paras[:alpha],digits=4)), Bond Dim = $(model_paras[:mdim])")

		starting = time()
		net = build_HH_net(model_paras)
		ham = pfaffian_ham((lx,ly); model_paras...)
		metadata_dict["ham"] = ham
		metadata_dict["net"] = net
		if es_count > 0
			all_states, hamilt, all_obs, all_densmats, all_runtimes = find_spectrum(model_paras,es_count,metadata_dict)
		else
			og_ttn, hamilt, gs_sp, gs_obs, gs_runtime, gs_dens = find_spectrum(model_paras,es_count,metadata_dict)
		end
		total_time = time() - starting
		println("Running time = $total_time")
		
	end

	if es_count < 1
		return gs_sp.ttn, hamilt, gs_obs, gs_dens, gs_runtime
	else
		return all_states, hamilt, all_obs, all_densmats, all_runtimes
	end
end

lx,ly,n = 4,4,4
layers = Int(log(2,lx*ly))
if_periodic = false

strens = range(0.0,1.0,length=11)

for (idx,stren) in enumerate(strens)
    pdict = Dict([("if_gpu",false),("particles",n),("nbody_inter",2),("es_count",0),("if_check_fluxes",false),("filling",1.0),("layers",layers),("expander_fraction",10),("mdim",100),("if_save_data",false),("interaction_strength",stren),("if_periodic_phys",if_periodic),("if_periodic_synth",if_periodic)])

    psis, hamilt, obs, rhos, rts = run_pfaffian_generic(pdict)

    c3mat,c3val = c3(psis)
    scatter(stren,real(c3val),c="b")
    xlabel("Interaction Strength")
    ylabel("C3 Value")
end








































"fin"
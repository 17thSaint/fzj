#####################################################
#=

This file contains the simple observable functions for ED

Depends on:
    execute-ed.jl

=#
######################################################


function find_center()
	all_folders = split(pwd(),"/")
	if "fzj" in all_folders
		return "fzj"
	elseif "local" in all_folders
		return all_folders[findfirst(x -> all_folders[x] == "local",1:length(all_folders))+1]
	elseif "Local" in all_folders
		return all_folders[findfirst(x -> all_folders[x] == "Local",1:length(all_folders))+1]
	else
		println("Not sure where the center is: $(pwd())")
	end
end

function include_other_files(all_files,output_level=0)
	center = find_center()
	get_to_fzj = split(pwd(),center)[1]
	if typeof(all_files) == String
		all_files = [all_files]
	end
	for file in all_files
		occursin("main-git",pwd()) ? include(get_to_fzj * center * "/main-git/" * file) : include(get_to_fzj * center * "/" * file)
		output_level > 0 ? println("Included $file") : nothing
	end
end

include_other_files(["exact-diag/execute-ed.jl"])

function pulse_function(nsteps::Int,dt::Float64; kwargs...)

    starting_value::Float64 = 0.01
    ending_value::Float64 = 1.0

    starting_time::Float64 = 0.0
    ending_time::Float64 = kwargs[:ending_time]

    steps_until_end::Int = Int(ceil(ending_time / dt))

    pulse_ramp = kwargs[:pulse_ramp][:]

    return vcat(pulse_ramp, ending_value .* ones(nsteps - steps_until_end + 1))
end

function compute_fidelity(pulses,parameters_dictionary)

    parameters_dictionary["output_level"] = 0
    
    # first find starting and final groundstates
    parameters_dictionary["tx"] = 0.01
    #startingGS_pdict = get_normal_model_params_ed(parameters_dictionary)
    startingGS_states,_,_,startingGS_filepath,startingGS_if_found,startingGS_lattice_params,startingGS_hamilt_params = run_normal_ed(parameters_dictionary; output_level=0)

    parameters_dictionary["tx"] = 1.0
    #finalGS_pdict = get_normal_model_params_ed(parameters_dictionary)
    finalGS_states,_,_,finalGS_filepath,finalGS_if_found,finalGS_lattice_params,finalGS_hamilt_params = run_normal_ed(parameters_dictionary; output_level=0)

    # then find the time-evolved state
    time_running_args = (nev=1,output_level=0,if_instant_gs=false,)
    tmax_global = 25.0
    dt_global = 0.05
    ramptime = 2.0
    rampsteps = Int(2*ceil(ramptime / dt_global)-1)

    @assert length(pulses[1]) == rampsteps "The length of the pulse must match the number of ramp steps"

    tevo_params = Dict([ ("tx",(pulse_function,ramptime,pulses[1])) ])
    tevo_gs,tevo_dict,intspec = run_timeevo(startingGS_states[1],tevo_params,startingGS_lattice_params,startingGS_hamilt_params; time_running_args...)
    println("Finished Time Evolution")
    
    return abs2(dot(tevo_gs[:,end-1],finalGS_states[1]))
end































"fin"
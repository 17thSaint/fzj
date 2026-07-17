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

# bins_number in config_txRamp.py must equal ceil(2*ramptime/dt) + 1 to match
# pulse_ramp's (time-evolution.jl) half-step sampling grid

# pulse-independent part of the tx-ramp figure of merit: the endpoint ground states.
# Run once per optimization (the QuOCS FoM object caches the returned tuple) instead
# of re-diagonalizing at every function evaluation.
function setup_tx_ramp(parameters_dictionary)

    parameters_dictionary["output_level"] = 0

    parameters_dictionary["tx"] = 0.001
    startingGS_states,_,_,_,_,startingGS_lattice_params,startingGS_hamilt_params = run_normal_ed(parameters_dictionary; output_level=0)

    parameters_dictionary["tx"] = 1.0
    finalGS_states,_,_,_,_,_,_ = run_normal_ed(parameters_dictionary; output_level=0)

    starting_states = [Vector{ComplexF64}(startingGS_states[1])]
    target_states = [Vector{ComplexF64}(finalGS_states[1])]

    return (starting_states,target_states,startingGS_lattice_params,startingGS_hamilt_params)
end

function compute_fidelity(pulses,parameters_dictionary,setup)

    starting_states,target_states,lattice_params,hamilt_params = setup

    # ramp settings, overridable from the config (defaults preserve the historical
    # hardcoded values this FoM was introduced with)
    dt = get(parameters_dictionary,"dt",0.05)
    ramptime = get(parameters_dictionary,"ramptime",2.0)
    time_running_args = (nev=1,output_level=0,if_instant_gs=false,if_save_data=false)

    final_states = run_ramp_stages(starting_states,[("tx",collect(pulses[1]),ramptime)],lattice_params,hamilt_params,dt; time_running_args...)

    return abs2(dot(final_states[1],target_states[1]))
end

compute_fidelity(pulses,parameters_dictionary) = compute_fidelity(pulses,parameters_dictionary,setup_tx_ramp(parameters_dictionary))

function groundstate_manifold_fidelity(comparison_states::Vector{Vector{ComplexF64}},target_states::Vector{Vector{ComplexF64}})
    fidelity_matrix = zeros(ComplexF64,length(comparison_states),length(target_states))
    for i in 1:length(comparison_states)
        for j in 1:length(target_states)
            fidelity_matrix[i,j] = adjoint(comparison_states[i]) * target_states[j]
        end
    end
    return 0.5 * tr(adjoint(fidelity_matrix)*fidelity_matrix)
end

function groundstate_manifold_fidelity(comparison_states::Vector{SparseVector{ComplexF64,Int64}},target_states::Vector{Vector{ComplexF64}})
    fidelity_matrix = zeros(ComplexF64,length(comparison_states),length(target_states))
    for i in 1:length(comparison_states)
        for j in 1:length(target_states)
            fidelity_matrix[i,j] = adjoint(vec(comparison_states[i])) * target_states[j]
        end
    end
    return 0.5 * tr(adjoint(fidelity_matrix)*fidelity_matrix)
end





























"fin"
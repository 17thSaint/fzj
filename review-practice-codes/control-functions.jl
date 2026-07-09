#####################################################
#=

This file contains the simple observable functions for TTNs

Depends on:
    review-practice-codes/ttn.jl
    review-practice-codes/observables.jl

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

include_other_files(["review-practice-codes/ttn.jl","review-practice-codes/observables.jl"])

function groundstate_manifold_fidelity(comparison_states::Vector{TTN.TreeTensorNetwork}, target_states::Vector{TTN.TreeTensorNetwork})
    fidelity_matrix = zeros(ComplexF64,length(comparison_states),length(target_states))
    for i in 1:length(comparison_states)
        for j in 1:length(target_states)
            fidelity_matrix[i,j] = TTN.inner(comparison_states[i], target_states[j])
        end
    end
    return 0.5 * tr(adjoint(fidelity_matrix)*fidelity_matrix)
end


















"fin"
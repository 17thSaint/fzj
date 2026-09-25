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

# the interaction-profile settings a figure of merit may carry, forwarded to the ED so a non-flat
# scaling (exp / gaussian / rydberg / dd) reaches get_normal_model_params_ed instead of silently
# falling back to flat
const OPTIMIZATION_PROFILE_KEYS = ("scaling_type","corr_length","sigma","blockade_radius","magnetic_spacing","trap_frequency","interaction_strength")

# The ED parameter dict an optimization setup diagonalizes with: the lattice from the figure of
# merit's parameters_dictionary (as quocs_common.JuliaFoM serializes it), the defaults every
# optimization here shares, the profile keys above when present, then the caller's overrides
# as "key"=>value pairs.
function optimization_ed_params(parameters_dictionary,overrides::Pair...)
    pd = parameters_dictionary
    ed = Dict{String,Any}(
        "output_level"=>0,"Lx"=>Int(pd["Lx"]),"Ly"=>Int(pd["Ly"]),"N"=>Int(pd["N"]),
        "lr"=>get(pd,"lr","all"),"if_periodic_x"=>get(pd,"if_periodic_x",true),"if_periodic_y"=>get(pd,"if_periodic_y",true),
        "hopping_anisotropy"=>1.0,"filling"=>0.5,"nev"=>Int(pd["speccount"]),"if_find_data"=>false,"if_save_data"=>false,
    )
    for k in OPTIMIZATION_PROFILE_KEYS
        haskey(pd,k) && (ed[k] = pd[k])
    end
    for (k,v) in overrides
        ed[k] = v
    end
    return ed
end

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

# Population of the target manifold averaged over the comparison states,
#     (1/n) sum_ij |<comparison_i|target_j>|^2 = tr(F^dag F)/n,   F_ij = <comparison_i|target_j>,
# for n comparison states (dense or sparse vectors). 1 when the comparison states lie inside the
# target manifold, whatever basis either manifold is returned in.
function groundstate_manifold_fidelity(comparison_states::AbstractVector{<:AbstractVector},target_states::AbstractVector{<:AbstractVector})
    return sum(abs2(dot(c,t)) for c in comparison_states, t in target_states) / length(comparison_states)
end





























"fin"
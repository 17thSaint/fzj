#####################################################
#=

This file contains the figure-of-merit function used to optimize, with QuOCS dCRAB,
the two sequential ramps used in tevo-daily-things.jl to connect a real-space
corner-pinned product state to the isotropic-hopping FCI ground-state manifold:

    pinned state --(ty: 0 -> end_ty)--> intermediate state --(tx: 0 -> end_tx)--> final state

Depends on:
    control-functions.jl
    time-evolution.jl (pulse_ramp, run_ramp_stages)

=#
######################################################

# pulse-independent part of the figure of merit: the pinned starting state and the
# target manifold. Run once per optimization (the QuOCS FoM object caches the
# returned tuple) instead of re-diagonalizing at every function evaluation.
function setup_pinned_ramp(parameters_dictionary)

    Lx::Int = Int(parameters_dictionary["Lx"])
    Ly::Int = Int(parameters_dictionary["Ly"])
    N::Int = Int(parameters_dictionary["N"])
    speccount::Int = Int(parameters_dictionary["speccount"])

    common_params = Dict{String,Any}(
        "output_level"=>0,"Lx"=>Lx,"Ly"=>Ly,"N"=>N,"lr"=>parameters_dictionary["lr"],
        "if_periodic_x"=>parameters_dictionary["if_periodic_x"],"if_periodic_y"=>parameters_dictionary["if_periodic_y"],
        "hopping_anisotropy"=>1.0,"interaction_strength"=>parameters_dictionary["interaction_strength"],
        "filling"=>0.5,"nev"=>speccount,"if_find_data"=>false,"if_save_data"=>false,
    )

    # starting_config arrives from Python as a flat [col1,row1,col2,row2,...] list,
    # one (col,row) pair per particle -- see pinnedRamp.starting_config
    flat_config = collect(parameters_dictionary["starting_config"])
    starting_config = [(Int(flat_config[2i-1]),Int(flat_config[2i])) for i in 1:N]

    # pinned starting state: particles confined to starting_config sites, no hopping
    states_starting,_,lattice_params,hamilt_params = position_state(starting_config,copy(common_params); output_level=0)
    hamilt_params["tx"] = 0.0

    # target isotropic-hopping ground-state manifold
    pdict_ending = merge(common_params,Dict("tx"=>parameters_dictionary["end_tx"],"ty"=>parameters_dictionary["end_ty"]))
    states_ending,_,_,_,_,_,_ = run_normal_ed(pdict_ending; output_level=0)

    starting_states = [Vector{ComplexF64}(states_starting[i]) for i in 1:speccount]
    target_states = [Vector{ComplexF64}(states_ending[i]) for i in 1:speccount]

    return (starting_states,target_states,lattice_params,hamilt_params)
end

function compute_fidelity_pinned_ramp(pulses,parameters_dictionary,setup)

    starting_states,target_states,lattice_params,hamilt_params = setup

    stages = [
        ("ty",collect(pulses[1]),parameters_dictionary["ramptime_firstramp"]),
        ("tx",collect(pulses[2]),parameters_dictionary["ramptime_secondramp"]),
    ]
    time_running_args = (nev=length(starting_states),output_level=0,if_instant_gs=false,if_save_data=false)
    final_states = run_ramp_stages(starting_states,stages,lattice_params,hamilt_params,parameters_dictionary["dt"]; time_running_args...)

    return real(groundstate_manifold_fidelity(final_states,target_states))
end

compute_fidelity_pinned_ramp(pulses,parameters_dictionary) = compute_fidelity_pinned_ramp(pulses,parameters_dictionary,setup_pinned_ramp(parameters_dictionary))

"fin"

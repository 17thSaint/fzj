#####################################################
#=

This file contains the figure-of-merit function used to optimize, with QuOCS dCRAB,
the interaction strength ramp used in tevo-daily-things.jl to connect the strongly
interacting ULR ground-state manifold to the FCI ground-state manifold:

    ULR manifold (U = intstren_start) --(interaction_strength ramp)--> FCI manifold (U = intstren_end)

Depends on:
    control-functions.jl (optimization_ed_params, groundstate_manifold_fidelity)
    time-evolution.jl (pulse_ramp, run_ramp_stages)

=#
######################################################

# pulse-independent part of the figure of merit: the endpoint ground-state manifolds.
# Run once per optimization (the QuOCS FoM object caches the returned tuple) instead
# of re-diagonalizing at every function evaluation.
function setup_intstren_ramp(parameters_dictionary)

    speccount::Int = Int(parameters_dictionary["speccount"])

    # strongly interacting ULR starting manifold
    pdict_starting = optimization_ed_params(parameters_dictionary,"interaction_strength"=>parameters_dictionary["intstren_start"])
    states_starting,_,_,_,_,lattice_params,hamilt_params = run_normal_ed(pdict_starting; output_level=0)

    # target FCI ground-state manifold
    pdict_ending = optimization_ed_params(parameters_dictionary,"interaction_strength"=>parameters_dictionary["intstren_end"])
    states_ending,_,_,_,_,_,_ = run_normal_ed(pdict_ending; output_level=0)

    starting_states = [Vector{ComplexF64}(states_starting[i]) for i in 1:speccount]
    target_states = [Vector{ComplexF64}(states_ending[i]) for i in 1:speccount]

    return (starting_states,target_states,lattice_params,hamilt_params)
end

function compute_fidelity_intstren_ramp(pulses,parameters_dictionary,setup)

    starting_states,target_states,lattice_params,hamilt_params = setup

    stages = [
        ("interaction_strength",collect(pulses[1]),parameters_dictionary["ramptime"]),
    ]
    time_running_args = (nev=length(starting_states),output_level=0,if_instant_gs=false,if_save_data=false)
    final_states = run_ramp_stages(starting_states,stages,lattice_params,hamilt_params,parameters_dictionary["dt"]; time_running_args...)

    return real(groundstate_manifold_fidelity(final_states,target_states))
end

"fin"

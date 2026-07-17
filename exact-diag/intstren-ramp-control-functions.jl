#####################################################
#=

This file contains the figure-of-merit function used to optimize, with QuOCS dCRAB,
the interaction strength ramp used in tevo-daily-things.jl to connect the strongly
interacting ULR ground-state manifold to the FCI ground-state manifold:

    ULR manifold (U = intstren_start) --(interaction_strength ramp)--> FCI manifold (U = intstren_end)

Depends on:
    control-functions.jl
    time-evolution.jl (pulse_ramp, run_ramp_stages)

=#
######################################################

# pulse-independent part of the figure of merit: the endpoint ground-state manifolds.
# Run once per optimization (the QuOCS FoM object caches the returned tuple) instead
# of re-diagonalizing at every function evaluation.
function setup_intstren_ramp(parameters_dictionary)

    Lx::Int = Int(parameters_dictionary["Lx"])
    Ly::Int = Int(parameters_dictionary["Ly"])
    N::Int = Int(parameters_dictionary["N"])
    speccount::Int = Int(parameters_dictionary["speccount"])

    common_params = Dict{String,Any}(
        "output_level"=>0,"Lx"=>Lx,"Ly"=>Ly,"N"=>N,"lr"=>parameters_dictionary["lr"],
        "if_periodic_x"=>parameters_dictionary["if_periodic_x"],"if_periodic_y"=>parameters_dictionary["if_periodic_y"],
        "hopping_anisotropy"=>1.0,"filling"=>0.5,"nev"=>speccount,"if_find_data"=>false,"if_save_data"=>false,
    )

    # forward the interaction profile settings so non-flat (exp / rydberg / dd / gaussian)
    # scalings reach get_normal_model_params_ed instead of silently falling back to flat
    for k in ("scaling_type","corr_length","sigma","blockade_radius","magnetic_spacing")
        haskey(parameters_dictionary,k) && (common_params[k] = parameters_dictionary[k])
    end

    # strongly interacting ULR starting manifold
    pdict_starting = merge(common_params,Dict("interaction_strength"=>parameters_dictionary["intstren_start"]))
    states_starting,_,_,_,_,lattice_params,hamilt_params = run_normal_ed(pdict_starting; output_level=0)

    # target FCI ground-state manifold
    pdict_ending = merge(common_params,Dict("interaction_strength"=>parameters_dictionary["intstren_end"]))
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

compute_fidelity_intstren_ramp(pulses,parameters_dictionary) = compute_fidelity_intstren_ramp(pulses,parameters_dictionary,setup_intstren_ramp(parameters_dictionary))

"fin"

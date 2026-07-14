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

function compute_fidelity_intstren_ramp(pulses,parameters_dictionary)

    Lx::Int = Int(parameters_dictionary["Lx"])
    Ly::Int = Int(parameters_dictionary["Ly"])
    N::Int = Int(parameters_dictionary["N"])
    speccount::Int = Int(parameters_dictionary["speccount"])
    dt = parameters_dictionary["dt"]

    common_params = Dict{String,Any}(
        "output_level"=>0,"Lx"=>Lx,"Ly"=>Ly,"N"=>N,"lr"=>parameters_dictionary["lr"],
        "if_periodic_x"=>parameters_dictionary["if_periodic_x"],"if_periodic_y"=>parameters_dictionary["if_periodic_y"],
        "hopping_anisotropy"=>1.0,"filling"=>0.5,"nev"=>speccount,"if_find_data"=>false,"if_save_data"=>false,
    )

    # strongly interacting ULR starting manifold
    pdict_starting = merge(common_params,Dict("interaction_strength"=>parameters_dictionary["intstren_start"]))
    states_starting,_,_,_,_,lattice_params,hamilt_params = run_normal_ed(pdict_starting; output_level=0)

    # target FCI ground-state manifold
    pdict_ending = merge(common_params,Dict("interaction_strength"=>parameters_dictionary["intstren_end"]))
    states_ending,_,_,_,_,_,_ = run_normal_ed(pdict_ending; output_level=0)

    starting_states = [Vector{ComplexF64}(states_starting[i]) for i in 1:speccount]
    stages = [
        ("interaction_strength",collect(pulses[1]),parameters_dictionary["ramptime"]),
    ]
    time_running_args = (nev=speccount,output_level=0,if_instant_gs=false,if_save_data=false)
    final_states = run_ramp_stages(starting_states,stages,lattice_params,hamilt_params,dt; time_running_args...)

    return real(groundstate_manifold_fidelity(final_states,[Vector{ComplexF64}(s) for s in states_ending[1:speccount]]))
end

"fin"

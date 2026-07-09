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

function compute_fidelity_pinned_ramp(pulses,parameters_dictionary)

    Lx::Int = Int(parameters_dictionary["Lx"])
    Ly::Int = Int(parameters_dictionary["Ly"])
    N::Int = Int(parameters_dictionary["N"])
    speccount::Int = Int(parameters_dictionary["speccount"])
    dt = parameters_dictionary["dt"]

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
    stages = [
        ("ty",collect(pulses[1]),parameters_dictionary["ramptime_firstramp"]),
        ("tx",collect(pulses[2]),parameters_dictionary["ramptime_secondramp"]),
    ]
    time_running_args = (nev=speccount,output_level=0,if_instant_gs=false,if_save_data=false)
    final_states = run_ramp_stages(starting_states,stages,lattice_params,hamilt_params,dt; time_running_args...)

    return real(groundstate_manifold_fidelity(final_states,states_ending[1:speccount]))
end

"fin"

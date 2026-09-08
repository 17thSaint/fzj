#####################################################
#=

This file contains any random functions written to do one-off tasks

Depends on:
    execute-ed.jl

=#
######################################################

include("execute-ed.jl")
include("time-evolution.jl")
include("control-functions.jl")
include("../other-funcs/basic-2d-plottings.jl")
include("plottings.jl")

using NPZ  # for reading QuOCS best_controls.npz files

## ramp from strongly interacting state to FCI

#=if_all::Bool = true
# model parameters
if false || if_all
    lx,ly,n = 4,4,2
    
    intstren = 10.0

    tmax_global = 10.0
    dt_global = 0.0001
    dataloc = "tevo-daily-things-data/"
end

#= plot the ULR for starting and ending states
corrlengths = range(0.1,50.0,length=12)
for xi in corrlengths
    us_starting = long_range_scaling(ly-1,ly,intstren; corr_length=xi,scaling="exp")
    plot(0:length(us_starting)-1,us_starting,label="$(round(xi, digits=1))")
end
xlabel("y distance")
ylabel("Interaction strength")
title("ULR for starting and ending states")
legend()=#

# define starting state
if false || if_all

    start_tx = 1e-3

    params_dict_starting = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("tx",start_tx),("ty",1.0),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",5),("if_find_data",false),("if_save_data",false)])

    states_i,nrgs_i,rhos_i,filepath_i,if_found_i,lattice_params_i,hamilt_params_i = run_normal_ed(params_dict_starting; output_level=0)
    println("Found starting state")

    psi0_1 = states_i[1]
    psi0_2 = states_i[2]

end

# define ending state
if false || if_all

    end_tx = 1.0

    params_dict_ending = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("ty",1.0),("tx",end_tx),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",5),("if_find_data",false),("if_save_data",false)])

    states_f,nrgs_f,rhos_f,filepath_f,if_found_f,lattice_params_f,hamilt_params_f = run_normal_ed(params_dict_ending; output_level=0)
    println("Found ending state")

    psif_1 = states_f[1]
    psif_2 = states_f[2]

end

# time evolution with linear ramp of interaction length for different ramp times
if false || if_all
    speccount = 2
    time_running_args = (nev=speccount,output_level=1,if_instant_gs=false,if_save_data=false,dataloc=dataloc,)

    #quench_fidelity = groundstate_manifold_fidelity(states_f[1:speccount],states_i[1:speccount])

    ramptimes = 10 .^ range(-0.9,1.0,length=5)
    for ramptime in ramptimes
        tmax_global = ramptime + 0.5
        tevo_params = Dict([ ("tx",(linear_ramp,start_tx,end_tx,ramptime)),("dt",dt_global),("tmax",tmax_global) ])
        tevo_gs,tevo_dict,intspec,saving_args = run_timeevo([psi0_1,psi0_2],tevo_params,lattice_params_i,hamilt_params_i; time_running_args...)

    #end

    # calculate final fidelity with target state
    #if true || if_all
        final_gs_manifold = [tevo_gs[1][:,end-1],tevo_gs[2][:,end-1]]
        final_fidelity = groundstate_manifold_fidelity(states_f[1:speccount],final_gs_manifold)
        println("Final fidelity for ramp time $(ramptime): $(final_fidelity)")
        scatter(ramptime,final_fidelity,c="b")
    end
    xlabel("Ramp time")
    ylabel("Fidelity with target manifold")
    title("Fidelity vs ramp time $(lx)x$(ly) N=$(n) U=$(intstren) ramp tx")
    xscale("log")

end=#

#= check finite size scaling of manifold overlap
# not much interesting result and can't go larger because TTN precision isn't high enough
if false
    lattices = [(6,3,3),(8,4,4),(10,5,5)]
    overlaps = Float64[]
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    for (lx,ly,n) in lattices
        intstren = 300.0
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren)])
        all_files_ulr = find_data_file(pdict,"ed",dataloc; file_type="jld2")
        filter!(f -> !occursin("twist_angle",f), all_files_ulr)
        display(all_files_ulr)
        length(all_files_ulr) != 1 && error("Expected exactly one file for ULR state, found $(length(all_files_ulr))")
        d,m = read_data(joinpath(dataloc,all_files_ulr[1]); output_level=0)

        pdict["interaction_strength"] = 0.0
        all_files_laugh = find_data_file(pdict,"ed",dataloc; file_type="jld2")
        filter!(f -> !occursin("twist_angle",f), all_files_laugh)
        length(all_files_laugh) != 1 && error("Expected exactly one file for Laughlin state, found $(length(all_files_laugh))")
        d_laugh,m_laugh = read_data(joinpath(dataloc,all_files_laugh[1]); output_level=0)

        manifold_overlap = groundstate_manifold_fidelity(d["state"][1:2],d_laugh["state"][1:2])
        #println("Manifold overlap between ULR and Laughlin states: $(manifold_overlap)")
        push!(overlaps, manifold_overlap)
    end
    scatter([lx for (lx,ly,n) in lattices], overlaps, c="b")
    xlabel("Lattice size (Lx)")
    ylabel("Manifold overlap between ULR and Laughlin states")
    title("Finite size scaling of manifold overlap")


end=#

#= check manifold overlap for periodic potential states, use this to benchmark dt
if false
    lx,ly,n = 6,3,3
    intstren = 300.0
    ppstren_start = 20.0
    ppstren_end = 0.0
    
    pdict_start = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("periodic_potential_strength",ppstren_start),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])
    states_start,nrgs_start,rhos_start,filepath_start,if_found_start,lattice_params_start,hamilt_params_start = run_normal_ed(pdict_start; output_level=1)

    pdict_end = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("periodic_potential_strength",ppstren_end),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])
    states_end,nrgs_end,rhos_end,filepath_end,if_found_end,lattice_params_end,hamilt_params_end = run_normal_ed(pdict_end; output_level=1)

    manifold_overlap = groundstate_manifold_fidelity(states_start[1:2],states_end[1:2])
    println("Manifold overlap between periodic potential states: $(manifold_overlap)")
end=#

#= benchmark dt: compare each dt against finest-dt reference run for periodic potential ramp
if false
    speccount = 2
    ramptime = 0.5
    tmax_global = ramptime + 0.5
    time_running_args = (nev=speccount, output_level=0, if_instant_gs=false, if_save_data=false)

    dts = sort(10 .^ range(log10(0.001), log10(0.05), length=8))  # 0.001 → 0.05, predicted cutoff ~0.003

    all_final_manifolds = []
    for dt_global in dts
        tevo_params = Dict([
            ("periodic_potential_strength", (linear_ramp, ppstren_start, ppstren_end, ramptime)),
            ("dt", dt_global),
            ("tmax", tmax_global)
        ])
        tevo_gs, _, _, _ = run_timeevo(
            [states_start[1], states_start[2]],
            tevo_params, lattice_params_start, hamilt_params_start;
            time_running_args...
        )
        push!(all_final_manifolds, [Vector(tevo_gs[1][:, end-1]), Vector(tevo_gs[2][:, end-1])])
        println("Finished dt=$(dt_global)")
    end

    ref_manifold = all_final_manifolds[1]  # finest dt is the reference
    figure()
    for (i, dt_global) in enumerate(dts)
        fid = real(groundstate_manifold_fidelity(ref_manifold, all_final_manifolds[i]))
        println("dt=$(dt_global): fidelity vs reference = $(fid)")
        scatter(dt_global, fid, c="b")
    end
    xscale("log")
    xlabel("Time step dt")
    ylabel("Fidelity vs finest-dt reference")
    title("dt convergence $(lx)x$(ly) N=$(n) U=$(intstren) ppstren $(ppstren_start)→$(ppstren_end) T=$(ramptime)")
end=#

#= benchmark dt: tx ramp
if false
    lx, ly, n = 6, 3, 3
    intstren = 300.0
    tx_start = 1e-3
    tx_end = 1.0

    pdict_tx_start = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("tx",tx_start),("ty",1.0),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("filling",0.5),("nev",5),("if_find_data",false),("if_save_data",false)])
    states_tx_start,_,_,_,_,lattice_params_tx,hamilt_params_tx = run_normal_ed(pdict_tx_start; output_level=0)

    pdict_tx_end = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("tx",tx_end),("ty",1.0),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("filling",0.5),("nev",5),("if_find_data",false),("if_save_data",false)])
    states_tx_end,_,_,_,_,_,_ = run_normal_ed(pdict_tx_end; output_level=0)

    speccount = 2
    ramptime = 0.5
    tmax_global = ramptime + 0.5
    dts = sort(10 .^ range(log10(0.001), log10(0.05), length=8))
    time_running_args = (nev=speccount, output_level=0, if_instant_gs=false, if_save_data=false)

    all_final_manifolds_tx = []
    for dt_global in dts
        tevo_params = Dict([
            ("tx", (linear_ramp, tx_start, tx_end, ramptime)),
            ("dt", dt_global),
            ("tmax", tmax_global)
        ])
        tevo_gs, _, _, _ = run_timeevo(
            [states_tx_start[1], states_tx_start[2]],
            tevo_params, lattice_params_tx, hamilt_params_tx;
            time_running_args...
        )
        push!(all_final_manifolds_tx, [Vector(tevo_gs[1][:, end-1]), Vector(tevo_gs[2][:, end-1])])
        println("tx ramp: Finished dt=$(dt_global)")
    end

    ref_manifold_tx = all_final_manifolds_tx[1]
    figure()
    for (i, dt_global) in enumerate(dts)
        fid = real(groundstate_manifold_fidelity(ref_manifold_tx, all_final_manifolds_tx[i]))
        println("tx ramp: dt=$(dt_global): fidelity vs reference = $(fid)")
        scatter(dt_global, fid, c="b")
    end
    xscale("log")
    xlabel("Time step dt")
    ylabel("Fidelity vs finest-dt reference")
    title("dt convergence $(lx)x$(ly) N=$(n) U=$(intstren) tx $(tx_start)→$(tx_end) T=$(ramptime)")
end=#

#= benchmark dt: dipole-dipole magnetic_spacing ramp
# for 6x3 there is a breakdown in the degenerate manifold description
if false
    lx, ly, n = 6, 3, 3
    intstren = 300.0
    ms_start = 0.5
    ms_end = 2.0

    pdict_dd_start = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("scaling_type","dd"),("magnetic_spacing",ms_start),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("filling",0.5),("nev",5),("if_find_data",false),("if_save_data",false)])
    states_dd_start,_,_,_,_,lattice_params_dd,hamilt_params_dd = run_normal_ed(pdict_dd_start; output_level=0)

    pdict_dd_end = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("scaling_type","dd"),("magnetic_spacing",ms_end),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("filling",0.5),("nev",5),("if_find_data",false),("if_save_data",false)])
    states_dd_end,_,_,_,_,_,_ = run_normal_ed(pdict_dd_end; output_level=0)

    speccount = 2
    ramptime = 0.5
    tmax_global = ramptime + 0.5
    dts = sort(10 .^ range(log10(0.001), log10(0.05), length=8))
    time_running_args = (nev=speccount, output_level=0, if_instant_gs=false, if_save_data=false)

    all_final_manifolds_dd = []
    for dt_global in dts
        tevo_params = Dict([
            ("magnetic_spacing", (linear_ramp, ms_start, ms_end, ramptime)),
            ("dt", dt_global),
            ("tmax", tmax_global)
        ])
        tevo_gs, _, _, _ = run_timeevo(
            [states_dd_start[1], states_dd_start[2]],
            tevo_params, lattice_params_dd, hamilt_params_dd;
            time_running_args...
        )
        push!(all_final_manifolds_dd, [Vector(tevo_gs[1][:, end-1]), Vector(tevo_gs[2][:, end-1])])
        println("dd ramp: Finished dt=$(dt_global)")
    end

    ref_manifold_dd = all_final_manifolds_dd[1]
    figure()
    for (i, dt_global) in enumerate(dts)
        fid = real(groundstate_manifold_fidelity(ref_manifold_dd, all_final_manifolds_dd[i]))
        println("dd ramp: dt=$(dt_global): fidelity vs reference = $(fid)")
        scatter(dt_global, fid, c="r")
    end
    xscale("log")
    xlabel("Time step dt")
    ylabel("Fidelity vs finest-dt reference")
    title("dt convergence $(lx)x$(ly) N=$(n) U=$(intstren) dd ms $(ms_start)→$(ms_end) T=$(ramptime)")
end=#



### Initialize single column full and then ramp tx into FCI state
# Strategy: start particles pinned to real-space sites (ty=0, no hopping),
#= then adiabatically ramp ty → 1 to connect to the FCI ground state manifold.
if false

    if_all::Bool = false

    # model parameters
    if false || if_all
        lx,ly,n = 4,4,2

        intstren = 0.0  # non-interacting: topology alone drives the FCI

        # pre-compute full Fock basis and cache to avoid rebuilding it in each sub-block
        lattice_params::Dict{String,Any} = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true)])
        full_basis = n_particle_basis(lattice_params; output_level=0,dataloc=get_folder_location("cluster-data/exact-diag"))
        lattice_params["full_basis"] = full_basis
    end

    # define starting state: particles pinned to specific real-space sites, no hopping
    if false || if_all
        # each tuple is a (column, row) site index for one of the n particles
        starting_config = [(1,1),(1,2)]

        pdict_starting = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])

        states_starting, nrgs_starting, lattice_params_starting, hamilt_params_starting = position_state(starting_config, pdict_starting; output_level=0)
        hamilt_params_starting["tx"] = 0.0

        #occs_starting = get_occupancy(states_starting[1], lattice_params_starting; plot_title="Starting state occupancy")
    end

    # define ending state: isotropic hopping target used for fidelity comparison
    if false || if_all
        end_tx = 1.0
        end_ty = 1.0

        pdict_ending = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("tx",end_tx),("ty",end_ty),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])

        states_ending, nrgs_ending, _, _, _, lattice_params_ending, hamilt_params_ending = run_normal_ed(pdict_ending; output_level=0)
    end

    firstramp_times = range(0.01, 0.1, length=11)
    final_fidelities = zeros(length(firstramp_times))
    for (idx,ramptime_firstramp) in enumerate(firstramp_times)

        # ramp ty from 0 to 1; tx stays at its default from hamilt_params_starting
        if true || if_all
            speccount_firstramp = 3
            #ramptime_firstramp = 0.5
            tmax_global_firstramp = ramptime_firstramp + 0.0  # extra hold time after ramp end to check convergence
            time_running_args_firstramp = (nev=speccount_firstramp, output_level=1, if_instant_gs=true, if_save_data=false, dataloc="tevo-daily-things-data/")

            tevo_params_firstramp = Dict([ ("ty",(linear_ramp,0.0,end_ty,ramptime_firstramp)),("tmax",tmax_global_firstramp) ])
            tevo_data_firstramp, tevo_dict_firstramp, instdata_firstramp, saving_args_firstramp = run_timeevo([states_starting[1],states_starting[2],states_starting[3]],tevo_params_firstramp,lattice_params_starting,hamilt_params_starting; time_running_args_firstramp...)

            # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
            #occs_midpoint = get_occupancy(tevo_gs_firstramp[1][:,end-1], lattice_params_starting; plot_title="Midpoint state occupancy")
        end

        # ramp tx from 0 to 1; ty stays at its default from hamilt_params_starting
        if true || if_all
            speccount_secondramp = 3
            ramptime_secondramp = ramptime_firstramp
            tmax_global_secondramp = ramptime_secondramp + 0.0  # extra hold time after ramp
            time_running_args_secondramp = (nev=speccount_secondramp, output_level=1, if_instant_gs=true, if_save_data=false, dataloc="tevo-daily-things-data/")

            initial_states = [Vector{ComplexF64}(tevo_data_firstramp[1][1][:,end-1]),Vector{ComplexF64}(tevo_data_firstramp[1][2][:,end-1]),Vector{ComplexF64}(tevo_data_firstramp[1][3][:,end-1])]
            tevo_params_secondramp = Dict([ ("tx",(linear_ramp,0.0,end_tx,ramptime_secondramp)),("tmax",tmax_global_secondramp) ])
            tevo_data_secondramp, tevo_dict_secondramp, instdata_secondramp, saving_args_secondramp = run_timeevo(initial_states,tevo_params_secondramp,lattice_params_starting,hamilt_params_starting; time_running_args_secondramp...)

        end

        # displaying and plotting stuff
        if true || if_all
            final_states = [Vector{ComplexF64}(tevo_data_secondramp[1][1][:,end-1]),Vector{ComplexF64}(tevo_data_secondramp[1][2][:,end-1]),Vector{ComplexF64}(tevo_data_secondramp[1][3][:,end-1])]
            #=final_nrgs = [real(adjoint(wavefunc) * hamilt_params_ending["H"] * wavefunc) for wavefunc in final_states]
            display(final_nrgs)

            overlap_matrix = zeros(Float64, speccount_secondramp, 2)
            for i in 1:speccount_secondramp
                for j in 1:2
                    overlap_matrix[i,j] = abs2(adjoint(final_states[i]) * states_ending[j])
                end
            end
            display(overlap_matrix)=#

            #=times_firstramp = range(0.0, tmax_global_firstramp, length=length(instdata_firstramp[2]["1"]))
            times_secondramp = range(0.0, tmax_global_secondramp, length=length(instdata_secondramp[2]["1"]))

            cols = ["b","g","r"]
            for i in 1:3
                scatter(times_firstramp,instdata_firstramp[2][string(i)],c=cols[i],label="E$(i)")
                scatter(times_firstramp,tevo_data_firstramp[2][i][1:end-1],c="k",label="E$(i)",marker="x")
                scatter(times_secondramp .+ tmax_global_firstramp,instdata_secondramp[2][string(i)],c=cols[i])
                scatter(times_secondramp .+ tmax_global_firstramp,tevo_data_secondramp[2][i][1:end-1],c="k",marker="x")
            end
            legend()
            xlabel("Time")
            ylabel("Energy")
            title("Energy vs time for two ramps $(lx)x$(ly) N=$(n) ramptime $(ramptime_firstramp) ty and tx")=#

            final_fidelity = groundstate_manifold_fidelity(final_states[1:2],states_ending[1:2])

            final_fidelities[idx] = final_fidelity
            println("Final fidelity for ramp time $(ramptime_firstramp): $(final_fidelity)")
            scatter(ramptime_firstramp,final_fidelity,c="b")
            xlabel("Ramp time")
            ylabel("Fidelity with target manifold")
            title("Fidelity vs ramp time $(lx)x$(ly) N=$(n) U=$(intstren) ramp tx and ty")
            xscale("log")


            # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
            #occs_final = get_occupancy(tevo_gs_secondramp[1][:,end-1], lattice_params_starting; plot_title="Final Fidelity = $(round(final_fidelity,digits=6))")
        end

    end
end=#


### Time evolution with the QuOCS-optimized pulses from optimal-control/config_pinnedRamp.py
# All parameters must match the ones the optimization ran with (see config_pinnedRamp.py):
#= 4x4 N=2 U=0 pbc, particles pinned at [(1,1),(1,2)], ty then tx ramped 0 -> 1 over 0.5 each
if false

    if_all::Bool = false

    # model parameters and starting/ending states
    if false || if_all
        lx,ly,n = 4,4,2
        intstren = 0.0
        end_tx, end_ty = 1.0, 1.0
        speccount_quocs = 2  # optimization used the 2-state groundstate manifold
        speccount_energy = 3  # low-lying states tracked in the energy-vs-time section below

        starting_config = [(1,1),(1,2)]
        pdict_quocs = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",speccount_energy),("if_find_data",false),("if_save_data",false)])

        states_starting, nrgs_starting, lattice_params_starting, hamilt_params_starting = position_state(starting_config, copy(pdict_quocs); output_level=0)
        hamilt_params_starting["tx"] = 0.0

        pdict_ending = merge(pdict_quocs, Dict("tx"=>end_tx,"ty"=>end_ty))
        states_ending,_,_,_,_,_,_ = run_normal_ed(pdict_ending; output_level=0)
    end

    # load the optimized pulses and run the two-stage time evolution
    if false || if_all
        quocs_folder = "../optimal-control/QuOCS_Results/20260709_162520_pinnedRamp_dCRAB"
        controls_file = filter(f -> endswith(f,"best_controls.npz"), readdir(quocs_folder))[1]
        # only read the numeric arrays: NPZ.jl cannot parse the numpy unicode-string arrays
        # (pulse_names etc.) that QuOCS also stores in the file
        best_controls = npzread(joinpath(quocs_folder,controls_file),["tyRamp","txRamp","time_grid_for_tyRamp","time_grid_for_txRamp"])

        # pulses are sampled on the RK4 half-step grid (spacing dt/2), so dt must match the
        # value used in config_pinnedRamp.py: pulse length = ceil(2*ramptime/dt) + 1
        dt_quocs = 0.005
        ramptime_ty = best_controls["time_grid_for_tyRamp"][end]
        ramptime_tx = best_controls["time_grid_for_txRamp"][end]

        starting_states_quocs = [Vector{ComplexF64}(states_starting[i]) for i in 1:speccount_quocs]
        stages_quocs = [
            ("ty",best_controls["tyRamp"],ramptime_ty),
            ("tx",best_controls["txRamp"],ramptime_tx),
        ]
        time_running_args_quocs = (nev=speccount_quocs,output_level=0,if_instant_gs=false,if_save_data=false)
        final_states_quocs = run_ramp_stages(starting_states_quocs,stages_quocs,lattice_params_starting,hamilt_params_starting,dt_quocs; time_running_args_quocs...)

        fidelity_quocs = real(groundstate_manifold_fidelity(final_states_quocs,[Vector{ComplexF64}(s) for s in states_ending[1:speccount_quocs]]))
        println("Fidelity with target manifold using QuOCS pulses: $(fidelity_quocs)")
    end

    # instantaneous groundstate energies vs time-evolved energies for a few low-lying states
    # along the QuOCS pulses, plotted like the linear-ramp version in the commented block above
    if false || if_all
        # work on a copy: run_timeevo's timeham writes each ramp's current value back into the
        # dict, which would leave tx=1.0 in hamilt_params_starting for any later section
        hamilt_params_energy = copy(hamilt_params_starting)
        hamilt_params_energy["tx"] = 0.0

        time_running_args_energy = (nev=speccount_energy, output_level=1, if_instant_gs=true, if_save_data=false, dataloc="tevo-daily-things-data/")

        starting_states_energy = [Vector{ComplexF64}(states_starting[i]) for i in 1:speccount_energy]
        tevo_params_tyramp = Dict([ ("ty",(pulse_ramp,ramptime_ty,best_controls["tyRamp"])),("tmax",ramptime_ty),("dt",dt_quocs) ])
        tevo_data_tyramp, tevo_dict_tyramp, instdata_tyramp, saving_args_tyramp = run_timeevo(starting_states_energy,tevo_params_tyramp,lattice_params_starting,hamilt_params_energy; time_running_args_energy...)

        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        midpoint_states_energy = [Vector{ComplexF64}(tevo_data_tyramp[1][i][:,end-1]) for i in 1:speccount_energy]
        tevo_params_txramp = Dict([ ("tx",(pulse_ramp,ramptime_tx,best_controls["txRamp"])),("tmax",ramptime_tx),("dt",dt_quocs) ])
        tevo_data_txramp, tevo_dict_txramp, instdata_txramp, saving_args_txramp = run_timeevo(midpoint_states_energy,tevo_params_txramp,lattice_params_starting,hamilt_params_energy; time_running_args_energy...)

    end

    if true || if_all
        times_tyramp = range(0.0, ramptime_ty, length=length(instdata_tyramp[2]["1"]))
        times_txramp = range(0.0, ramptime_tx, length=length(instdata_txramp[2]["1"]))

        figure()
        cols = ["b","g","r"]
        for i in 1:speccount_energy
            plot(times_tyramp,instdata_tyramp[2][string(i)],c=cols[i],"-p",label="E$(i)")
            plot(times_tyramp,tevo_data_tyramp[2][i][1:end-1],c="k",marker="x")
            plot(times_txramp .+ ramptime_ty,instdata_txramp[2][string(i)],"-p",c=cols[i])
            plot(times_txramp .+ ramptime_ty,tevo_data_txramp[2][i][1:end-1],c="k",marker="x")
        end
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time for QuOCS pulses $(lx)x$(ly) N=$(n) ramptimes $(ramptime_ty) ty and $(ramptime_tx) tx")
    end
    
end=#


### Ramp from strongly interacting state to FCI with linear ramp of interaction strength
#= Strategy: start with strongly interacting state (ULR) and ramp interaction strength to 0.0 to connect to the FCI ground state manifold.
if false
    if_all::Bool = true

    # define starting state: strongly interacting ULR state
    if false || if_all
        lx,ly,n = 4,4,2
        intstren_start = 10.0

        pdict_starting = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren_start),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])

        states_starting, nrgs_starting,_,_,_, lattice_params_starting, hamilt_params_starting = run_normal_ed(pdict_starting; output_level=0)
    end

    # define ending state: FCI
    if false || if_all
        intstren_end = 0.0

        pdict_ending = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren_end),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])

        states_ending, nrgs_ending,_,_,_, lattice_params_ending, hamilt_params_ending = run_normal_ed(pdict_ending; output_level=0)
    end

    # define reference fidelity
    if false || if_all
        speccount = 2
        reference_fidelity = groundstate_manifold_fidelity(states_ending[1:speccount],states_starting[1:speccount])
        println("Reference fidelity between ULR and FCI states: $(reference_fidelity)")
    end

    # time evolution with linear ramp of interaction strength for various ramp times
    if false || if_all
        speccount = 2
        dataloc = get_folder_location("cluster-data/exact-diag/time-evo")
        time_running_args = (nev=speccount,output_level=1,if_instant_gs=false,if_save_data=false,dataloc=dataloc,)

        ramptimes = 10 .^ range(-2.0,1.5,length=21)
        for ramptime in ramptimes
            tmax_global = ramptime
            tevo_params = Dict([ ("interaction_strength",(linear_ramp,intstren_start,intstren_end,ramptime)),("tmax",tmax_global) ])
            tevo_data,tevo_dict,_,saving_args = run_timeevo([states_starting[1],states_starting[2]],tevo_params,lattice_params_starting,hamilt_params_starting; time_running_args...)

            # calculate final fidelity with target state
            final_gs_manifold = [tevo_data[1][1][:,end-1],tevo_data[1][2][:,end-1]]
            final_fidelity = groundstate_manifold_fidelity(final_gs_manifold,states_ending[1:speccount])
            println("Final fidelity for ramp time $(ramptime): $(final_fidelity)")
            scatter(ramptime,final_fidelity,c="b")
        end
        plot([0.0,10.0],[reference_fidelity,reference_fidelity],"--",c="r")
        xlabel("Ramp time")
        ylabel("Fidelity with target manifold")
        title("Fidelity vs ramp time $(lx)x$(ly) N=$(n) U=$(intstren_start)→$(intstren_end) ramp interaction strength")
        xscale("log")

    end

    #= time evolution with linear ramp looking at instantaneous energies
    if false || if_all
        speccount = 3
        dataloc = get_folder_location("cluster-data/exact-diag/time-evo")
        time_running_args = (nev=speccount,output_level=1,if_instant_gs=true,if_save_data=false,dataloc=dataloc,)

        ramptime = 10.0
        tmax_global = ramptime
        tevo_params = Dict([ ("interaction_strength",(linear_ramp,intstren_start,intstren_end,ramptime)),("tmax",tmax_global) ])
        tevo_data,tevo_dict,instdata,saving_args = run_timeevo([states_starting[1],states_starting[2],states_starting[3]],tevo_params,lattice_params_starting,hamilt_params_starting; time_running_args...)

        times = range(0.0,tmax_global,length=length(instdata[2]["1"]))
        cols = ["b","g","r"]
        for i in 1:speccount
            plot(times,instdata[2][string(i)],c=cols[i],"-p",label="E$(i)")
            plot(times,tevo_data[2][i][1:end-1],c="k",marker="x")
        end
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time for interaction strength ramp $(lx)x$(ly) N=$(n) ramptime $(ramptime) U $(intstren_start)→$(intstren_end)")
    end=#
end=#


### Time evolution with the QuOCS-optimized interaction strength ramp from optimal-control/config_intstrenRamp.py
# All parameters must match the ones the optimization ran with (see config_intstrenRamp.py):
# 4x4 N=2 pbc, U ramped 10.0 -> 0.0 over ramptime 1.0, dt 0.005
#= Reference: linear ramp fidelity 0.8891, QuOCS-optimized fidelity 0.9101 (run 20260714_114412) but Claude seems to think more superiterations of DCRAB would improve it further
if false

    if_all::Bool = true

    # model parameters and starting/ending states
    if false || if_all
        lx,ly,n = 4,4,2
        intstren_start, intstren_end = 10.0, 0.0
        speccount_intquocs = 2

        pdict_intquocs = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren_start),("filling",0.5),("nev",3),("if_find_data",false),("if_save_data",false)])
        states_starting_intquocs,_,_,_,_,lattice_params_intquocs,hamilt_params_intquocs = run_normal_ed(pdict_intquocs; output_level=0)

        pdict_ending_intquocs = merge(pdict_intquocs,Dict("interaction_strength"=>intstren_end))
        states_ending_intquocs,_,_,_,_,_,_ = run_normal_ed(pdict_ending_intquocs; output_level=0)
    end

    # load the optimized pulse and run the time evolution with instantaneous energies
    if false || if_all
        quocs_folder = "../optimal-control/QuOCS_Results/20260714_114412_intstrenRamp_dCRAB"
        controls_file = filter(f -> endswith(f,"best_controls.npz"), readdir(quocs_folder))[1]
        # only read the numeric arrays: NPZ.jl cannot parse the numpy unicode-string arrays
        # (pulse_names etc.) that QuOCS also stores in the file
        best_controls = npzread(joinpath(quocs_folder,controls_file),["intstrenRamp","time_grid_for_intstrenRamp"])

        # pulse is sampled on the RK4 half-step grid (spacing dt/2), so dt must match the
        # value used in config_intstrenRamp.py: pulse length = ceil(2*ramptime/dt) + 1
        dt_intquocs = 0.005
        ramptime_intquocs = best_controls["time_grid_for_intstrenRamp"][end]

        # work on a copy: run_timeevo's timeham writes the ramp's current value back into
        # the dict, which would leave interaction_strength=0.0 for any later section
        hamilt_params_energy_intquocs = copy(hamilt_params_intquocs)

        time_running_args_intquocs = (nev=speccount_intquocs,output_level=1,if_instant_gs=true,if_save_data=false,dataloc="tevo-daily-things-data/")
        starting_states_intquocs = [Vector{ComplexF64}(states_starting_intquocs[i]) for i in 1:speccount_intquocs]
        tevo_params_intquocs = Dict([ ("interaction_strength",(pulse_ramp,ramptime_intquocs,best_controls["intstrenRamp"])),("tmax",ramptime_intquocs),("dt",dt_intquocs) ])
        tevo_data_intquocs,tevo_dict_intquocs,instdata_intquocs,saving_args_intquocs = run_timeevo(starting_states_intquocs,tevo_params_intquocs,lattice_params_intquocs,hamilt_params_energy_intquocs; time_running_args_intquocs...)

        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        final_manifold_intquocs = [Vector{ComplexF64}(tevo_data_intquocs[1][i][:,end-1]) for i in 1:speccount_intquocs]
        fidelity_intquocs = real(groundstate_manifold_fidelity(final_manifold_intquocs,[Vector{ComplexF64}(s) for s in states_ending_intquocs[1:speccount_intquocs]]))
        println("Fidelity with target manifold using QuOCS pulse: $(fidelity_intquocs)")
    end

    # plot the optimized pulse and the instantaneous vs time-evolved energies along it
    if true || if_all
        times_intquocs = range(0.0,ramptime_intquocs,length=length(instdata_intquocs[2]["1"]))

        figure()
        cols = ["b","g","r"]
        for i in 1:speccount_intquocs
            plot(times_intquocs,instdata_intquocs[2][string(i)],"-p",c=cols[i],label="E$(i)")
            plot(times_intquocs,tevo_data_intquocs[2][i][1:end-1],c="k",marker="x")
        end
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time for QuOCS interaction strength ramp $(lx)x$(ly) N=$(n) ramptime $(ramptime_intquocs) U $(intstren_start)→$(intstren_end)")

        figure()
        plot(best_controls["time_grid_for_intstrenRamp"],best_controls["intstrenRamp"])
        xlabel("Time")
        ylabel("Interaction strength")
        title("QuOCS optimized pulse, fidelity = $(round(fidelity_intquocs,digits=6))")
    end

end=#


### Ramp from strongly interacting state to FCI with linear ramp of interaction strength for 8x4 lattice
#= Strategy: start with strongly interacting state (ULR) and ramp interaction strength to 0.0 to connect to the FCI ground state manifold.
if false
    if_all::Bool = true

    # define starting state: strongly interacting ULR state
    if false || if_all
        lx,ly,n = 8,4,4
        intstren_start = 10.0

        pdict_starting = Dict([("output_level",1),("if_reading",true),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren_start),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])

        states_starting, nrgs_starting,_,_,_, lattice_params_starting, hamilt_params_starting = run_normal_ed(pdict_starting; output_level=0)
    end

    # define ending state: FCI
    if false || if_all
        intstren_end = 0.0

        pdict_ending = Dict([("output_level",1),("if_reading",true),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren_end),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])

        states_ending, nrgs_ending,_,_,_, lattice_params_ending, hamilt_params_ending = run_normal_ed(pdict_ending; output_level=0)
    end

    # define reference fidelity
    if false || if_all
        speccount = 2
        reference_fidelity = groundstate_manifold_fidelity(states_ending[1:speccount],states_starting[1:speccount])
        println("Reference fidelity between ULR and FCI states: $(reference_fidelity)")
    end

    # time evolution with linear ramp of interaction strength for various ramp times
    if true || if_all
        speccount = 2
        dataloc = get_folder_location("cluster-data/exact-diag/time-evo/ulr-flat")
        time_running_args = (nev=speccount,output_level=1,if_instant_gs=false,if_reading=true,if_save_data=true,dataloc=dataloc,)

        tevo_hamilt_params = copy(hamilt_params_starting)

        #ramptimes = 10 .^ range(-2.0,1.5,length=21)
        #for ramptime in ramptimes
        ramptime = 0.1
        tmax_global = ramptime
        tevo_params = Dict([ ("interaction_strength",(linear_ramp,intstren_start,intstren_end,ramptime)),("tmax",tmax_global) ])
        tevo_data,tevo_dict,_,saving_args = run_timeevo([states_starting[1],states_starting[2]],tevo_params,lattice_params_starting,tevo_hamilt_params; time_running_args...)

        # calculate final fidelity with target state
        final_gs_manifold = [tevo_data[1][1][:,end-1],tevo_data[1][2][:,end-1]]
        final_fidelity = groundstate_manifold_fidelity(final_gs_manifold,states_ending[1:speccount])
        println("Final fidelity for ramp time $(ramptime): $(final_fidelity)")
        
        #=scatter(ramptime,final_fidelity,c="b")
        end
        plot([0.0,10.0],[reference_fidelity,reference_fidelity],"--",c="r")
        xlabel("Ramp time")
        ylabel("Fidelity with target manifold")
        title("Fidelity vs ramp time $(lx)x$(ly) N=$(n) U=$(intstren_start)→$(intstren_end) ramp interaction strength")
        xscale("log")=#

    end

    #= time evolution with linear ramp looking at instantaneous energies
    if false || if_all
        speccount = 3
        dataloc = get_folder_location("cluster-data/exact-diag/time-evo")
        time_running_args = (nev=speccount,output_level=1,if_instant_gs=true,if_save_data=false,dataloc=dataloc,)

        ramptime = 10.0
        tmax_global = ramptime
        tevo_params = Dict([ ("interaction_strength",(linear_ramp,intstren_start,intstren_end,ramptime)),("tmax",tmax_global) ])
        tevo_data,tevo_dict,instdata,saving_args = run_timeevo([states_starting[1],states_starting[2],states_starting[3]],tevo_params,lattice_params_starting,hamilt_params_starting; time_running_args...)

        times = range(0.0,tmax_global,length=length(instdata[2]["1"]))
        cols = ["b","g","r"]
        for i in 1:speccount
            plot(times,instdata[2][string(i)],c=cols[i],"-p",label="E$(i)")
            plot(times,tevo_data[2][i][1:end-1],c="k",marker="x")
        end
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time for interaction strength ramp $(lx)x$(ly) N=$(n) ramptime $(ramptime) U $(intstren_start)→$(intstren_end)")
    end=#
end=#

#= look at 8x4 ulr flat linear ramp data
if false
    lx,ly,n = 8,4,4
    intstren_start = 10.0
    dataloc = get_folder_location("cluster-data/exact-diag/time-evo/ulr-flat")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("ramptype","linear"),("rampparam","interaction_strength")])
    all_files = find_data_file(pdict,"tevo",dataloc; file_type="jld2")
    
    for f in all_files
        rez = read_data(joinpath(dataloc,f); output_level=0)
        if haskey(rez,"final_fidelity")
            ff = rez["final_fidelity"]
            ramptime = get_params_dict_from_filename(f)["ramptime"]
            scatter(ramptime,ff,c="b")
            xlabel("Ramp time")
            ylabel("Fidelity with target manifold")
            title("Fidelity vs ramp time $(lx)x$(ly) N=$(n) U=$(intstren_start)→0.0 ramp ULR")
            xscale("log")
        end
    end

end=#


### Time evolution with the QuOCS AD-optimized interaction strength ramp from optimal-control/config_intstrenRamp_AD.py
# All parameters must match the ones the optimization ran with (see config_intstrenRamp_AD.py):
# 4x4 N=2 pbc, U ramped 10.0 -> 0.0 over ramptime 1.0, dt 0.005, intstren_max 12.0
#= Reference: dCRAB-optimized fidelity 0.9101 (run 20260714_114412), AD-optimized fidelity 0.9158 (run 20260717_114522, best of three AD runs / dCRAB fidelity is lower on the same 200-eval budget)
if false

    if_all::Bool = true

    # model parameters and starting/ending states
    if false || if_all
        lx,ly,n = 4,4,2
        intstren_start, intstren_end = 10.0, 0.0
        speccount_intad = 2

        pdict_intad = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren_start),("filling",0.5),("nev",3),("if_find_data",false),("if_save_data",false)])
        states_starting_intad,_,_,_,_,lattice_params_intad,hamilt_params_intad = run_normal_ed(pdict_intad; output_level=0)

        pdict_ending_intad = merge(pdict_intad,Dict("interaction_strength"=>intstren_end))
        states_ending_intad,_,_,_,_,_,_ = run_normal_ed(pdict_ending_intad; output_level=0)
    end

    # load the AD-optimized pulse and run the time evolution with instantaneous energies
    if false || if_all
        quocs_folder_intad = "../optimal-control/QuOCS_Results/20260717_114522_intstrenRamp_AD"
        controls_file_intad = filter(f -> endswith(f,"best_controls.npz"), readdir(quocs_folder_intad))[1]
        # only read the numeric arrays: NPZ.jl cannot parse the numpy unicode-string arrays
        # (pulse_names etc.) that QuOCS also stores in the file
        best_controls_intad = npzread(joinpath(quocs_folder_intad,controls_file_intad),["intstrenRamp","time_grid_for_intstrenRamp"])

        # AD mode hands get_FoM (and therefore dumps) the pulse as complex64 (quocslib
        # AD pitfall: Controls._get_controls_jax_obj uses complex jax tracers); the
        # interaction strength itself is real, so drop the (zero) imaginary part
        pulse_intad = Float64.(real.(best_controls_intad["intstrenRamp"]))

        # pulse is sampled on the RK4 half-step grid (spacing dt/2), so dt must match the
        # value used in config_intstrenRamp_AD.py: pulse length = ceil(2*ramptime/dt) + 1
        dt_intad = 0.005
        ramptime_intad = Float64(real(best_controls_intad["time_grid_for_intstrenRamp"][end]))

        # work on a copy: run_timeevo's timeham writes the ramp's current value back into
        # the dict, which would leave interaction_strength=0.0 for any later section
        hamilt_params_energy_intad = copy(hamilt_params_intad)

        time_running_args_intad = (nev=speccount_intad,output_level=1,if_instant_gs=true,if_save_data=false,dataloc="tevo-daily-things-data/")
        starting_states_intad = [Vector{ComplexF64}(states_starting_intad[i]) for i in 1:speccount_intad]
        tevo_params_intad = Dict([ ("interaction_strength",(pulse_ramp,ramptime_intad,pulse_intad)),("tmax",ramptime_intad),("dt",dt_intad) ])
        tevo_data_intad,tevo_dict_intad,instdata_intad,saving_args_intad = run_timeevo(starting_states_intad,tevo_params_intad,lattice_params_intad,hamilt_params_energy_intad; time_running_args_intad...)

        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        final_manifold_intad = [Vector{ComplexF64}(tevo_data_intad[1][i][:,end-1]) for i in 1:speccount_intad]
        fidelity_intad = real(groundstate_manifold_fidelity(final_manifold_intad,[Vector{ComplexF64}(s) for s in states_ending_intad[1:speccount_intad]]))
        println("Fidelity with target manifold using AD-optimized pulse: $(fidelity_intad)")
    end

    # plot the instantaneous vs transported (time-evolved) state energies along the AD ramp
    if true || if_all
        times_intad = range(0.0,ramptime_intad,length=length(instdata_intad[2]["1"]))

        figure()
        cols = ["b","g","r"]
        for i in 1:speccount_intad
            plot(times_intad,instdata_intad[2][string(i)],"-p",c=cols[i],label="E$(i) instantaneous")
            plot(times_intad,tevo_data_intad[2][i][1:end-1],c="k",marker="x",label=(i==1 ? "transported" : nothing))
        end
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time for AD interaction strength ramp $(lx)x$(ly) N=$(n) ramptime $(ramptime_intad) U $(intstren_start)→$(intstren_end), fidelity = $(round(fidelity_intad,digits=6))")
    end

end=#


### Time-dependent magnetic field gradient: pulse the gradient, then switch it off
# The gradient displaces the synthetic states, so the spacing entering the dipolar tail
# follows the driven response a(t) = a0 + int_0^t B(t') sin(t-t') dt' (accumulated by
# timeham, applied in the "magnetic_gradient" branch of long_range_scaling). B is ramped
# linearly 0.1 -> 10.0 and then switched straight off; because of the retarded sin kernel
# the displacement keeps ringing after the pulse ends rather than freezing at whatever it
# reached, so the interaction profile is still moving through the hold. At t = 0 the
# integral is zero, so the starting Hamiltonian is exactly the ordinary "dd" profile at
# spacing a0 -- the starting state is that dd groundstate, so the evolution begins in an
# eigenstate.
if true

    if_all::Bool = true

    # define starting state: dipole-dipole groundstate at the initial spacing a0
    if false || if_all
        lx,ly,n = 4,4,2
        intstren_maggrad = 300.0
        a0_maggrad = 0.1            # initial synthetic spacing, i.e. a(t=0)
        speccount_maggrad = 2       # low-lying states tracked through the pulse

        pdict_maggrad = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("scaling_type","dd"),("magnetic_spacing",a0_maggrad),("interaction_strength",intstren_maggrad),("filling",0.5),("nev",speccount_maggrad),("if_find_data",false),("if_save_data",false)])
        states_starting_maggrad,nrgs_starting_maggrad,_,_,_,lattice_params_maggrad,hamilt_params_maggrad = run_normal_ed(pdict_maggrad; output_level=0)

        # the t=0 gradient profile (zero accumulated integral) must reproduce the dd profile
        # the starting state was found with, otherwise the run does not start in an eigenstate
        us_at_zero_maggrad = long_range_scaling(ly-1,ly,intstren_maggrad; scaling="magnetic_gradient",magnetic_spacing=a0_maggrad,magnetic_gradient_integral=0.0)
        @assert us_at_zero_maggrad ≈ hamilt_params_maggrad["U"] "t=0 gradient profile $(us_at_zero_maggrad) does not match the dd profile $(hamilt_params_maggrad["U"])"
        println("Starting U profile: $(hamilt_params_maggrad["U"])")
    end

    # build the gradient pulse: linear ramp 0.1 -> 10.0 over ramptime, then off
    if false || if_all
        bstart_maggrad, bend_maggrad = 0.1, 10.0
        ramptime_maggrad = 1.0
        holdtime_maggrad = 0.5      # gradient sits at zero here, but a keeps ringing on the sin kernel
        tmax_maggrad = ramptime_maggrad + holdtime_maggrad
        # the step comes from the RK4 stability limit for the interaction scale rather than being
        # picked by hand: max(U) here is set by the tightest spacing the run visits, which is the
        # initial a0 since the gradient only pushes the states apart, so the t=0 profile is the
        # worst case and the same step is safe for the whole pulse. run_timeevo recomputes this
        # and errors if it is handed anything larger
        dt_maggrad = get_critical_dt(tmax_maggrad,lattice_params_maggrad,hamilt_params_maggrad)
        println("Critical time step: $(dt_maggrad) ($(Int(ceil(tmax_maggrad/dt_maggrad))) RK4 steps for tmax $(tmax_maggrad))")

        # pulse_ramp samples on the RK4 half-step grid (spacing dt/2) and timeham indexes it
        # directly, so the pulse must have exactly ceil(tmax/(dt/2)) + 1 entries
        bgrid_maggrad = collect(0:Int(ceil(tmax_maggrad/(dt_maggrad/2)))) .* (dt_maggrad/2)
        bpulse_maggrad = [t <= ramptime_maggrad ? bstart_maggrad + (bend_maggrad-bstart_maggrad)*(t/ramptime_maggrad) : 0.0 for t in bgrid_maggrad]

        # same retarded trapezoid rule as get_magnetic_gradient_integral, for the plots and
        # checks below: entry k is int_0^t_k B(t') sin(t_k-t') dt', which has to be redone at
        # every k because the kernel depends on the upper limit, so no running sum here. The
        # step size is uniform over the whole run (run_timeevo sets when_change_dt past the
        # last step), so bgrid_maggrad is exactly the sample-time vector that function builds.
        integral_maggrad = [k < 2 ? 0.0 :
            let integrand = bpulse_maggrad[1:k] .* sin.(bgrid_maggrad[k] .- bgrid_maggrad[1:k])
                0.5 * (dt_maggrad/2) * sum(integrand[1:k-1] .+ integrand[2:k])
            end for k in 1:length(bgrid_maggrad)]
        spacings_maggrad = a0_maggrad .+ integral_maggrad
        println("Spacing a: $(a0_maggrad) -> $(spacings_maggrad[end]) (peak $(maximum(spacings_maggrad)), still ringing at tmax)")
    end

    # time evolution along the pulse, tracking the instantaneous spectrum
    if false || if_all
        # work on a copy: timeham writes the ramped values, the accumulated integral and the
        # rebuilt U back into the dict, which would otherwise leave the starting-state
        # parameters overwritten for any later section
        hamilt_params_tevo_maggrad = copy(hamilt_params_maggrad)
        hamilt_params_tevo_maggrad["scaling_type"] = "magnetic_gradient"

        time_running_args_maggrad = (nev=speccount_maggrad,output_level=1,if_instant_gs=true,if_save_data=false,dataloc="tevo-daily-things-data/")
        starting_states_maggrad = [Vector{ComplexF64}(states_starting_maggrad[i]) for i in 1:speccount_maggrad]
        tevo_params_maggrad = Dict([ ("magnetic_gradient_time",(pulse_ramp,tmax_maggrad,bpulse_maggrad)),("tmax",tmax_maggrad),("dt",dt_maggrad) ])
        tevo_data_maggrad,tevo_dict_maggrad,instdata_maggrad,saving_args_maggrad = run_timeevo(starting_states_maggrad,tevo_params_maggrad,lattice_params_maggrad,hamilt_params_tevo_maggrad; time_running_args_maggrad...)

        # the U the run ended on must be the one the final accumulated spacing gives
        us_final_maggrad = long_range_scaling(ly-1,ly,intstren_maggrad; scaling="magnetic_gradient",magnetic_spacing=a0_maggrad,magnetic_gradient_integral=integral_maggrad[end])
        println("Final U profile: $(hamilt_params_tevo_maggrad["U"]) (expected $(us_final_maggrad))")
        @assert hamilt_params_tevo_maggrad["U"] ≈ us_final_maggrad "final U does not match the spacing the retarded integral ends on"

        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        final_states_maggrad = [Vector{ComplexF64}(tevo_data_maggrad[1][i][:,end-1]) for i in 1:speccount_maggrad]
        # the second starting state comes out of a near-degenerate pair, so which vector the ED
        # returns for it changes from run to run and a manifold-vs-manifold fidelity is not
        # reproducible here; the population of the final instantaneous manifold by the lowest
        # transported state is invariant under rotations inside that manifold, so track that
        final_instant_maggrad = [Vector{ComplexF64}(instdata_maggrad[1][string(i)][:,end]) for i in 1:2]
        manifold_population_maggrad = sum(abs2(dot(final_states_maggrad[1],instant_state)) for instant_state in final_instant_maggrad)
        println("Starting energies: $(nrgs_starting_maggrad)")
        println("Final instantaneous energies: $([instdata_maggrad[2][string(i)][end] for i in 1:speccount_maggrad])")
        println("Population of the final instantaneous manifold by the lowest transported state: $(manifold_population_maggrad)")
    end

    # the gradient and the spacing it accumulates
    if false || if_all
        figure()
        plot(bgrid_maggrad,bpulse_maggrad,c="b",label="gradient B")
        plot(bgrid_maggrad,spacings_maggrad,c="r",label="spacing a = a0 + ∫B(t')sin(t-t')dt'")
        legend()
        xlabel("Time")
        ylabel("Value")
        title("Gradient pulse $(bstart_maggrad)→$(bend_maggrad) over $(ramptime_maggrad), then off")
    end

    # the interaction profile at a few points along the pulse
    if false || if_all
        figure()
        # start / mid-ramp / end-of-ramp / end-of-hold, on the half-step grid
        snapshot_indices_maggrad = [1, Int(ceil(0.5*ramptime_maggrad/(dt_maggrad/2))), Int(ceil(ramptime_maggrad/(dt_maggrad/2))), length(bgrid_maggrad)]
        for idx in snapshot_indices_maggrad
            us_snapshot = long_range_scaling(ly-1,ly,intstren_maggrad; scaling="magnetic_gradient",magnetic_spacing=a0_maggrad,magnetic_gradient_integral=integral_maggrad[idx])
            plot(0:length(us_snapshot)-1,us_snapshot,"-p",label="t=$(round(bgrid_maggrad[idx],digits=3)), a=$(round(spacings_maggrad[idx],digits=3))")
        end
        legend()
        xlabel("y distance")
        ylabel("Interaction strength")
        yscale("log")
        title("ULR profile along the gradient pulse $(lx)x$(ly) N=$(n) U=$(intstren_maggrad)")
    end

    # instantaneous vs transported energies
    if false || if_all
        times_maggrad = range(0.0,tmax_maggrad,length=length(instdata_maggrad[2]["1"]))

        # once the gradient is off the spacing, and therefore H, is time independent, so both
        # the instantaneous spectrum and the transported energies have to be flat over the hold
        hold_start_maggrad = findfirst(t -> t >= ramptime_maggrad, times_maggrad)
        instant_drift_maggrad = maximum(abs.(instdata_maggrad[2]["1"][hold_start_maggrad:end] .- instdata_maggrad[2]["1"][end]))
        transported_drift_maggrad = maximum(abs.(tevo_data_maggrad[2][1][hold_start_maggrad:end-1] .- tevo_data_maggrad[2][1][end-1]))
        println("Energy drift over the hold: instantaneous $(instant_drift_maggrad), transported $(transported_drift_maggrad)")

        figure()
        cols = ["b","g","r"]
        for i in 1:speccount_maggrad
            plot(times_maggrad,instdata_maggrad[2][string(i)],"-p",c=cols[i],label="E$(i) instantaneous")
            plot(times_maggrad,tevo_data_maggrad[2][i][1:end-1],c="k",marker="x",label=(i==1 ? "transported" : nothing))
        end
        axvline(ramptime_maggrad,ls="--",c="k")
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time for gradient pulse $(lx)x$(ly) N=$(n) B $(bstart_maggrad)→$(bend_maggrad) over $(ramptime_maggrad), manifold population = $(round(manifold_population_maggrad,digits=6))")
    end

    #= animation of the density profile of the lowest transported state
    if false
        # column k of the transported data is the state after k RK4 steps, i.e. t = k*dt, and
        # the last column is the save point that never gets written, hence end-1; the starting
        # state is prepended so the animation actually opens at t = 0
        anim_states_maggrad = vcat([Vector{ComplexF64}(states_starting_maggrad[1])],[Vector{ComplexF64}(tevo_data_maggrad[1][1][:,k]) for k in 1:size(tevo_data_maggrad[1][1],2)-1])
        anim_times_maggrad = [(k-1)*dt_maggrad for k in 1:length(anim_states_maggrad)]

        # every frame_stride'th step only, otherwise this is a few hundred frames
        frame_stride_maggrad = 5
        frame_indices_maggrad = 1:frame_stride_maggrad:length(anim_states_maggrad)

        occs_frames_maggrad = [get_occupancy(anim_states_maggrad[k],lattice_params_maggrad; if_plot=false) for k in frame_indices_maggrad]
        # one colour scale for the whole animation, so frames can be compared by eye
        vmax_maggrad = maximum(maximum.(occs_frames_maggrad))

        # spacing at each frame time, to show it freezing when the gradient switches off
        # (spacings_maggrad lives on the half-step grid, so time t sits at index 2t/dt + 1)
        spacing_frames_maggrad = [spacings_maggrad[min(Int(round(2*anim_times_maggrad[k]/dt_maggrad))+1,length(spacings_maggrad))] for k in frame_indices_maggrad]

        animation_maggrad = PyPlot.PyCall.pyimport("matplotlib.animation")
        fig_maggrad = figure()
        img_maggrad = imshow(occs_frames_maggrad[1],origin="lower",vmin=0.0,vmax=vmax_maggrad)
        ax_maggrad = gca()
        colorbar()
        xlabel("Physical")
        ylabel("Synthetic")

        # matplotlib counts frames from zero
        function update_density_maggrad(frame)
            i = frame + 1
            img_maggrad.set_data(occs_frames_maggrad[i])
            ax_maggrad.set_title("Density of transported gs, t = $(round(anim_times_maggrad[frame_indices_maggrad[i]],digits=3)), a = $(round(spacing_frames_maggrad[i],digits=3))")
            return (img_maggrad,)
        end

        anim_maggrad = animation_maggrad.FuncAnimation(fig_maggrad,update_density_maggrad,frames=length(occs_frames_maggrad),interval=100)
        gifpath_maggrad = joinpath(get_folder_location("local-plots"),"tevo-density-maggrad-$(lx)x$(ly)-N-$(n)-B-$(bstart_maggrad)-to-$(bend_maggrad)-ramptime-$(ramptime_maggrad).gif")
        anim_maggrad.save(gifpath_maggrad,writer="pillow",fps=10)
        println("Saved density animation to $(gifpath_maggrad)")
    end=#

end






































"fin"
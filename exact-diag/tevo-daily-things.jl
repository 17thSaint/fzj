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
# Strategy: start with strongly interacting state (ULR) and ramp interaction strength to 0.0 to connect to the FCI ground state manifold.
if true
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
end#


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















































"fin"
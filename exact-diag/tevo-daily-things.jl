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
# follows the driven response a(t) = a0 + int_0^t B(t') sin(w(t-t'))/w dt' of a trapped mode
# of frequency w (accumulated by timeham, applied in the "magnetic_gradient" branch of
# long_range_scaling). B is ramped
# linearly 0.1 -> 10.0 and then switched straight off; because of the retarded sin kernel
# the displacement keeps ringing after the pulse ends rather than freezing at whatever it
# reached, so the interaction profile is still moving through the hold. At t = 0 the
#= integral is zero, so the starting Hamiltonian is exactly the ordinary "dd" profile at
# spacing a0 -- the starting state is that dd groundstate, so the evolution begins in an
# eigenstate.
if false

    if_all::Bool = true

    # define starting state: dipole-dipole groundstate at the initial spacing a0
    if true || if_all
        lx,ly,n = 4,4,2
        intstren_maggrad = 10.0
        a0_maggrad = 0.1            # initial synthetic spacing, i.e. a(t=0)
        speccount_maggrad = 2       # low-lying states tracked through the pulse
        omega = 10.0

        pdict_maggrad = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("scaling_type","dd"),("trap_frequency",omega),("magnetic_spacing",a0_maggrad),("interaction_strength",intstren_maggrad),("filling",0.5),("nev",speccount_maggrad),("if_find_data",false),("if_save_data",false)])
        states_starting_maggrad,nrgs_starting_maggrad,_,_,_,lattice_params_maggrad,hamilt_params_maggrad = run_normal_ed(pdict_maggrad; output_level=0)

        # the t=0 gradient profile (zero accumulated integral) must reproduce the dd profile
        # the starting state was found with, otherwise the run does not start in an eigenstate
        us_at_zero_maggrad = long_range_scaling(ly-1,ly,intstren_maggrad; scaling="magnetic_gradient",magnetic_spacing=a0_maggrad,magnetic_gradient_integral=0.0)
        @assert us_at_zero_maggrad ≈ hamilt_params_maggrad["U"] "t=0 gradient profile $(us_at_zero_maggrad) does not match the dd profile $(hamilt_params_maggrad["U"])"
        println("Starting U profile: $(hamilt_params_maggrad["U"])")
    end

    # build the gradient pulse: linear ramp 0.1 -> 10.0 over ramptime, then off
    if true || if_all
        bstart_maggrad, bend_maggrad = 0.1, 30.0
        # trap frequency of the displaced mode: sets the ringing period 2pi/w and scales the
        # response as 1/w^2. This is the only place to set it -- it feeds both the integral
        # below and tevo_params_maggrad, and it deliberately does not live in pdict_maggrad
        # because it has no effect on the t=0 Hamiltonian, so changing it never invalidates
        # the starting state and never needs the ED block above re-run
        w_maggrad = omega
        ramptime_maggrad = 0.1
        holdtime_maggrad = 2.0      # gradient sits at zero here, but a keeps ringing on the sin kernel
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
            let integrand = bpulse_maggrad[1:k] .* sin.(w_maggrad .* (bgrid_maggrad[k] .- bgrid_maggrad[1:k]))
                0.5 * (dt_maggrad/2) * sum(integrand[1:k-1] .+ integrand[2:k]) / w_maggrad
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
        tevo_params_maggrad = Dict([ ("magnetic_gradient_time",(pulse_ramp,tmax_maggrad,bpulse_maggrad)),("trap_frequency",w_maggrad),("tmax",tmax_maggrad),("dt",dt_maggrad) ])
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
    if true || if_all
        figure()
        plot(bgrid_maggrad,bpulse_maggrad,c="b",label="gradient B")
        plot(bgrid_maggrad,spacings_maggrad,c="r",label="spacing a = a0 + ∫B(t')sin(ω(t-t'))dt'/ω")
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

    # density in real space: synthetic index turned into a physical y position
    # get_occupancy resolves the density on (physical x, synthetic index m) only. The gradient
    # is what gives the synthetic direction a real extent: the displacement the retarded kernel
    # accumulates is exactly the ladder spacing a(t) = a0 + int_0^t B(t')sin(t-t')dt' that
    # enters the dipolar tail (U ~ 1/(a*x)^3 for synthetic separation x), so row m sits at
    # physical y_m(t) = (m-1)*a(t). The ladder therefore breathes with a(t): it stretches
    # through the ramp and keeps ringing over the hold even though the gradient is off, which
    # is motion of the cloud in real space that the synthetic-index picture cannot show.
    if false || if_all
        # column k of the transported data is the state after k RK4 steps, i.e. t = k*dt, and
        # the last column is the save point that never gets written, hence end-1; the starting
        # state is prepended so the series actually opens at t = 0
        states_physy_maggrad = vcat([Vector{ComplexF64}(states_starting_maggrad[1])],[Vector{ComplexF64}(tevo_data_maggrad[1][1][:,k]) for k in 1:size(tevo_data_maggrad[1][1],2)-1])
        times_physy_maggrad = [(k-1)*dt_maggrad for k in 1:length(states_physy_maggrad)]
        nframes_physy_maggrad = length(times_physy_maggrad)

        # occs[m,x] at every step: rows are the synthetic index, columns the physical x
        occs_physy_maggrad = [get_occupancy(s,lattice_params_maggrad; if_plot=false) for s in states_physy_maggrad]
        @assert all(occ -> sum(occ) ≈ n, occs_physy_maggrad) "occupancies do not sum to the particle number, so the transported states are not normalised"

        # cell edges in time (midpoints between steps) and the spacing evaluated there, so the
        # mesh below has one more edge than it has cells in each direction
        tedges_physy_maggrad = vcat(0.0,0.5 .* (times_physy_maggrad[1:end-1] .+ times_physy_maggrad[2:end]),times_physy_maggrad[end])
        # spacings_maggrad lives on the half-step grid, so time t sits at index 2t/dt + 1
        a_physy_maggrad = [spacings_maggrad[min(Int(round(2*t/dt_maggrad))+1,length(spacings_maggrad))] for t in times_physy_maggrad]
        aedges_physy_maggrad = [spacings_maggrad[min(Int(round(2*t/dt_maggrad))+1,length(spacings_maggrad))] for t in tedges_physy_maggrad]
        @assert a_physy_maggrad[1] ≈ a0_maggrad "the spacing at t=0 is not the initial spacing the starting state was found at"

        println("Physical y of the top synthetic row: $((ly-1)*a_physy_maggrad[1]) -> $((ly-1)*a_physy_maggrad[end]) (peak $((ly-1)*maximum(a_physy_maggrad)))")

        # space-time map: density summed over the physical x direction, against the physical y
        # the rows actually sit at. The cell edges move with a(t), so the rows fan out as the
        # gradient pushes them apart instead of sitting on fixed synthetic-index lines
        figure()
        dens_y_physy_maggrad = reduce(hcat,[vec(sum(occ,dims=2)) for occ in occs_physy_maggrad])
        tmesh_physy_maggrad = [tedges_physy_maggrad[k] for j in 1:ly+1, k in 1:nframes_physy_maggrad+1]
        ymesh_physy_maggrad = [(j-1.5)*aedges_physy_maggrad[k] for j in 1:ly+1, k in 1:nframes_physy_maggrad+1]
        pcolormesh(tmesh_physy_maggrad,ymesh_physy_maggrad,dens_y_physy_maggrad)
        colorbar(label="Density (summed over physical x)")
        axvline(ramptime_maggrad,ls="--",c="w")
        xlabel("Time")
        ylabel("Physical y")
        title("Real-space y density of the transported gs $(lx)x$(ly) N=$(n) B $(bstart_maggrad)→$(bend_maggrad) over $(ramptime_maggrad)")

        # the full real-space density at a few points along the pulse, on a common colour scale
        # and common axes so the panels can be compared by eye
        snapshot_times_physy_maggrad = [0.0, 0.5*ramptime_maggrad, ramptime_maggrad, tmax_maggrad]
        snapshot_frames_physy_maggrad = [argmin(abs.(times_physy_maggrad .- t)) for t in snapshot_times_physy_maggrad]
        vmax_physy_maggrad = maximum(maximum(occs_physy_maggrad[k]) for k in snapshot_frames_physy_maggrad)
        # physical x is the ordinary lattice direction, so its cells are one site wide
        xedges_physy_maggrad = collect(-0.5:1.0:lx-0.5)
        ylim_physy_maggrad = (ly-0.5) * maximum(a_physy_maggrad[k] for k in snapshot_frames_physy_maggrad)

        figure(figsize=(4*length(snapshot_frames_physy_maggrad),4))
        for (i,k) in enumerate(snapshot_frames_physy_maggrad)
            subplot(1,length(snapshot_frames_physy_maggrad),i)
            yedges_physy_maggrad = [(j-1.5)*a_physy_maggrad[k] for j in 1:ly+1]
            pcolormesh(xedges_physy_maggrad,yedges_physy_maggrad,occs_physy_maggrad[k],vmin=0.0,vmax=vmax_physy_maggrad)
            ylim(-0.5*a_physy_maggrad[1],ylim_physy_maggrad)
            xlabel("Physical x")
            i == 1 ? ylabel("Physical y") : nothing
            title("t=$(round(times_physy_maggrad[k],digits=3)), a=$(round(a_physy_maggrad[k],digits=3))")
        end
        # the loop leaves the last panel current, so this picks up its mesh; every panel is on
        # the same vmin/vmax, so the one bar reads for all of them
        colorbar(label="Density")
        suptitle("Real-space density of the transported gs $(lx)x$(ly) N=$(n) B $(bstart_maggrad)→$(bend_maggrad) over $(ramptime_maggrad)")
    end

end=#



### Posicast / zero-vibration two-step gradient pulse: settle the synthetic spacing exactly
# The retarded kernel a(t) = a0 + int_0^t B(t') sin(w(t-t'))/w dt' that long_range_scaling's
# "magnetic_gradient" branch is driven by is the Green's function of an undamped oscillator:
# x = a - a0 obeys xddot + w^2 x = B(t) from rest. Nothing removes energy from that mode, so
# whatever ringing a pulse leaves behind stays forever -- which is exactly what the previous
# (ramp-then-off) block suffers from, and why its interaction profile never stops moving.
#
# A step to B rings about the static displacement B/w^2 and so overshoots it by exactly a
# factor of two, reaching 2B/w^2 at wt = pi where the mode is momentarily at rest. A point at
# rest stays at rest if it sits at the equilibrium of the current dynamics, and we control
# where that equilibrium is, so switching at tsw = pi/w to the gradient whose equilibrium is
# the present position collapses the trajectory onto a fixed point:
#
#     B(t) = B1 = w^2*da/2   for 0 <= t < tsw,      B(t) = B2 = w^2*da = 2*B1   for t >= tsw,
#     tsw = pi/w,            da = atarg - a0.
#
# The second step is just the standing gradient that statically holds the target displacement;
# the first is precisely half of it, applied for precisely half a trap period. Three properties
# matter here specifically because U ~ 1/a^3:
#   - xdot = (B1/w)*sin(wt) >= 0 on [0,tsw], so a climbs monotonically from a0 to atarg with no
#     overshoot and no zero crossing, and U never diverges;
#   - min_t a(t) = a0 at t = 0 only, so max_t U(t) = U(a0) and the RK4 step computed from the
#     t = 0 Hamiltonian bounds the whole run -- which is the assumption get_critical_dt makes;
#   - after tsw the Hamiltonian is genuinely static, not merely correct on average.
# This is Smith's 1957 posicast controller / the ZV input shaper, and in the cold-atom setting
# the same logic as shortcuts to adiabaticity for transport: rather than moving slowly enough
# never to excite the mode, excite it deliberately and then cancel the excitation exactly.
#
#= Caveat, with the measured outcome: with atarg = 2.0 and a0 = 0.1 the nearest-neighbour
# coupling falls 10/(0.1)^3 = 10000 -> 10/(2.0)^3 = 1.25 in half a trap period. That is a
# violent quench, but it costs less than a first look suggests -- the lowest transported state
# still ends with 0.857 of its weight in the final instantaneous manifold. The pulse is
# designed for a(t), not for fidelity; lower atarg if adiabaticity is what is wanted.
if false

    if_all::Bool = true

    # define starting state: dipole-dipole groundstate at the initial spacing a0
    if true || if_all
        lx,ly,n = 4,4,2
        intstren_posicast = 10.0
        a0_posicast = 0.1            # initial synthetic spacing, i.e. a(t=0)
        speccount_posicast = 2       # low-lying states tracked through the pulse
        omega_posicast = 10.0        # trap frequency of the displaced mode

        pdict_posicast = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("scaling_type","dd"),("trap_frequency",omega_posicast),("magnetic_spacing",a0_posicast),("interaction_strength",intstren_posicast),("filling",0.5),("nev",speccount_posicast),("if_find_data",false),("if_save_data",false)])
        states_starting_posicast,nrgs_starting_posicast,_,_,_,lattice_params_posicast,hamilt_params_posicast = run_normal_ed(pdict_posicast; output_level=0)

        # the t=0 gradient profile (zero accumulated integral) must reproduce the dd profile
        # the starting state was found with, otherwise the run does not start in an eigenstate
        us_at_zero_posicast = long_range_scaling(ly-1,ly,intstren_posicast; scaling="magnetic_gradient",magnetic_spacing=a0_posicast,magnetic_gradient_integral=0.0)
        @assert us_at_zero_posicast ≈ hamilt_params_posicast["U"] "t=0 gradient profile $(us_at_zero_posicast) does not match the dd profile $(hamilt_params_posicast["U"])"
        println("Starting U profile: $(hamilt_params_posicast["U"])")
    end

    # build the two-step pulse and the step it lives on
    if true || if_all
        atarg_posicast = 2.0                                    # where the spacing is to park
        da_posicast = atarg_posicast - a0_posicast
        b1_posicast = omega_posicast^2 * da_posicast / 2        # first step: half the holding gradient
        b2_posicast = omega_posicast^2 * da_posicast            # second step: the standing gradient
        tsw_posicast = pi / omega_posicast                      # switch at half a trap period

        # the binding amplitude constraint is the second step, not the first: holding a displaced
        # spacing needs a standing gradient, so atarg <= a0 + Bmax/w^2 for any non-ringing
        # solution, shaped or not. The two-step pulse buys a perfect trajectory, not a larger
        # reachable target. Bmax = 200 is the ceiling the dCRAB runs in config_maggradPulse.py use
        bmax_posicast = 200.0
        @assert b2_posicast <= bmax_posicast "holding gradient B2 = $(b2_posicast) exceeds Bmax = $(bmax_posicast): the largest spacing that can settle at all is a0 + Bmax/w^2 = $(a0_posicast + bmax_posicast/omega_posicast^2)"

        # settling is complete at tsw, so the remaining half-periods are there purely to
        # demonstrate that nothing rings afterwards and to give the state time to evolve under
        # the now-static Hamiltonian
        nhalfperiods_posicast = 4
        tmax_posicast = nhalfperiods_posicast * tsw_posicast

        # the step comes from the RK4 stability limit for the interaction scale rather than being
        # picked by hand. max(U) over the run is set by the tightest spacing visited, and the
        # pulse rises monotonically from a0, so the t=0 profile is the worst case and the same
        # step is safe throughout -- this is the property the shape was chosen for
        dtcrit_posicast = get_critical_dt(tmax_posicast,lattice_params_posicast,hamilt_params_posicast)

        # ... then snapped down so that tsw is an exact multiple of the half-step spacing dt/2.
        # The catch only works if the switch is applied at the instant the mode is at rest, and
        # timeham can only see the pulse on its sample grid, so the switch must land on a sample
        nsw_posicast = Int(ceil(tsw_posicast / (dtcrit_posicast/2)))
        dt_posicast = 2 * tsw_posicast / nsw_posicast
        @assert dt_posicast <= dtcrit_posicast "snapped dt $(dt_posicast) is above the RK4 stability limit $(dtcrit_posicast)"

        # pulse_ramp samples on the RK4 half-step grid (spacing dt/2) and timeham indexes it
        # directly, so the pulse must have exactly ceil(tmax/(dt/2)) + 1 entries -- the same
        # expression pulse_ramp recomputes, so the two agree whichever way the ceil rounds
        nhalf_posicast = Int(ceil(tmax_posicast/(dt_posicast/2)))
        bgrid_posicast = collect(0:nhalf_posicast) .* (dt_posicast/2)
        isw_posicast = nsw_posicast + 1     # 1-based index of the sample sitting on tsw
        @assert bgrid_posicast[isw_posicast] ≈ tsw_posicast "the switch sample sits at $(bgrid_posicast[isw_posicast]) rather than tsw = $(tsw_posicast)"

        bpulse_posicast = fill(b2_posicast,nhalf_posicast+1)
        bpulse_posicast[1:isw_posicast-1] .= b1_posicast
        # the one sample sitting on the discontinuity carries the mean of the two steps. This is
        # not a fudge: for the trapezoid rule the two intervals either side of tsw then integrate
        # to exactly B1*int_{tsw-h}^{tsw} + B2*int_{tsw}^{tsw+h}, i.e. the true step function,
        # instead of ramping B1 -> B2 across a half-step. Leaving the full B2 here is worse than
        # not snapping dt at all, because snapping puts a sample exactly where the error is
        # largest: measured residual over the parked stretch is 3e-4 with the plain rule against
        # 3e-16 with this one (and 6e-5 with neither fix, where the switch misses the grid by luck)
        bpulse_posicast[isw_posicast] = 0.5*(b1_posicast + b2_posicast)

        println("Two-step pulse: B1 = $(b1_posicast), B2 = $(b2_posicast), tsw = $(round(tsw_posicast,digits=5))")
        println("Critical time step: $(dtcrit_posicast), snapped to $(dt_posicast) ($(Int(ceil(tmax_posicast/dt_posicast))) RK4 steps for tmax $(round(tmax_posicast,digits=5)))")
    end

    # the spacing the pulse accumulates, and the checks that it is the intended one
    if true || if_all
        # same retarded trapezoid rule as get_magnetic_gradient_integral, evaluated here for the
        # plots and checks below. Written O(n) rather than O(n^2) by splitting the kernel,
        # sin(w(t-t')) = sin(wt)cos(wt') - cos(wt)sin(wt'), which turns the retarded integral into
        # two cumulative trapezoids over integrands that do not depend on the upper limit. Same
        # samples and same weights as the literal rule, so it is the same quadrature, not an
        # approximation to it -- the spot checks below confirm that against the function itself
        cumtrap_posicast = f -> pushfirst!(cumsum(0.5*(dt_posicast/2) .* (f[1:end-1] .+ f[2:end])),0.0)
        cosgrid_posicast, singrid_posicast = cos.(omega_posicast .* bgrid_posicast), sin.(omega_posicast .* bgrid_posicast)
        ic_posicast = cumtrap_posicast(bpulse_posicast .* cosgrid_posicast)
        is_posicast = cumtrap_posicast(bpulse_posicast .* singrid_posicast)
        integral_posicast = (singrid_posicast .* ic_posicast .- cosgrid_posicast .* is_posicast) ./ omega_posicast
        spacings_posicast = a0_posicast .+ integral_posicast

        # check it against get_magnetic_gradient_integral itself, i.e. against the quadrature the
        # evolution actually runs, rather than against a second copy of my own arithmetic. A
        # uniform step throughout is what run_timeevo sets up (when_change_dt lands past the last
        # step), so when_dt_ends is put beyond the end here to match
        checkparams_posicast = Dict{String,Any}([("magnetic_gradient_time",bpulse_posicast),("when_dt_ends",[nhalf_posicast+1,nhalf_posicast+1]),("dt",[dt_posicast,dt_posicast]),("trap_frequency",omega_posicast)])
        for k in [1, 2, isw_posicast÷2, isw_posicast, isw_posicast+1, nhalf_posicast+1]
            @assert isapprox(get_magnetic_gradient_integral(k,checkparams_posicast),integral_posicast[k]; atol=1e-12) "the O(n) kernel split disagrees with get_magnetic_gradient_integral at sample $k"
        end

        # against the closed form: x = (da/2)(1 - cos wt) before the switch, x = da after it
        analytic_posicast = [t < tsw_posicast ? a0_posicast + 0.5*da_posicast*(1-cos(omega_posicast*t)) : atarg_posicast for t in bgrid_posicast]
        quaderr_posicast = maximum(abs.(spacings_posicast .- analytic_posicast))
        # the residual ringing: peak-to-peak spread of a over everything past the switch. Zero in
        # exact arithmetic, and what any imperfect pulse would show as a surviving oscillation
        ringing_posicast = maximum(spacings_posicast[isw_posicast:end]) - minimum(spacings_posicast[isw_posicast:end])
        println("Spacing a: $(a0_posicast) -> $(spacings_posicast[end]) (target $(atarg_posicast), min over the run $(minimum(spacings_posicast)))")
        println("Deviation from the closed form: $(quaderr_posicast); residual ringing after the switch (peak-to-peak): $(ringing_posicast)")
        @assert quaderr_posicast < 1e-5 "the sampled pulse does not reproduce the closed-form response to 1e-5"
        @assert ringing_posicast < 1e-9 "the spacing still rings by $(ringing_posicast) after the switch, so the catch is landing at the wrong phase"
        # the monotone rise is what makes the t=0 RK4 step valid for the whole run
        @assert minimum(spacings_posicast) ≈ a0_posicast "the spacing dips below a0, so max(U) is not U(a0) and the critical dt no longer bounds the run"
    end

    # time evolution along the pulse, tracking the instantaneous spectrum
    if false || if_all
        # work on a copy: timeham writes the ramped values, the accumulated integral and the
        # rebuilt U back into the dict, which would otherwise leave the starting-state
        # parameters overwritten for any later section
        hamilt_params_tevo_posicast = copy(hamilt_params_posicast)
        hamilt_params_tevo_posicast["scaling_type"] = "magnetic_gradient"

        time_running_args_posicast = (nev=speccount_posicast,output_level=1,if_instant_gs=true,if_save_data=false,dataloc="tevo-daily-things-data/")
        starting_states_posicast = [Vector{ComplexF64}(states_starting_posicast[i]) for i in 1:speccount_posicast]
        tevo_params_posicast = Dict([ ("magnetic_gradient_time",(pulse_ramp,tmax_posicast,bpulse_posicast)),("trap_frequency",omega_posicast),("tmax",tmax_posicast),("dt",dt_posicast) ])
        tevo_data_posicast,tevo_dict_posicast,instdata_posicast,saving_args_posicast = run_timeevo(starting_states_posicast,tevo_params_posicast,lattice_params_posicast,hamilt_params_tevo_posicast; time_running_args_posicast...)

        # the U the run ended on must be the profile the target spacing gives. Comparing against
        # atarg rather than against integral_posicast[end] is the point of the whole construction:
        # the final Hamiltonian is the static one for the target, not whatever the ringing left
        us_final_posicast = long_range_scaling(ly-1,ly,intstren_posicast; scaling="magnetic_gradient",magnetic_spacing=a0_posicast,magnetic_gradient_integral=da_posicast)
        println("Final U profile: $(hamilt_params_tevo_posicast["U"]) (expected $(us_final_posicast))")
        @assert hamilt_params_tevo_posicast["U"] ≈ us_final_posicast "final U is not the static profile at the target spacing"

        # end-1 skips the final save column, which is allocated but never written
        final_states_posicast = [Vector{ComplexF64}(tevo_data_posicast[1][i][:,end-1]) for i in 1:speccount_posicast]
        # the second starting state comes out of a near-degenerate pair, so which vector the ED
        # returns for it changes from run to run and a manifold-vs-manifold fidelity is not
        # reproducible here; the population of the final instantaneous manifold by the lowest
        # transported state is invariant under rotations inside that manifold, so track that
        final_instant_posicast = [Vector{ComplexF64}(instdata_posicast[1][string(i)][:,end]) for i in 1:2]
        manifold_population_posicast = sum(abs2(dot(final_states_posicast[1],instant_state)) for instant_state in final_instant_posicast)
        println("Starting energies: $(nrgs_starting_posicast)")
        println("Final instantaneous energies: $([instdata_posicast[2][string(i)][end] for i in 1:speccount_posicast])")
        println("Population of the final instantaneous manifold by the lowest transported state: $(manifold_population_posicast)")
    end

    # the two-step gradient and the spacing it accumulates
    if false || if_all
        figure()
        plot(bgrid_posicast,bpulse_posicast,c="b",label="gradient B")
        plot(bgrid_posicast,spacings_posicast,c="r",label="spacing a = a0 + ∫B(t')sin(ω(t-t'))dt'/ω")
        axhline(atarg_posicast,ls=":",c="r",label="target a = $(atarg_posicast)")
        axvline(tsw_posicast,ls="--",c="k",label="switch tsw = π/ω")
        legend()
        xlabel("Time")
        ylabel("Value")
        title("Two-step pulse B1=$(b1_posicast)→B2=$(b2_posicast) at tsw=$(round(tsw_posicast,digits=4)), ω=$(omega_posicast), ringing $(round(ringing_posicast,sigdigits=2))")
    end

    # the interaction profile at a few points along the pulse
    if false || if_all
        figure()
        # start / mid-rise / switch / end, on the half-step grid
        snapshot_indices_posicast = [1, Int(round(0.5*tsw_posicast/(dt_posicast/2)))+1, isw_posicast, nhalf_posicast+1]
        for idx in snapshot_indices_posicast
            us_snapshot = long_range_scaling(ly-1,ly,intstren_posicast; scaling="magnetic_gradient",magnetic_spacing=a0_posicast,magnetic_gradient_integral=integral_posicast[idx])
            plot(0:length(us_snapshot)-1,us_snapshot,"-p",label="t=$(round(bgrid_posicast[idx],digits=3)), a=$(round(spacings_posicast[idx],digits=3))")
        end
        legend()
        xlabel("y distance")
        ylabel("Interaction strength")
        yscale("log")
        title("ULR profile along the two-step pulse $(lx)x$(ly) N=$(n) U=$(intstren_posicast)")
    end

    # instantaneous vs transported energies
    if false || if_all
        # instdata holds one entry per RK4 step, written at t = k*dt for k = 1..nsteps, so the
        # series starts one step in rather than at t = 0
        times_posicast = [k*dt_posicast for k in 1:length(instdata_posicast[2]["1"])]

        # past the switch the spacing, and therefore H, is time independent, so both the
        # instantaneous spectrum and the transported energies have to be flat from tsw onward.
        # This is the check that distinguishes "settled" from "correct on average"
        settled_start_posicast = findfirst(t -> t >= tsw_posicast, times_posicast)
        instant_drift_posicast = maximum(abs.(instdata_posicast[2]["1"][settled_start_posicast:end] .- instdata_posicast[2]["1"][end]))
        transported_drift_posicast = maximum(abs.(tevo_data_posicast[2][1][settled_start_posicast:end-1] .- tevo_data_posicast[2][1][end-1]))
        println("Energy drift after the switch: instantaneous $(instant_drift_posicast), transported $(transported_drift_posicast)")

        figure()
        cols = ["b","g","r"]
        for i in 1:speccount_posicast
            plot(times_posicast,instdata_posicast[2][string(i)],"-p",c=cols[i],label="E$(i) instantaneous")
            plot(times_posicast,tevo_data_posicast[2][i][1:end-1],c="k",marker="x",label=(i==1 ? "transported" : nothing))
        end
        axvline(tsw_posicast,ls="--",c="k")
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time, two-step pulse $(lx)x$(ly) N=$(n) a $(a0_posicast)→$(atarg_posicast), manifold population = $(round(manifold_population_posicast,digits=6))")
    end

    # density in real space: synthetic index turned into a physical y position
    # get_occupancy resolves the density on (physical x, synthetic index m) only. The gradient
    # is what gives the synthetic direction a real extent: the displacement the retarded kernel
    # accumulates is exactly the ladder spacing a(t) that enters the dipolar tail, so row m sits
    # at physical y_m(t) = (m-1)*a(t). Unlike the ramp-then-off pulse, this ladder stretches
    # monotonically to its target and then stops dead -- the rows fan out over the first half
    # period and sit still for the rest of the run
    if false || if_all
        # column k of the transported data is the state after k RK4 steps, i.e. t = k*dt, and
        # the last column is the save point that never gets written, hence end-1; the starting
        # state is prepended so the series actually opens at t = 0
        states_physy_posicast = vcat([Vector{ComplexF64}(states_starting_posicast[1])],[Vector{ComplexF64}(tevo_data_posicast[1][1][:,k]) for k in 1:size(tevo_data_posicast[1][1],2)-1])
        times_physy_posicast = [(k-1)*dt_posicast for k in 1:length(states_physy_posicast)]
        nframes_physy_posicast = length(times_physy_posicast)

        # occs[m,x] at every step: rows are the synthetic index, columns the physical x
        occs_physy_posicast = [get_occupancy(s,lattice_params_posicast; if_plot=false) for s in states_physy_posicast]
        @assert all(occ -> sum(occ) ≈ n, occs_physy_posicast) "occupancies do not sum to the particle number, so the transported states are not normalised"

        # cell edges in time (midpoints between steps) and the spacing evaluated there, so the
        # mesh below has one more edge than it has cells in each direction
        tedges_physy_posicast = vcat(0.0,0.5 .* (times_physy_posicast[1:end-1] .+ times_physy_posicast[2:end]),times_physy_posicast[end])
        # spacings_posicast lives on the half-step grid, so time t sits at index 2t/dt + 1
        a_physy_posicast = [spacings_posicast[min(Int(round(2*t/dt_posicast))+1,length(spacings_posicast))] for t in times_physy_posicast]
        aedges_physy_posicast = [spacings_posicast[min(Int(round(2*t/dt_posicast))+1,length(spacings_posicast))] for t in tedges_physy_posicast]
        @assert a_physy_posicast[1] ≈ a0_posicast "the spacing at t=0 is not the initial spacing the starting state was found at"

        println("Physical y of the top synthetic row: $((ly-1)*a_physy_posicast[1]) -> $((ly-1)*a_physy_posicast[end])")

        # space-time map: density summed over the physical x direction, against the physical y
        # the rows actually sit at
        figure()
        dens_y_physy_posicast = reduce(hcat,[vec(sum(occ,dims=2)) for occ in occs_physy_posicast])
        tmesh_physy_posicast = [tedges_physy_posicast[k] for j in 1:ly+1, k in 1:nframes_physy_posicast+1]
        ymesh_physy_posicast = [(j-1.5)*aedges_physy_posicast[k] for j in 1:ly+1, k in 1:nframes_physy_posicast+1]
        pcolormesh(tmesh_physy_posicast,ymesh_physy_posicast,dens_y_physy_posicast)
        colorbar(label="Density (summed over physical x)")
        axvline(tsw_posicast,ls="--",c="w")
        xlabel("Time")
        ylabel("Physical y")
        title("Real-space y density of the transported gs $(lx)x$(ly) N=$(n), two-step a $(a0_posicast)→$(atarg_posicast)")

        # the full real-space density at a few points along the pulse, on a common colour scale
        # and common axes so the panels can be compared by eye
        snapshot_times_physy_posicast = [0.0, 0.5*tsw_posicast, tsw_posicast, tmax_posicast]
        snapshot_frames_physy_posicast = [argmin(abs.(times_physy_posicast .- t)) for t in snapshot_times_physy_posicast]
        vmax_physy_posicast = maximum(maximum(occs_physy_posicast[k]) for k in snapshot_frames_physy_posicast)
        # physical x is the ordinary lattice direction, so its cells are one site wide
        xedges_physy_posicast = collect(-0.5:1.0:lx-0.5)
        ylim_physy_posicast = (ly-0.5) * maximum(a_physy_posicast[k] for k in snapshot_frames_physy_posicast)

        figure(figsize=(4*length(snapshot_frames_physy_posicast),4))
        for (i,k) in enumerate(snapshot_frames_physy_posicast)
            subplot(1,length(snapshot_frames_physy_posicast),i)
            yedges_physy_posicast = [(j-1.5)*a_physy_posicast[k] for j in 1:ly+1]
            pcolormesh(xedges_physy_posicast,yedges_physy_posicast,occs_physy_posicast[k],vmin=0.0,vmax=vmax_physy_posicast)
            ylim(-0.5*a_physy_posicast[1],ylim_physy_posicast)
            xlabel("Physical x")
            i == 1 ? ylabel("Physical y") : nothing
            title("t=$(round(times_physy_posicast[k],digits=3)), a=$(round(a_physy_posicast[k],digits=3))")
        end
        # the loop leaves the last panel current, so this picks up its mesh; every panel is on
        # the same vmin/vmax, so the one bar reads for all of them
        colorbar(label="Density")
        suptitle("Real-space density of the transported gs $(lx)x$(ly) N=$(n), two-step a $(a0_posicast)→$(atarg_posicast)")
    end

end=#




### ZVD three-step gradient pulse: settle the spacing and tolerate a mis-measured trap frequency
# Same retarded kernel as the two-step block -- x = a - a0 obeys xddot + w^2 x = B(t) from rest,
# undamped, so any ringing left at the end of the pulse is permanent. The two-step (posicast /
# ZV) catch is exact but fragile: it relies on pi/w being the true half period, and at an actual
# frequency w(1+eps) the switch lands off the rest point, leaving a residual linear in eps.
#
# Input shaping puts this in the right frame: a shaped command is the step convolved with an
# impulse train A_i at times t_i, and the residual oscillation it leaves is proportional to
# sum_i A_i exp(i w t_i). The ZV shaper (1/2, 1/2 at 0, T/2) zeroes that sum -- that is exactly
# the two-step pulse. The ZVD shaper (1/4, 1/2, 1/4 at 0, T/2, T with T = 2pi/w) zeroes the sum
# *and* its derivative with respect to w, so the residual is O(eps^2) instead of O(eps). It buys
# that with one extra half period of settling time.
#
# Convolving a step of final amplitude Bf = w^2*da with that train gives the staircase
#
#     B(t) = Bf * [1/4, 3/4, 1]   on   [0, T/2), [T/2, T), [T, inf),
#
# the cumulative sums of the impulse amplitudes. The response is the matching superposition of
# shifted step responses, x(t) = da * sum_{t_i <= t} A_i*(1 - cos w(t-t_i)), which parks at da/2
# at T/2 and at da from T onward. Both segments rise monotonically, so as with the two-step pulse
# min_t a(t) = a0 at t = 0 only: max_t U(t) = U(a0) and the RK4 step from the t = 0 Hamiltonian
# bounds the whole run. Everything below is written over the impulse train rather than over
# spelled-out branches, so the same code builds either shaper from its (A_i, t_i) list.
#
# Note the O(eps^2) robustness is the reason for the extra step but is asserted here, not
#= demonstrated -- nothing in this block drives the pulse at a detuned frequency.
if false

    if_all::Bool = true

    # starting state: dipole-dipole groundstate at the initial spacing a0, which the t=0 gradient
    # profile (zero accumulated integral) has to reproduce or the run does not start in an eigenstate
    if true || if_all
        lx,ly,n = 4,4,2
        intstren_zvd = 10.0
        a0_zvd, atarg_zvd = 0.1, 2.0
        omega_zvd = 10.0
        speccount_zvd = 2

        pdict_zvd = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("scaling_type","dd"),("trap_frequency",omega_zvd),("magnetic_spacing",a0_zvd),("interaction_strength",intstren_zvd),("filling",0.5),("nev",speccount_zvd),("if_find_data",false),("if_save_data",false)])
        states_starting_zvd,nrgs_starting_zvd,_,_,_,lattice_params_zvd,hamilt_params_zvd = run_normal_ed(pdict_zvd; output_level=0)
        @assert long_range_scaling(ly-1,ly,intstren_zvd; scaling="magnetic_gradient",magnetic_spacing=a0_zvd,magnetic_gradient_integral=0.0) ≈ hamilt_params_zvd["U"] "t=0 gradient profile does not match the dd profile the starting state was found with"
        println("Starting U profile: $(hamilt_params_zvd["U"])")
    end

    # the impulse train, the staircase it convolves to, and the spacing that comes back
    if true || if_all
        da_zvd = atarg_zvd - a0_zvd
        bf_zvd = omega_zvd^2 * da_zvd               # final standing gradient, = peak |B|
        amps_zvd = [0.25, 0.5, 0.25]                # ZVD impulse amplitudes, as fractions of Bf
        timps_zvd = [0.0, 1.0, 2.0] .* (pi/omega_zvd)

        # the two shaper conditions, checked rather than trusted: the first kills the residual
        # amplitude, the second its first derivative in w (which is what ZV lacks)
        @assert abs(sum(a*cis(omega_zvd*t) for (a,t) in zip(amps_zvd,timps_zvd))) < 1e-12 "impulse train does not satisfy the zero-vibration condition"
        @assert abs(sum(a*t*cis(omega_zvd*t) for (a,t) in zip(amps_zvd,timps_zvd))) < 1e-12 "impulse train does not satisfy the zero-derivative condition"
        # holding a displaced spacing needs a standing gradient, so atarg <= a0 + Bmax/w^2 caps
        # any non-ringing solution, shaped or not; ZVD does not raise that ceiling
        @assert bf_zvd <= 200.0 "standing gradient Bf = $(bf_zvd) exceeds Bmax = 200.0; largest spacing that can settle is $(a0_zvd + 200.0/omega_zvd^2)"

        # settling completes at the last impulse, T = 2pi/w; the remaining half period is there to
        # show nothing rings afterwards. dt comes from the RK4 stability limit for the interaction
        # scale, then is snapped down so that pi/w sits an exact whole number of half-steps away --
        # every impulse time is a multiple of pi/w, so one snap puts all of them on the grid, and
        # the switches have to land on samples or the catches are applied at the wrong phase
        tmax_zvd = 3 * (pi/omega_zvd)
        dtcrit_zvd = get_critical_dt(tmax_zvd,lattice_params_zvd,hamilt_params_zvd)
        nsw_zvd = Int(ceil((pi/omega_zvd) / (dtcrit_zvd/2)))
        dt_zvd = 2 * (pi/omega_zvd) / nsw_zvd
        @assert dt_zvd <= dtcrit_zvd "snapped dt $(dt_zvd) is above the RK4 stability limit $(dtcrit_zvd)"

        # pulse_ramp samples on the half-step grid (spacing dt/2) and timeham indexes it directly,
        # so the pulse needs exactly ceil(tmax/(dt/2)) + 1 entries -- the same expression pulse_ramp
        # recomputes, so the two agree whichever way the ceil rounds
        nhalf_zvd = Int(ceil(tmax_zvd/(dt_zvd/2)))
        bgrid_zvd = collect(0:nhalf_zvd) .* (dt_zvd/2)
        iimps_zvd = Int.(round.(timps_zvd ./ (dt_zvd/2))) .+ 1
        @assert bgrid_zvd[iimps_zvd] ≈ timps_zvd "the impulse samples sit at $(bgrid_zvd[iimps_zvd]) rather than $(timps_zvd)"

        # accumulate the staircase impulse by impulse. Each interior jump sample carries half the
        # amplitude it is jumping by, which is not a fudge: for the trapezoid rule the two intervals
        # either side then integrate to exactly the true step function rather than ramping across a
        # half-step, and without it the residual is O(w*da*h) -- worse than not snapping dt at all,
        # since snapping puts a sample exactly where that error is largest. t = 0 is the exception
        # and takes the full amplitude: there is no "before" side inside the integration domain, so
        # halving it would just lose half of the first interval
        bpulse_zvd = zeros(Float64,nhalf_zvd+1)
        for (amp,idx) in zip(amps_zvd,iimps_zvd)
            bpulse_zvd[idx] += bf_zvd * amp * (idx == 1 ? 1.0 : 0.5)
            bpulse_zvd[idx+1:end] .+= bf_zvd * amp
        end

        # same retarded trapezoid as get_magnetic_gradient_integral, written O(n) instead of O(n^2)
        # by splitting sin(w(t-t')) = sin(wt)cos(wt') - cos(wt)sin(wt'), which turns the retarded
        # integral into two cumulative trapezoids over integrands independent of the upper limit.
        # Same samples and weights, so it is the same quadrature rather than an approximation to it
        cumtrap_zvd = f -> pushfirst!(cumsum(0.5*(dt_zvd/2) .* (f[1:end-1] .+ f[2:end])),0.0)
        cosg_zvd, sing_zvd = cos.(omega_zvd .* bgrid_zvd), sin.(omega_zvd .* bgrid_zvd)
        integral_zvd = (sing_zvd .* cumtrap_zvd(bpulse_zvd .* cosg_zvd) .- cosg_zvd .* cumtrap_zvd(bpulse_zvd .* sing_zvd)) ./ omega_zvd
        spacings_zvd = a0_zvd .+ integral_zvd

        # check it against get_magnetic_gradient_integral itself, i.e. against the quadrature the
        # evolution actually runs, not a second copy of the same arithmetic. A uniform step
        # throughout is what run_timeevo sets up, so when_dt_ends goes past the end here to match
        checkparams_zvd = Dict{String,Any}([("magnetic_gradient_time",bpulse_zvd),("when_dt_ends",[nhalf_zvd+1,nhalf_zvd+1]),("dt",[dt_zvd,dt_zvd]),("trap_frequency",omega_zvd)])
        for k in vcat(1,2,iimps_zvd,iimps_zvd.+1,nhalf_zvd+1)
            @assert isapprox(get_magnetic_gradient_integral(k,checkparams_zvd),integral_zvd[k]; atol=1e-12) "the O(n) kernel split disagrees with get_magnetic_gradient_integral at sample $k"
        end

        # closed form: the superposition of step responses from the impulses that have fired
        analytic_zvd = [a0_zvd + da_zvd*sum(amp*(1-cos(omega_zvd*(t-ti))) for (amp,ti) in zip(amps_zvd,timps_zvd) if t >= ti; init=0.0) for t in bgrid_zvd]
        quaderr_zvd = maximum(abs.(spacings_zvd .- analytic_zvd))
        # residual ringing: peak-to-peak spread of a past the last impulse, zero in exact arithmetic
        settled_zvd = iimps_zvd[end]
        ringing_zvd = maximum(spacings_zvd[settled_zvd:end]) - minimum(spacings_zvd[settled_zvd:end])
        println("ZVD staircase: Bf*[1/4,3/4,1] = $(bf_zvd .* cumsum(amps_zvd)) at t = $(round.(timps_zvd,digits=5))")
        println("Critical time step: $(dtcrit_zvd), snapped to $(dt_zvd) ($(Int(ceil(tmax_zvd/dt_zvd))) RK4 steps for tmax $(round(tmax_zvd,digits=5)))")
        println("Spacing a: $(a0_zvd) -> $(spacings_zvd[end]) (target $(atarg_zvd), halfway park $(spacings_zvd[iimps_zvd[2]]), min over the run $(minimum(spacings_zvd)))")
        println("Deviation from the closed form: $(quaderr_zvd); residual ringing after the last impulse (peak-to-peak): $(ringing_zvd)")

        @assert quaderr_zvd < 1e-5 "the sampled staircase does not reproduce the closed-form response to 1e-5"
        @assert ringing_zvd < 1e-9 "the spacing still rings by $(ringing_zvd), so a catch is landing at the wrong phase"
        # atol, not the default rtol of ~1.5e-8: the halfway park carries about half the O(h^2)
        # quadrature error bounded by quaderr_zvd above, so it lands ~8e-8 off the exact value
        @assert isapprox(spacings_zvd[iimps_zvd[2]],a0_zvd + da_zvd/2; atol=1e-5) "the first catch parks at $(spacings_zvd[iimps_zvd[2]]) rather than halfway to the target at $(a0_zvd + da_zvd/2)"
        # the monotone rise over both segments is what makes the t=0 RK4 step valid for the whole run
        @assert minimum(spacings_zvd) ≈ a0_zvd "the spacing dips below a0, so max(U) is not U(a0) and the critical dt no longer bounds the run"
        @assert minimum(diff(spacings_zvd[1:settled_zvd])) > -1e-15 "the rise to the target is not monotonic"
    end

    # time evolution along the staircase, tracking the instantaneous spectrum
    if false || if_all
        # work on a copy: timeham writes the ramped values, the accumulated integral and the rebuilt
        # U back into the dict, which would otherwise leave the starting-state parameters overwritten
        hamilt_params_tevo_zvd = copy(hamilt_params_zvd)
        hamilt_params_tevo_zvd["scaling_type"] = "magnetic_gradient"

        starting_states_zvd = [Vector{ComplexF64}(states_starting_zvd[i]) for i in 1:speccount_zvd]
        tevo_params_zvd = Dict([ ("magnetic_gradient_time",(pulse_ramp,tmax_zvd,bpulse_zvd)),("trap_frequency",omega_zvd),("tmax",tmax_zvd),("dt",dt_zvd) ])
        tevo_data_zvd,tevo_dict_zvd,instdata_zvd,saving_args_zvd = run_timeevo(starting_states_zvd,tevo_params_zvd,lattice_params_zvd,hamilt_params_tevo_zvd; nev=speccount_zvd,output_level=1,if_instant_gs=true,if_save_data=false,dataloc="tevo-daily-things-data/")

        # the run must end on the static profile the target spacing gives -- comparing against atarg
        # rather than against integral_zvd[end] is the point of the construction
        us_final_zvd = long_range_scaling(ly-1,ly,intstren_zvd; scaling="magnetic_gradient",magnetic_spacing=a0_zvd,magnetic_gradient_integral=da_zvd)
        println("Final U profile: $(hamilt_params_tevo_zvd["U"]) (expected $(us_final_zvd))")
        @assert hamilt_params_tevo_zvd["U"] ≈ us_final_zvd "final U is not the static profile at the target spacing"

        # end-1 skips the final transported column, which is allocated but never written; which
        # vector ED returns for the near-degenerate partner is not reproducible run to run, so track
        # the population of the final instantaneous manifold instead of a state-to-state fidelity
        final_state_zvd = Vector{ComplexF64}(tevo_data_zvd[1][1][:,end-1])
        manifold_population_zvd = sum(abs2(dot(final_state_zvd,Vector{ComplexF64}(instdata_zvd[1][string(i)][:,end]))) for i in 1:speccount_zvd)
        println("Starting energies: $(nrgs_starting_zvd)")
        println("Final instantaneous energies: $([instdata_zvd[2][string(i)][end] for i in 1:speccount_zvd])")
        println("Population of the final instantaneous manifold by the lowest transported state: $(manifold_population_zvd)")
    end

    # the staircase and the spacing it accumulates
    if false || if_all
        figure()
        plot(bgrid_zvd,bpulse_zvd,c="b",label="gradient B")
        plot(bgrid_zvd,spacings_zvd,c="r",label="spacing a = a0 + ∫B(t')sin(ω(t-t'))dt'/ω")
        axhline(atarg_zvd,ls=":",c="r",label="target a = $(atarg_zvd)")
        for (i,t) in enumerate(timps_zvd[2:end])
            axvline(t,ls="--",c="k",label=(i==1 ? "impulses at T/2, T" : nothing))
        end
        legend()
        xlabel("Time")
        ylabel("Value")
        title("ZVD staircase Bf=$(bf_zvd)*[1/4,3/4,1], ω=$(omega_zvd), ringing $(round(ringing_zvd,sigdigits=2))")
    end

    # instantaneous vs transported energies: past the last impulse the spacing, and so H, is time
    # independent, and both series have to be flat. This is what separates "settled" from
    # "correct on average"
    if false || if_all
        times_zvd = [k*dt_zvd for k in 1:length(instdata_zvd[2]["1"])]
        settled_start_zvd = findfirst(t -> t >= timps_zvd[end], times_zvd)
        println("Energy drift after the last impulse: instantaneous $(maximum(abs.(instdata_zvd[2]["1"][settled_start_zvd:end] .- instdata_zvd[2]["1"][end]))), transported $(maximum(abs.(tevo_data_zvd[2][1][settled_start_zvd:end-1] .- tevo_data_zvd[2][1][end-1])))")

        figure()
        for (i,col) in enumerate(["b","g","r"][1:speccount_zvd])
            plot(times_zvd,instdata_zvd[2][string(i)],"-p",c=col,label="E$(i) instantaneous")
            plot(times_zvd,tevo_data_zvd[2][i][1:end-1],c="k",marker="x",label=(i==1 ? "transported" : nothing))
        end
        for t in timps_zvd[2:end]; axvline(t,ls="--",c="k"); end
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time, ZVD pulse $(lx)x$(ly) N=$(n) a $(a0_zvd)→$(atarg_zvd), manifold population = $(round(manifold_population_zvd,digits=6))")
    end

    # real-space density: get_occupancy resolves the density on (physical x, synthetic index m)
    # only, and it is the gradient that gives the synthetic direction a real extent -- the
    # displacement the kernel accumulates is the ladder spacing a(t) entering the dipolar tail, so
    # row m sits at physical y_m(t) = (m-1)*a(t). The ZVD ladder therefore fans out in two stages,
    # pausing at half the target displacement over T/2 <= t <= T, and then stops dead
    if false || if_all
        # column k is the state after k RK4 steps, t = k*dt, and the last column is never written,
        # hence end-1; the starting state is prepended so the series opens at t = 0
        states_py_zvd = vcat([Vector{ComplexF64}(states_starting_zvd[1])],[Vector{ComplexF64}(tevo_data_zvd[1][1][:,k]) for k in 1:size(tevo_data_zvd[1][1],2)-1])
        times_py_zvd = [(k-1)*dt_zvd for k in 1:length(states_py_zvd)]
        occs_py_zvd = [get_occupancy(s,lattice_params_zvd; if_plot=false) for s in states_py_zvd]
        @assert all(occ -> sum(occ) ≈ n, occs_py_zvd) "occupancies do not sum to the particle number, so the transported states are not normalised"

        # cell edges in time (midpoints between steps) and the spacing evaluated there, so the mesh
        # carries one more edge than it has cells; spacings_zvd lives on the half-step grid, so time
        # t sits at index 2t/dt + 1
        tedges_py_zvd = vcat(0.0,0.5 .* (times_py_zvd[1:end-1] .+ times_py_zvd[2:end]),times_py_zvd[end])
        aat_py_zvd = t -> spacings_zvd[min(Int(round(2*t/dt_zvd))+1,length(spacings_zvd))]
        @assert aat_py_zvd(0.0) ≈ a0_zvd "the spacing at t=0 is not the initial spacing the starting state was found at"
        println("Physical y of the top synthetic row: $((ly-1)*aat_py_zvd(0.0)) -> $((ly-1)*aat_py_zvd(times_py_zvd[end]))")

        figure()
        pcolormesh([tedges_py_zvd[k] for j in 1:ly+1, k in 1:length(tedges_py_zvd)],
                   [(j-1.5)*aat_py_zvd(tedges_py_zvd[k]) for j in 1:ly+1, k in 1:length(tedges_py_zvd)],
                   reduce(hcat,[vec(sum(occ,dims=2)) for occ in occs_py_zvd]))
        colorbar(label="Density (summed over physical x)")
        for t in timps_zvd[2:end]; axvline(t,ls="--",c="w"); end
        xlabel("Time")
        ylabel("Physical y")
        title("Real-space y density of the transported gs $(lx)x$(ly) N=$(n), ZVD a $(a0_zvd)→$(atarg_zvd)")
    end

end=#



### Ramp from strongly interacting dd state to weak dd state by ramping the magnetic spacing
# Baseline for a QuOCS optimization of the magnetic_spacing ramp: the spacing analogue of the
# interaction strength ramp blocks above. Instead of lowering the overall strength, the dd
# profile U_r = intstren/(r*a)^3 is weakened by pulling the synthetic states apart, a 0.5 -> 2.0
# at fixed intstren = 10, i.e. nearest-neighbour coupling 80 -> 1.25. Both endpoint manifolds are
# dd groundstates, so an adiabatic ramp can reach fidelity 1.
# The ramp is specified through the gradient rather than the spacing: B(t) is the ZVD staircase
# of the three-step block above, and a(t) is its closed-form response, handed to the evolution as
# a magnetic_spacing pulse_ramp. Going this way round keeps a(t) physical -- a finite gradient
# can only start the mode from rest, so a(t) leaves a0 with zero slope, whereas the linear ramp
# this block used first (fidelity 0.9253 at T = 1) needed delta kicks at both ends.
# dt is pinned at 0.005 (the intstren optimizations' value) rather than left to get_critical_dt,
# so a later dCRAB pulse lives on the same half-step grid as config_intstrenRamp.py. That is only
#= safe while a(t) >= a_start: U grows as 1/a^3, and dt_crit comes from the t = 0 Hamiltonian.
if false

    if_all::Bool = true

    # starting and ending states: dd groundstate manifolds at the two spacings
    if false || if_all
        lx,ly,n = 4,4,2
        intstren_ms = 10.0
        a_start_ms, a_end_ms = 0.5, 2.0
        omega_ms = 10.0         # trap frequency of the displaced mode, as in the posicast/ZVD blocks
        speccount_ms = 2

        pdict_ms = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("scaling_type","dd"),("trap_frequency",omega_ms),("magnetic_spacing",a_start_ms),("interaction_strength",intstren_ms),("filling",0.5),("nev",speccount_ms+2),("if_find_data",false),("if_save_data",false)])
        _,_,_,_,_,lattice_params_ms,hamilt_params_ms = run_normal_ed(pdict_ms; output_level=0)

        pdict_ending_ms = merge(pdict_ms,Dict("magnetic_spacing"=>a_end_ms))
        _,_,_,_,_,_,hamilt_params_ending_ms = run_normal_ed(pdict_ending_ms; output_level=0)

        # both endpoints have a 2-fold degenerate groundstate, and a fresh Lanczos draw sometimes
        # returns only one of the pair (seen at a = 2.0: target "manifold" = gs + first excited,
        # fidelity capped near 0.5). The 120-dim space is tiny, so diagonalize densely instead
        using LinearAlgebra
        eig_starting_ms = eigen(Hermitian(Matrix(hamilt_params_ms["H"])))
        eig_ending_ms = eigen(Hermitian(Matrix(hamilt_params_ending_ms["H"])))
        nrgs_starting_ms, nrgs_ending_ms = eig_starting_ms.values[1:speccount_ms+2], eig_ending_ms.values[1:speccount_ms+2]
        states_starting_ms = [eig_starting_ms.vectors[:,i] for i in 1:speccount_ms+2]
        states_ending_ms = [eig_ending_ms.vectors[:,i] for i in 1:speccount_ms+2]
        for nrgs in (nrgs_starting_ms,nrgs_ending_ms)
            @assert nrgs[speccount_ms] - nrgs[1] < 1e-8 && nrgs[speccount_ms+1] - nrgs[speccount_ms] > 1e-3 "groundstate manifold is not an isolated $(speccount_ms)-fold degenerate level: $(nrgs)"
        end

        println("Starting U profile: $(hamilt_params_ms["U"])")
        println("Starting energies: $(nrgs_starting_ms)")
        println("Ending energies: $(nrgs_ending_ms)")

        # sudden-quench floor: what the ramp has to beat
        reference_fidelity_ms = real(groundstate_manifold_fidelity(states_ending_ms[1:speccount_ms],states_starting_ms[1:speccount_ms]))
        println("Sudden-quench fidelity a $(a_start_ms)→$(a_end_ms): $(reference_fidelity_ms)")
    end

    # ZVD gradient pulse and the spacing it drives. x = a - a0 obeys xddot + w^2 x = B(t) from
    # rest; the staircase Bf*[1/4,3/4,1] switching at 0, pi/w, 2pi/w is the impulse train
    # (1/4,1/2,1/4) convolved with a step, so x is the matching sum of shifted step responses,
    # x(t) = da * sum_{t_i <= t} A_i (1 - cos w(t - t_i)), which parks at da from 2pi/w onward.
    # a(t) is sampled from that closed form, not integrated from B, so the switches need not sit
    # on the half-step grid (a is C1 across them; only addot jumps)
    if false || if_all
        ramptime_ms = 1.0
        dt_ms = 0.005

        da_ms = a_end_ms - a_start_ms
        bf_ms = omega_ms^2 * da_ms                  # final standing gradient, = peak B
        amps_ms = [0.25, 0.5, 0.25]                 # ZVD impulse amplitudes, as fractions of Bf
        timps_ms = [0.0, 1.0, 2.0] .* (pi/omega_ms)
        @assert timps_ms[end] <= ramptime_ms "ZVD pulse settles at $(timps_ms[end]), after the end of the ramp $(ramptime_ms)"

        gradient_ms(t) = bf_ms * sum((A for (A,ti) in zip(amps_ms,timps_ms) if ti <= t); init=0.0)
        spacing_ms(t) = a_start_ms + da_ms * sum((A*(1 - cos(omega_ms*(t - ti))) for (A,ti) in zip(amps_ms,timps_ms) if ti <= t); init=0.0)

        # same half-step grid and length pulse_ramp expects: ceil(ramptime/(dt/2)) + 1 samples
        h_ms = dt_ms / 2
        times_a_ms = (0:Int(ceil(ramptime_ms/h_ms))) .* h_ms
        a_ms = spacing_ms.(times_a_ms)

        # a starts exactly at a0 and only rises, so the t = 0 RK4 step bounds the whole run, and
        # from the last impulse on it sits at a_end (the cos terms cancel identically)
        @assert minimum(a_ms) == a_start_ms "spacing dips below a_start, the pinned dt is no longer stable"
        settled_ms = times_a_ms .>= timps_ms[end]
        @assert maximum(abs, a_ms[settled_ms] .- a_end_ms) < 1e-12 "spacing has not settled at a_end after the last impulse"

        # check the closed form against the equation of motion by inverting it back to B, at every
        # sample whose stencil does not straddle a switch (O(h^2 w^4 da) truncation error, hence
        # the explicit tolerance), and that it starts from rest, i.e. needs no kick at t = 0
        b_back_ms, kick_back_ms = get_magnetic_gradient_from_spacing(a_ms,h_ms,omega_ms)
        smooth_ms = [i for i in 1:length(a_ms) if !any(ti -> abs(times_a_ms[i] - ti) < 1.5h_ms, timps_ms[2:end])]
        err_b_back_ms = maximum(abs, b_back_ms[smooth_ms] .- gradient_ms.(times_a_ms[smooth_ms]))
        println("ZVD spacing ramp: settles at t=$(timps_ms[end]), B recovered from a(t) to $(err_b_back_ms) (Bf = $(bf_ms)), kick at t=0 $(kick_back_ms)")
        @assert err_b_back_ms < 1e-3*bf_ms "a(t) does not satisfy xddot + w^2 x = B(t)"
        @assert abs(kick_back_ms) < 1e-3*da_ms "a(t) leaves a0 with a slope, which no finite gradient can produce"
    end

    # ZVD spacing ramp with instantaneous energies
    if false || if_all
        # work on a copy: run_timeevo's timeham writes the ramp's current value back into
        # the dict, which would leave magnetic_spacing=a_end for any later section
        hamilt_params_energy_ms = copy(hamilt_params_ms)

        time_running_args_ms = (nev=speccount_ms,output_level=1,if_instant_gs=true,if_save_data=false,dataloc="tevo-daily-things-data/")
        starting_states_ms = [Vector{ComplexF64}(states_starting_ms[i]) for i in 1:speccount_ms]
        tevo_params_ms = Dict([ ("magnetic_spacing",(pulse_ramp,ramptime_ms,a_ms)),("tmax",ramptime_ms),("dt",dt_ms) ])
        tevo_data_ms,tevo_dict_ms,instdata_ms,saving_args_ms = run_timeevo(starting_states_ms,tevo_params_ms,lattice_params_ms,hamilt_params_energy_ms; time_running_args_ms...)

        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        final_manifold_ms = [Vector{ComplexF64}(tevo_data_ms[1][i][:,end-1]) for i in 1:speccount_ms]
        fidelity_ms = real(groundstate_manifold_fidelity(final_manifold_ms,[Vector{ComplexF64}(s) for s in states_ending_ms[1:speccount_ms]]))
        println("Fidelity with target manifold using ZVD spacing ramp T=$(ramptime_ms): $(fidelity_ms)")
    end

    # plot the instantaneous vs transported energies along the ramp
    if true || if_all
        times_ms = range(0.0,ramptime_ms,length=length(instdata_ms[2]["1"]))

        figure()
        cols = ["b","g","r"]
        for i in 1:speccount_ms
            plot(times_ms,instdata_ms[2][string(i)],"-p",c=cols[i],label="E$(i) instantaneous")
            plot(times_ms,tevo_data_ms[2][i][1:end-1],c="k",marker="x",label=(i==1 ? "transported" : nothing))
        end
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time for ZVD spacing ramp $(lx)x$(ly) N=$(n) U=$(intstren_ms) T=$(ramptime_ms) a $(a_start_ms)→$(a_end_ms), fidelity = $(round(fidelity_ms,digits=6))")
    end

    # plot the applied gradient and the spacing it drives
    if true || if_all
        fig, axs = subplots(2,1,sharex=true)
        levels_ms = bf_ms .* cumsum(amps_ms)
        axs[1].step(vcat(timps_ms,ramptime_ms),vcat(levels_ms,levels_ms[end]),where="post",c="b")
        axs[1].axhline(bf_ms,ls="--",c="gray",label="holding gradient w²(a_end-a0)")
        axs[1].set_ylabel("Magnetic gradient B(t)")
        axs[1].set_title("ZVD gradient pulse and the spacing it drives, w=$(omega_ms)")
        axs[1].legend()
        axs[1].set_ylim(0.0,1.1*bf_ms)
        axs[2].plot(times_a_ms,a_ms,c="k")
        for ti in timps_ms
            axs[2].axvline(ti,ls=":",c="gray")
        end
        axs[2].set_xlabel("Time")
        axs[2].set_ylabel("Spacing a(t)")
        axs[2].set_ylim(0.0,1.1*a_end_ms)
    end

end=#



### Time evolution with the QuOCS dCRAB-optimized magnetic spacing ramp from optimal-control/config_magspacRamp.py
# All parameters must match the ones the optimization ran with (see config_magspacRamp.py):
# 4x4 N=2 pbc, dd intstren 10, a 0.5 -> 2.0 over ramptime 1.0, pulse on the dt = 0.005 half-step
# grid, a >= 0.1. The initial guess was the ZVD spacing ramp of the block above, and the update
# was scaled by 16 (t/T)^2 (1 - t/T)^2, which pins value and slope of a(t) at both ends, so the
# optimized a(t) still leaves a_start and reaches a_end at rest: B(t) stays finite, no kicks.
# The fidelity is computed twice: with the fast split-operator evolution the optimizer used
# (magspac-ramp-control-functions.jl), and independently with run_timeevo's RK4, which is only
# stable down to the smallest spacing the pulse reaches -- so dt is refined by an integer factor
# from dt_crit at min(a), with the pulse linearly interpolated onto the finer half-step grid
#= (the fast evolution also treats a(t) as linear between samples).
if false

    if_all::Bool = true

    include("magspac-ramp-control-functions.jl")

    # endpoint manifolds and fast-evolution setup, same parameters as the optimization
    if false || if_all
        lx,ly,n = 4,4,2
        speccount_msopt = 2
        intstren_msopt = 10.0
        a_start_msopt, a_end_msopt = 0.5, 2.0
        ramptime_msopt = 1.0
        dt_msopt = 0.005
        omega_msopt = 10.0      # of the ZVD guess, and of the gradient read back from a(t)

        params_msopt = Dict("Lx"=>lx,"Ly"=>ly,"N"=>n,"speccount"=>speccount_msopt,"intstren"=>intstren_msopt,"a_start"=>a_start_msopt,"a_end"=>a_end_msopt,"ramptime"=>ramptime_msopt,"dt"=>dt_msopt)
        setup_msopt = setup_magspac_ramp(params_msopt)
    end

    # load the optimized pulse; fidelities of it and of the ZVD guess with the fast evolution
    if false || if_all
        quocs_folder_msopt = "../optimal-control/QuOCS_Results/20260924_150154_magspacRamp_dCRAB"
        controls_file_msopt = filter(f -> endswith(f,"best_controls.npz"), readdir(quocs_folder_msopt))[1]
        # only read the numeric arrays: NPZ.jl cannot parse the numpy unicode-string arrays
        # (pulse_names etc.) that QuOCS also stores in the file
        best_controls_msopt = npzread(joinpath(quocs_folder_msopt,controls_file_msopt),["magspacRamp","time_grid_for_magspacRamp"])
        a_msopt = Float64.(real.(best_controls_msopt["magspacRamp"]))
        times_a_msopt = Float64.(real.(best_controls_msopt["time_grid_for_magspacRamp"]))
        a_guess_msopt = zvd_spacing_pulse(a_start_msopt,a_end_msopt,omega_msopt,ramptime_msopt,dt_msopt)

        fidelity_msopt = compute_fidelity_magspac_ramp([a_msopt],params_msopt,setup_msopt)
        fidelity_guess_msopt = compute_fidelity_magspac_ramp([a_guess_msopt],params_msopt,setup_msopt)
        println("Fast evolution: optimized fidelity $(fidelity_msopt), ZVD guess $(fidelity_guess_msopt), min a $(minimum(a_msopt))")
    end

    # independent check with run_timeevo, with instantaneous energies for the plot below
    if false || if_all
        pdict_msopt = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("scaling_type","dd"),("trap_frequency",omega_msopt),("magnetic_spacing",a_start_msopt),("interaction_strength",intstren_msopt),("filling",0.5),("nev",speccount_msopt),("if_find_data",false),("if_save_data",false)])
        _,_,_,_,_,lattice_params_msopt,hamilt_params_msopt = run_normal_ed(pdict_msopt; output_level=0)

        # RK4 step for the strongest couplings the pulse reaches, as an integer refinement of dt
        u_at_min_a_msopt = long_range_scaling(length(hamilt_params_msopt["U"])-1,ly,intstren_msopt; scaling="dd",magnetic_spacing=minimum(a_msopt))
        dt_crit_msopt = get_critical_dt(ramptime_msopt,lattice_params_msopt,Dict("U"=>u_at_min_a_msopt))
        refine_msopt = max(1,ceil(Int,dt_msopt/dt_crit_msopt))
        dt_run_msopt = dt_msopt / refine_msopt

        # the pulse on the finer half-step grid, at the length pulse_ramp expects
        times_run_msopt = (0:Int(ceil(ramptime_msopt/(dt_run_msopt/2)))) .* (dt_run_msopt/2)
        a_run_msopt = map(times_run_msopt) do t
            i = clamp(searchsortedlast(times_a_msopt,t),1,length(times_a_msopt)-1)
            w = (t - times_a_msopt[i]) / (times_a_msopt[i+1] - times_a_msopt[i])
            (1 - w)*a_msopt[i] + w*a_msopt[i+1]
        end
        println("run_timeevo: min a $(minimum(a_msopt)) -> dt_crit $(dt_crit_msopt), dt refined x$(refine_msopt) to $(dt_run_msopt)")

        # work on a copy: run_timeevo's timeham writes the ramp's current value back into
        # the dict, which would leave magnetic_spacing=a_end for any later section
        hamilt_params_energy_msopt = copy(hamilt_params_msopt)
        starting_states_msopt = setup_msopt[1]
        time_running_args_msopt = (nev=speccount_msopt,output_level=0,if_instant_gs=true,if_save_data=false,dataloc="tevo-daily-things-data/")
        tevo_params_msopt = Dict([ ("magnetic_spacing",(pulse_ramp,ramptime_msopt,a_run_msopt)),("tmax",ramptime_msopt),("dt",dt_run_msopt) ])
        tevo_data_msopt,tevo_dict_msopt,instdata_msopt,_ = run_timeevo(starting_states_msopt,tevo_params_msopt,lattice_params_msopt,hamilt_params_energy_msopt; time_running_args_msopt...)

        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        final_manifold_msopt = [Vector{ComplexF64}(tevo_data_msopt[1][i][:,end-1]) for i in 1:speccount_msopt]
        fidelity_rk4_msopt = real(groundstate_manifold_fidelity(final_manifold_msopt,setup_msopt[2]))
        println("Fidelity with target manifold using the optimized spacing ramp: $(fidelity_rk4_msopt) (run_timeevo), $(fidelity_msopt) (fast evolution)")
    end

    # plot the instantaneous vs transported energies along the optimized ramp
    if true || if_all
        times_msopt = range(0.0,ramptime_msopt,length=length(instdata_msopt[2]["1"]))

        figure()
        cols = ["b","g","r"]
        for i in 1:speccount_msopt
            plot(times_msopt,instdata_msopt[2][string(i)],"-p",c=cols[i],label="E$(i) instantaneous")
            plot(times_msopt,tevo_data_msopt[2][i][1:end-1],c="k",marker="x",label=(i==1 ? "transported" : nothing))
        end
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time for dCRAB spacing ramp $(lx)x$(ly) N=$(n) a $(a_start_msopt)→$(a_end_msopt), fidelity = $(round(fidelity_rk4_msopt,digits=6))")
    end

    # plot the optimized a(t) against the ZVD guess, and the gradient each needs
    if true || if_all
        h_msopt = dt_msopt / 2
        b_msopt, kick_msopt = get_magnetic_gradient_from_spacing(a_msopt,h_msopt,omega_msopt)
        b_guess_msopt, kick_guess_msopt = get_magnetic_gradient_from_spacing(a_guess_msopt,h_msopt,omega_msopt)
        println("Gradient read back from a(t): optimized max|B| $(maximum(abs,b_msopt)), kick at t=0 $(kick_msopt); ZVD guess max|B| $(maximum(abs,b_guess_msopt)), kick $(kick_guess_msopt)")

        fig, axs = subplots(2,1,sharex=true)
        axs[1].plot(times_a_msopt,a_guess_msopt,"--",c="gray",label="ZVD guess, F = $(round(fidelity_guess_msopt,digits=4))")
        axs[1].plot(times_a_msopt,a_msopt,c="k",label="dCRAB, F = $(round(fidelity_msopt,digits=4))")
        axs[1].set_ylabel("Spacing a(t)")
        axs[1].set_title("dCRAB-optimized spacing ramp and its gradient, w=$(omega_msopt)")
        axs[1].legend()
        axs[2].plot(times_a_msopt,b_guess_msopt,"--",c="gray",label="ZVD guess")
        axs[2].plot(times_a_msopt,b_msopt,c="b",label="dCRAB")
        axs[2].set_xlabel("Time")
        axs[2].set_ylabel("Magnetic gradient B(t)")
        axs[2].legend()
    end

end=#


### Time evolution with the QuOCS dCRAB-optimized magnetic GRADIENT from optimal-control/config_magspacGrad.py
# Same transport as the spacing-ramp blocks above (4x4 N=2 pbc, dd intstren 10, a 0.5 -> 2.0), but
# over the ZVD pulse's own duration T = 2pi/w rather than 1.0, on a grid of exactly 256 half-steps
# (dt = 2T/256, so both ZVD switches sit on samples), and the control is the gradient B(t) itself, with the
# hard cap |B| <= 300 (2x the ZVD pulse). Optimizing a(t) directly (block above) reached F = 0.939
# only by shaking the lattice, read-back |B| ~ 3300; here a(t) is the trapped mode's response to B,
#     a(t) = a_start + (1/w) int_0^t B(t') sin(w(t - t')) dt',
# so every pulse is physical by construction. The end condition comes from the optimization: B(T)
# is pinned at the holding gradient w^2 (a_end - a_start) = 150, and the FoM subtracts the amplitude
# R = sqrt((a(T) - a_end)^2 + (adot(T)/w)^2) the mode would ring with after T.
# The fidelity is computed twice: with the fast split-operator evolution of a(t) the optimizer used,
# and independently with run_timeevo driven through the "magnetic_gradient" scaling, i.e. B goes
# in and timeham accumulates the spacing itself via get_magnetic_gradient_integral
#= (dt refined from dt_crit at min(a) if needed, with B linearly interpolated onto the finer grid).
if false

    if_all::Bool = true

    include("magspac-ramp-control-functions.jl")

    # endpoint manifolds and fast-evolution setup, same parameters as the optimization
    if false || if_all
        lx,ly,n = 4,4,2
        speccount_mg = 2
        intstren_mg = 10.0
        a_start_mg, a_end_mg = 0.5, 2.0
        omega_mg = 10.0
        ramptime_mg = 2pi/omega_mg      # the ZVD pulse's settling time
        dt_mg = 2*(ramptime_mg/256)     # power of two: ramptime/(dt/2) is exactly 256.0, so no stray sample
        bmax_mg = 300.0

        params_mg = Dict("Lx"=>lx,"Ly"=>ly,"N"=>n,"speccount"=>speccount_mg,"intstren"=>intstren_mg,"a_start"=>a_start_mg,"a_end"=>a_end_mg,"ramptime"=>ramptime_mg,"dt"=>dt_mg,"trap_frequency"=>omega_mg,"a_floor"=>0.1)
        setup_mg = setup_magspac_ramp(params_mg)
    end

    # load the optimized gradient, the spacing it drives, and the FoM terms of it and of the ZVD guess
    if false || if_all
        quocs_folder_mg = "../optimal-control/QuOCS_Results/20260925_102851_magspacGrad_dCRAB"
        controls_file_mg = filter(f -> endswith(f,"best_controls.npz"), readdir(quocs_folder_mg))[1]
        # only read the numeric arrays: NPZ.jl cannot parse the numpy unicode-string arrays
        # (pulse_names etc.) that QuOCS also stores in the file
        best_controls_mg = npzread(joinpath(quocs_folder_mg,controls_file_mg),["magspacGrad","time_grid_for_magspacGrad"])
        b_mg = Float64.(real.(best_controls_mg["magspacGrad"]))
        times_b_mg = Float64.(real.(best_controls_mg["time_grid_for_magspacGrad"]))
        b_guess_mg = zvd_gradient_pulse(a_start_mg,a_end_mg,omega_mg,ramptime_mg,dt_mg)

        h_mg = dt_mg / 2
        a_mg = a_start_mg .+ magnetic_gradient_response(b_mg,h_mg,omega_mg)[1]
        a_guess_mg = a_start_mg .+ magnetic_gradient_response(b_guess_mg,h_mg,omega_mg)[1]
        fom_mg, fidelity_mg, ringing_mg, min_a_mg = magspac_gradient_fom(b_mg,params_mg,setup_mg)
        _, fidelity_guess_mg, ringing_guess_mg, _ = magspac_gradient_fom(b_guess_mg,params_mg,setup_mg)
        println("Fast evolution: optimized fidelity $(fidelity_mg), ringing $(ringing_mg) (FoM $(fom_mg)); ZVD guess fidelity $(fidelity_guess_mg), ringing $(ringing_guess_mg)")
        println("Optimized gradient: max|B| $(maximum(abs,b_mg)) (cap $(bmax_mg)), B(0) $(b_mg[1]), B(T) $(b_mg[end]); a in [$(minimum(a_mg)), $(maximum(a_mg))], a(T) $(a_mg[end])")
        @assert maximum(abs,b_mg) <= bmax_mg "the loaded gradient exceeds the cap it was optimized under"
        @assert b_mg[end] ≈ omega_mg^2*(a_end_mg - a_start_mg) "B(T) is not the holding gradient"

        # the O(n) response against get_magnetic_gradient_integral itself, the quadrature run_timeevo
        # uses below, at a handful of samples (uniform step, so when_dt_ends goes past the end)
        checkparams_mg = Dict{String,Any}([("magnetic_gradient_time",b_mg),("when_dt_ends",[length(b_mg)+1,length(b_mg)+1]),("dt",[dt_mg,dt_mg]),("trap_frequency",omega_mg)])
        for k in (1,2,length(b_mg)÷3,length(b_mg)÷2,length(b_mg))
            @assert isapprox(get_magnetic_gradient_integral(k,checkparams_mg),a_mg[k] - a_start_mg; atol=1e-12) "the O(n) kernel split disagrees with get_magnetic_gradient_integral at sample $k"
        end
    end

    # independent check with run_timeevo driven by the gradient, with instantaneous energies for the plot
    if false || if_all
        pdict_mg = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("scaling_type","dd"),("trap_frequency",omega_mg),("magnetic_spacing",a_start_mg),("interaction_strength",intstren_mg),("filling",0.5),("nev",speccount_mg),("if_find_data",false),("if_save_data",false)])
        _,_,_,_,_,lattice_params_mg,hamilt_params_mg = run_normal_ed(pdict_mg; output_level=0)
        @assert long_range_scaling(length(hamilt_params_mg["U"])-1,ly,intstren_mg; scaling="magnetic_gradient",magnetic_spacing=a_start_mg,magnetic_gradient_integral=0.0) ≈ hamilt_params_mg["U"] "t=0 gradient profile does not match the dd profile the starting state was found with"

        # RK4 step for the strongest couplings the spacing reaches, as an integer refinement of dt
        u_at_min_a_mg = long_range_scaling(length(hamilt_params_mg["U"])-1,ly,intstren_mg; scaling="dd",magnetic_spacing=minimum(a_mg))
        dt_crit_mg = get_critical_dt(ramptime_mg,lattice_params_mg,Dict("U"=>u_at_min_a_mg))
        refine_mg = max(1,ceil(Int,dt_mg/dt_crit_mg))
        dt_run_mg = dt_mg / refine_mg

        # the gradient on the finer half-step grid, at the length pulse_ramp expects
        times_run_mg = (0:Int(ceil(ramptime_mg/(dt_run_mg/2)))) .* (dt_run_mg/2)
        b_run_mg = map(times_run_mg) do t
            i = clamp(searchsortedlast(times_b_mg,t),1,length(times_b_mg)-1)
            w = (t - times_b_mg[i]) / (times_b_mg[i+1] - times_b_mg[i])
            (1 - w)*b_mg[i] + w*b_mg[i+1]
        end
        println("run_timeevo: min a $(minimum(a_mg)) -> dt_crit $(dt_crit_mg), dt refined x$(refine_mg) to $(dt_run_mg)")

        # work on a copy: timeham writes the ramped values, the accumulated integral and the rebuilt
        # U back into the dict, which would otherwise leave the starting-state parameters overwritten
        hamilt_params_tevo_mg = copy(hamilt_params_mg)
        hamilt_params_tevo_mg["scaling_type"] = "magnetic_gradient"
        starting_states_mg = setup_mg[1]
        time_running_args_mg = (nev=speccount_mg,output_level=0,if_instant_gs=true,if_save_data=false,dataloc="tevo-daily-things-data/")
        tevo_params_mg = Dict([ ("magnetic_gradient_time",(pulse_ramp,ramptime_mg,b_run_mg)),("trap_frequency",omega_mg),("tmax",ramptime_mg),("dt",dt_run_mg) ])
        tevo_data_mg,tevo_dict_mg,instdata_mg,_ = run_timeevo(starting_states_mg,tevo_params_mg,lattice_params_mg,hamilt_params_tevo_mg; time_running_args_mg...)

        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        final_manifold_mg = [Vector{ComplexF64}(tevo_data_mg[1][i][:,end-1]) for i in 1:speccount_mg]
        fidelity_rk4_mg = real(groundstate_manifold_fidelity(final_manifold_mg,setup_mg[2]))
        println("Final U profile: $(hamilt_params_tevo_mg["U"])")
        println("Fidelity with target manifold using the optimized gradient: $(fidelity_rk4_mg) (run_timeevo), $(fidelity_mg) (fast evolution)")
    end

    # plot the instantaneous vs transported energies along the optimized ramp
    if true || if_all
        times_mg = range(0.0,ramptime_mg,length=length(instdata_mg[2]["1"]))

        figure()
        cols = ["b","g","r"]
        for i in 1:speccount_mg
            plot(times_mg,instdata_mg[2][string(i)],"-p",c=cols[i],label="E$(i) instantaneous")
            plot(times_mg,tevo_data_mg[2][i][1:end-1],c="k",marker="x",label=(i==1 ? "transported" : nothing))
        end
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time for dCRAB gradient $(lx)x$(ly) N=$(n) a $(a_start_mg)→$(a_end_mg), |B|≤$(bmax_mg), fidelity = $(round(fidelity_rk4_mg,digits=6))")
    end

    # plot the optimized gradient against the ZVD guess, and the spacing each drives
    if true || if_all
        fig, axs = subplots(2,1,sharex=true)
        axs[1].plot(times_b_mg,b_guess_mg,"--",c="gray",label="ZVD guess")
        axs[1].plot(times_b_mg,b_mg,c="b",label="dCRAB")
        for bl in (-bmax_mg,bmax_mg)
            axs[1].axhline(bl,ls=":",c="r")
        end
        axs[1].axhline(omega_mg^2*(a_end_mg - a_start_mg),ls="--",c="k",lw=0.8,label="holding gradient w²(a_end-a0)")
        axs[1].set_ylabel("Magnetic gradient B(t)")
        axs[1].set_title("dCRAB-optimized gradient (cap $(bmax_mg)) and the spacing it drives, w=$(omega_mg)")
        axs[1].legend()
        axs[2].plot(times_b_mg,a_guess_mg,"--",c="gray",label="ZVD guess, F = $(round(fidelity_guess_mg,digits=4))")
        axs[2].plot(times_b_mg,a_mg,c="k",label="dCRAB, F = $(round(fidelity_mg,digits=4)), R = $(round(ringing_mg,sigdigits=2))")
        axs[2].set_xlabel("Time")
        axs[2].set_ylabel("Spacing a(t)")
        axs[2].legend()
    end

end=#


### Regular CRAB / dCRAB for the magnetic gradient: optimal-control/config_magspacGradCRAB.py, then replay
# Same problem as the dCRAB gradient block above (4x4 N=2 pbc, dd intstren 10, a 0.5 -> 2.0,
# T = 2pi/w on 256 half-steps, |B| <= 300, FoM = fidelity - ringing R), optimizer as advised:
# regular CRAB (one super-iteration, one Nelder-Mead search), 4 randomized principal harmonics
# (one frequency in each [k - 1/2, k + 1/2] cycles per T), and both ends pinned at the gradients that
# hold the lattice at rest, B(0) = 0 and B(T) = w^2 (a_end - a_start) = 150, with a smooth guess
# B = 150 (1 - cos(pi t/T))/2 between them (fidelity 0.9048 but ringing R 0.50).
# The first sub-block runs the optimization in the quocs-env Python for seeds 1-4 (a few s each, writes
# QuOCS_Results/<stamp>_magspacGrad_CRAB_seed<k> folders and local-figs pngs; with si_crab > 1 it is dCRAB
# from the same guess and pins, folders <stamp>_magspacGradSmooth_dCRAB_seed<k>, ~30 s each); the rest
#= replay the best dCRAB run (2026-09-30: 10 SI seed 2, F 0.9208).
if false

    if_all::Bool = true

    include("magspac-ramp-control-functions.jl")

    # run the CRAB optimization once per seed (plain CRAB's result depends on its one frequency draw);
    # the config boots its own Julia and must run from optimal-control/. Off by default: seeds 1-4 exist
    # (2026-09-30, CRAB and 10-SI dCRAB) and every run adds a results folder
    if false
        si_crab = 10    # super-iterations: 1 = plain CRAB
        for seed in 1:4
            run(Cmd(`$(abspath("../optimal-control/quocs-env/bin/python")) config_magspacGradCRAB.py $(seed) $(si_crab)`; dir=abspath("../optimal-control")))
        end
    end

    # endpoint manifolds and fast-evolution setup, same parameters as the optimization
    if false || if_all
        lx,ly,n = 4,4,2
        speccount_crab = 2
        intstren_crab = 10.0
        a_start_crab, a_end_crab = 0.5, 2.0
        omega_crab = 10.0
        ramptime_crab = 2pi/omega_crab
        dt_crab = 2*(ramptime_crab/256)
        bmax_crab = 300.0
        bhold_crab = omega_crab^2*(a_end_crab - a_start_crab)

        params_crab = Dict("Lx"=>lx,"Ly"=>ly,"N"=>n,"speccount"=>speccount_crab,"intstren"=>intstren_crab,"a_start"=>a_start_crab,"a_end"=>a_end_crab,"ramptime"=>ramptime_crab,"dt"=>dt_crab,"trap_frequency"=>omega_crab,"a_floor"=>0.1)
        setup_crab = setup_magspac_ramp(params_crab)
    end

    # compare the dCRAB runs (seeds), keep the best one, and the FoM terms of it and of the smooth guess
    if false || if_all
        results_crab = "../optimal-control/QuOCS_Results"
        # only the dCRAB (10 SI) runs; plain CRAB runs are "_magspacGrad_CRAB"
        runs_crab = sort(filter(f -> occursin("_magspacGradSmooth_dCRAB",f),readdir(results_crab)))
        # only read the numeric arrays: NPZ.jl cannot parse the numpy unicode-string arrays
        # (pulse_names etc.) that QuOCS also stores in the file
        load_crab(folder) = npzread(joinpath(results_crab,folder,filter(f -> endswith(f,"best_controls.npz"),readdir(joinpath(results_crab,folder)))[1]),["magspacGrad","time_grid_for_magspacGrad"])
        foms_crab = [magspac_gradient_fom(Float64.(real.(load_crab(r)["magspacGrad"])),params_crab,setup_crab) for r in runs_crab]
        for (r,f) in zip(runs_crab,foms_crab)
            println("$(r): FoM $(f[1]), fidelity $(f[2]), ringing $(f[3])")
        end
        quocs_folder_crab = runs_crab[argmax(first.(foms_crab))]
        best_controls_crab = load_crab(quocs_folder_crab)
        b_crab = Float64.(real.(best_controls_crab["magspacGrad"]))
        times_b_crab = Float64.(real.(best_controls_crab["time_grid_for_magspacGrad"]))
        b_guess_crab = bhold_crab .* (1 .- cos.(pi .* times_b_crab ./ ramptime_crab)) ./ 2
        method_crab = occursin("dCRAB",quocs_folder_crab) ? "dCRAB" : "CRAB"
        println("Replaying the best run, $(quocs_folder_crab)")

        h_crab = dt_crab / 2
        a_crab = a_start_crab .+ magnetic_gradient_response(b_crab,h_crab,omega_crab)[1]
        a_guess_crab = a_start_crab .+ magnetic_gradient_response(b_guess_crab,h_crab,omega_crab)[1]
        fom_crab, fidelity_crab, ringing_crab, min_a_crab = magspac_gradient_fom(b_crab,params_crab,setup_crab)
        _, fidelity_guess_crab, ringing_guess_crab, _ = magspac_gradient_fom(b_guess_crab,params_crab,setup_crab)
        println("Fast evolution: $(method_crab) fidelity $(fidelity_crab), ringing $(ringing_crab) (FoM $(fom_crab)); smooth guess fidelity $(fidelity_guess_crab), ringing $(ringing_guess_crab)")
        println("$(method_crab) gradient: max|B| $(maximum(abs,b_crab)) (cap $(bmax_crab)), B(0) $(b_crab[1]), B(T) $(b_crab[end]); a in [$(minimum(a_crab)), $(maximum(a_crab))], a(T) $(a_crab[end])")
        @assert maximum(abs,b_crab) <= bmax_crab "the loaded gradient exceeds the cap it was optimized under"
        @assert b_crab[1] == 0.0 && b_crab[end] ≈ bhold_crab "the ends are not pinned at the holding gradients"
    end

    # independent check with run_timeevo driven by the gradient, with instantaneous energies for the plot
    if false || if_all
        pdict_crab = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("scaling_type","dd"),("trap_frequency",omega_crab),("magnetic_spacing",a_start_crab),("interaction_strength",intstren_crab),("filling",0.5),("nev",speccount_crab),("if_find_data",false),("if_save_data",false)])
        _,_,_,_,_,lattice_params_crab,hamilt_params_crab = run_normal_ed(pdict_crab; output_level=0)

        # RK4 step for the strongest couplings the spacing reaches, as an integer refinement of dt
        u_at_min_a_crab = long_range_scaling(length(hamilt_params_crab["U"])-1,ly,intstren_crab; scaling="dd",magnetic_spacing=minimum(a_crab))
        dt_crit_crab = get_critical_dt(ramptime_crab,lattice_params_crab,Dict("U"=>u_at_min_a_crab))
        refine_crab = max(1,ceil(Int,dt_crab/dt_crit_crab))
        dt_run_crab = dt_crab / refine_crab

        # the gradient on the finer half-step grid, at the length pulse_ramp expects
        times_run_crab = (0:Int(ceil(ramptime_crab/(dt_run_crab/2)))) .* (dt_run_crab/2)
        b_run_crab = map(times_run_crab) do t
            i = clamp(searchsortedlast(times_b_crab,t),1,length(times_b_crab)-1)
            w = (t - times_b_crab[i]) / (times_b_crab[i+1] - times_b_crab[i])
            (1 - w)*b_crab[i] + w*b_crab[i+1]
        end
        println("run_timeevo: min a $(minimum(a_crab)) -> dt_crit $(dt_crit_crab), dt refined x$(refine_crab) to $(dt_run_crab)")

        # work on a copy: timeham writes the ramped values back into the dict
        hamilt_params_tevo_crab = copy(hamilt_params_crab)
        hamilt_params_tevo_crab["scaling_type"] = "magnetic_gradient"
        # nev only sets how many instantaneous eigenstates are tracked; the evolved states are the 2 of the manifold
        nevinst_crab = 3
        time_running_args_crab = (nev=nevinst_crab,output_level=0,if_instant_gs=true,if_save_data=false,dataloc="tevo-daily-things-data/")
        tevo_params_crab = Dict([ ("magnetic_gradient_time",(pulse_ramp,ramptime_crab,b_run_crab)),("trap_frequency",omega_crab),("tmax",ramptime_crab),("dt",dt_run_crab) ])
        tevo_data_crab,tevo_dict_crab,instdata_crab,_ = run_timeevo(setup_crab[1],tevo_params_crab,lattice_params_crab,hamilt_params_tevo_crab; time_running_args_crab...)

        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        final_manifold_crab = [Vector{ComplexF64}(tevo_data_crab[1][i][:,end-1]) for i in 1:speccount_crab]
        fidelity_rk4_crab = real(groundstate_manifold_fidelity(final_manifold_crab,setup_crab[2]))
        println("Fidelity with target manifold using the $(method_crab) gradient: $(fidelity_rk4_crab) (run_timeevo), $(fidelity_crab) (fast evolution)")
    end

    # plot the instantaneous energies of the lowest 3 states vs the transported energies of the two
    # evolved manifold states, both under the same best pulse
    if false || if_all
        times_crab = range(0.0,ramptime_crab,length=length(instdata_crab[2]["1"]))

        figure()
        cols = ["b","g","r"]
        for i in 1:nevinst_crab
            plot(times_crab,instdata_crab[2][string(i)],"-p",c=cols[i],label="E$(i) instantaneous")
        end
        for (i,mk) in zip(1:speccount_crab,("x","+"))
            plot(times_crab,tevo_data_crab[2][i][1:end-1],c="k",marker=mk,label="transported manifold state $(i)")
        end
        legend()
        xlabel("Time")
        ylabel("Energy")
        title("Energy vs time for $(method_crab) gradient $(lx)x$(ly) N=$(n) a $(a_start_crab)→$(a_end_crab)\n|B|≤$(bmax_crab), fidelity = $(round(fidelity_rk4_crab,digits=6))")
        tight_layout()
    end

    # plot the best gradient against the smooth guess, and the spacing each drives
    if false || if_all
        fig, axs = subplots(2,1,sharex=true)
        axs[1].plot(times_b_crab,b_guess_crab,"--",c="gray",label="smooth guess")
        axs[1].plot(times_b_crab,b_crab,c="b",label=method_crab)
        for bl in (-bmax_crab,bmax_crab)
            axs[1].axhline(bl,ls=":",c="r")
        end
        axs[1].axhline(bhold_crab,ls="--",c="k",lw=0.8,label="holding gradient w²(a_end-a0)")
        axs[1].set_ylabel("Magnetic gradient B(t)")
        axs[1].set_title("$(method_crab)-optimized gradient (4 freqs/SI, NM, cap $(bmax_crab))\nand the spacing it drives, w=$(omega_crab)")
        axs[1].legend()
        axs[2].plot(times_b_crab,a_guess_crab,"--",c="gray",label="smooth guess, F = $(round(fidelity_guess_crab,digits=4)), R = $(round(ringing_guess_crab,sigdigits=2))")
        axs[2].plot(times_b_crab,a_crab,c="k",label="$(method_crab), F = $(round(fidelity_crab,digits=4)), R = $(round(ringing_crab,sigdigits=2))")
        axs[2].set_xlabel("Time")
        axs[2].set_ylabel("Spacing a(t)")
        axs[2].legend()
        tight_layout()
    end

end=#


### ZVD gradient pulse on the magspac transport: instantaneous vs transported energies
# The ZVD staircase B = Bf [1/4, 3/4, 1] (Bf = w^2 (a_end - a_start)) as a plain pulse, no
# optimization: dd groundstate manifold at a_start -> the one at a_end over T = 2pi/w, the transport
# of the magspacGrad / CRAB blocks above (dd intstren 10, a 0.5 -> 2.0, w = 10), run with run_timeevo
# through the "magnetic_gradient" scaling. Set up for the larger lattice (6x3 N=3, dim 816, 2026-10-02):
# setup_magspac_ramp's 1e-8 degeneracy assert fails there (splitting 1e-3 at a = 0.5, 0.051 at a = 2.0,
# next level 0.335 / 0.163 above), so the endpoint manifolds are diagonalized densely here and only
# required to sit below a gap. Fast-evolution reference at 6x3 N=3 (tolerance relaxed): F 0.9464,
#= quench 0.9313; at 4x4 N=2 ZVD gives 0.9111.
if true

    if_all::Bool = true

    include("magspac-ramp-control-functions.jl")

    # lattice, transport and the half-step grid: a power-of-two count of half-steps over T (so
    # T/(dt/2) is exact and both ZVD switches sit on samples), doubled until dt is below the RK4
    # limit at a_start; the ZVD a(t) rises monotonically, so a_start bounds the whole run
    if false || if_all
        lx,ly,n = 6,3,3
        speccount_zvdE = 2
        nevinst_zvdE = 4         # instantaneous levels to plot
        intstren_zvdE = 10.0
        a_start_zvdE, a_end_zvdE = 0.5, 2.0
        omega_zvdE = 10.0
        ramptime_zvdE = 2pi/omega_zvdE
        # dense diagonalization (endpoints and the instantaneous spectrum) while it fits in memory,
        # Lanczos above: 8x4 N=4 is 35960 states, ~20 GB per dense complex matrix
        if_exact_zvdE = binomial(lx*ly,n) <= 5000

        pdict_zvdE = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("scaling_type","dd"),("trap_frequency",omega_zvdE),("magnetic_spacing",a_start_zvdE),("interaction_strength",intstren_zvdE),("filling",0.5),("nev",speccount_zvdE),("if_reading",true),("if_find_data",false),("if_save_data",false)])
        _,_,_,_,_,lattice_params_zvdE,hamilt_params_zvdE = run_normal_ed(pdict_zvdE; output_level=0)

        dt_crit_zvdE = get_critical_dt(ramptime_zvdE,lattice_params_zvdE,hamilt_params_zvdE)
        nhalf_zvdE = 256 * 2^max(0,ceil(Int,log2(2*ramptime_zvdE/(256*dt_crit_zvdE))))
        dt_zvdE = 2*(ramptime_zvdE/nhalf_zvdE)
        println("dt_crit $(dt_crit_zvdE) -> $(nhalf_zvdE) half-steps, dt $(dt_zvdE)")
    end

    # endpoint manifolds: the speccount lowest states, which must be split by less than the gap above
    # them. Dense where possible, since Lanczos can return only part of a degenerate pair; with
    # Lanczos a missed partner shows up as a failed gap check
    if false || if_all
        manifolds_zvdE = map((a_start_zvdE,a_end_zvdE)) do a
            pd = copy(pdict_zvdE)
            pd["magnetic_spacing"] = a
            if if_exact_zvdE
                _,_,_,_,_,_,hp = run_normal_ed(pd; output_level=0)
                eig = eigen(Hermitian(Matrix(hp["H"])))
                nrgs = eig.values
                vecs = [eig.vectors[:,i] for i in 1:speccount_zvdE]
            else
                pd["nev"] = nevinst_zvdE
                vecs,nrgs,_,_,_,_,_ = run_normal_ed(pd; output_level=0)
                nrgs = real.(nrgs)
            end
            println("a = $(a): lowest levels relative to E0 $(round.(nrgs[1:nevinst_zvdE] .- nrgs[1],digits=4))")
            @assert nrgs[speccount_zvdE] - nrgs[1] < nrgs[speccount_zvdE+1] - nrgs[speccount_zvdE] "the lowest $(speccount_zvdE) levels at a = $(a) are not separated from the rest by a gap"
            [Vector{ComplexF64}(vecs[i]) for i in 1:speccount_zvdE]
        end
        starting_states_zvdE, target_states_zvdE = manifolds_zvdE
        println("Quench (overlap of the endpoint manifolds): $(real(groundstate_manifold_fidelity(starting_states_zvdE,target_states_zvdE)))")
    end

    # the ZVD gradient on the run grid and the spacing it drives
    if false || if_all
        b_zvdE = zvd_gradient_pulse(a_start_zvdE,a_end_zvdE,omega_zvdE,ramptime_zvdE,dt_zvdE)
        times_b_zvdE = halfstep_times(ramptime_zvdE,dt_zvdE)
        a_zvdE = a_start_zvdE .+ magnetic_gradient_response(b_zvdE,dt_zvdE/2,omega_zvdE)[1]
        @assert minimum(a_zvdE) >= a_start_zvdE - 1e-12 "a(t) dips below a_start, so dt_crit from the t = 0 Hamiltonian does not bound the run"
    end

    # time evolution with the instantaneous spectrum along the way
    if false || if_all
        # work on a copy: timeham writes the ramped values back into the dict
        hamilt_params_tevo_zvdE = copy(hamilt_params_zvdE)
        hamilt_params_tevo_zvdE["scaling_type"] = "magnetic_gradient"
        time_running_args_zvdE = (nev=nevinst_zvdE,output_level=0,if_instant_gs=true,if_instant_exact=if_exact_zvdE,if_reading=true,if_save_data=false,dataloc="tevo-daily-things-data/")
        tevo_params_zvdE = Dict([ ("magnetic_gradient_time",(pulse_ramp,ramptime_zvdE,b_zvdE)),("trap_frequency",omega_zvdE),("tmax",ramptime_zvdE),("dt",dt_zvdE) ])
        tevo_data_zvdE,tevo_dict_zvdE,instdata_zvdE,_ = run_timeevo(starting_states_zvdE,tevo_params_zvdE,lattice_params_zvdE,hamilt_params_tevo_zvdE; time_running_args_zvdE...)

        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        final_manifold_zvdE = [Vector{ComplexF64}(tevo_data_zvdE[1][i][:,end-1]) for i in 1:speccount_zvdE]
        fidelity_zvdE = real(groundstate_manifold_fidelity(final_manifold_zvdE,target_states_zvdE))
        println("Fidelity with target manifold after the ZVD pulse: $(fidelity_zvdE)")
    end

    # spacing and gradient on top, instantaneous levels vs the transported energies of the two
    # evolved manifold states below
    if false || if_all
        fig, axs = subplots(2,1,sharex=true,figsize=(8,8),gridspec_kw=Dict("height_ratios"=>[1,2]))
        axs[1].plot(times_b_zvdE,a_zvdE,c="k")
        axs[1].set_ylabel("Spacing a(t)")
        axb = axs[1].twinx()
        axb.plot(times_b_zvdE,b_zvdE,c="b",lw=0.8)
        axb.set_ylabel("Gradient B(t)",color="b")
        axs[1].set_title("ZVD gradient pulse $(lx)x$(ly) N=$(n), a $(a_start_zvdE)→$(a_end_zvdE), w=$(omega_zvdE)\nfidelity = $(round(fidelity_zvdE,digits=6))")

        times_zvdE = range(0.0,ramptime_zvdE,length=length(instdata_zvdE[2]["1"]))
        cols = ["b","g","r","m","c","y"]
        for i in 1:nevinst_zvdE
            axs[2].plot(times_zvdE,instdata_zvdE[2][string(i)],"-",c=cols[i],label="E$(i) instantaneous")
        end
        for (i,mk) in zip(1:speccount_zvdE,("x","+"))
            axs[2].plot(times_zvdE,tevo_data_zvdE[2][i][1:end-1],c="k",ls="none",marker=mk,ms=5,markevery=max(1,length(times_zvdE)÷40),label="transported manifold state $(i)")
        end
        axs[2].legend()
        axs[2].set_xlabel("Time")
        axs[2].set_ylabel("Energy")
        tight_layout()
    end

end=#



### ZVD gradient pulse on the magspac transport: instantaneous vs transported energies, 8x4 N=4
# Same as the block above, run at 8x4 N=4 (35960 states): quench overlap 0.848 (vs 0.896 at 4x4,
# 0.931 at 6x3), the first size where it drops. Too big for dense diagonalization (~20 GB per
# matrix), so the endpoints and the instantaneous spectrum use Lanczos; at a = 0.5 Lanczos returns
# the exactly degenerate pair (next level 0.366 above), at a = 2.0 a pair split by 0.004 (next 0.286).
# 512 half-steps (dt 0.0025). H is read and dressed (dressed hopping cached), ~50-70 ms per build;
# the per-step Lanczos for the instantaneous levels dominates, expect ~30-60 min. Prints progress.
if true

    if_all::Bool = true

    include("magspac-ramp-control-functions.jl")

    # lattice, transport and the half-step grid: a power-of-two count of half-steps over T (so
    # T/(dt/2) is exact and both ZVD switches sit on samples), doubled until dt is below the RK4
    # limit at a_start; the ZVD a(t) rises monotonically, so a_start bounds the whole run
    if false || if_all
        lx,ly,n = 8,4,4
        speccount_zvd84 = 2
        nevinst_zvd84 = 4         # instantaneous levels to plot
        intstren_zvd84 = 10.0
        a_start_zvd84, a_end_zvd84 = 0.5, 2.0
        omega_zvd84 = 10.0
        ramptime_zvd84 = 2pi/omega_zvd84
        # dense diagonalization (endpoints and the instantaneous spectrum) while it fits in memory,
        # Lanczos above: 8x4 N=4 is 35960 states, ~20 GB per dense complex matrix
        if_exact_zvd84 = binomial(lx*ly,n) <= 5000

        pdict_zvd84 = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("scaling_type","dd"),("trap_frequency",omega_zvd84),("magnetic_spacing",a_start_zvd84),("interaction_strength",intstren_zvd84),("filling",0.5),("nev",speccount_zvd84),("if_reading",true),("if_find_data",false),("if_save_data",false)])
        _,_,_,_,_,lattice_params_zvd84,hamilt_params_zvd84 = run_normal_ed(pdict_zvd84; output_level=0)

        dt_crit_zvd84 = get_critical_dt(ramptime_zvd84,lattice_params_zvd84,hamilt_params_zvd84)
        nhalf_zvd84 = 256 * 2^max(0,ceil(Int,log2(2*ramptime_zvd84/(256*dt_crit_zvd84))))
        dt_zvd84 = 2*(ramptime_zvd84/nhalf_zvd84)
        println("dt_crit $(dt_crit_zvd84) -> $(nhalf_zvd84) half-steps, dt $(dt_zvd84)")
    end

    # endpoint manifolds: the speccount lowest states, which must be split by less than the gap above
    # them. Dense where possible, since Lanczos can return only part of a degenerate pair; with
    # Lanczos a missed partner shows up as a failed gap check
    if false || if_all
        manifolds_zvd84 = map((a_start_zvd84,a_end_zvd84)) do a
            pd = copy(pdict_zvd84)
            pd["magnetic_spacing"] = a
            if if_exact_zvd84
                _,_,_,_,_,_,hp = run_normal_ed(pd; output_level=0)
                eig = eigen(Hermitian(Matrix(hp["H"])))
                nrgs = eig.values
                vecs = [eig.vectors[:,i] for i in 1:speccount_zvd84]
            else
                pd["nev"] = nevinst_zvd84
                vecs,nrgs,_,_,_,_,_ = run_normal_ed(pd; output_level=0)
                nrgs = real.(nrgs)
            end
            println("a = $(a): lowest levels relative to E0 $(round.(nrgs[1:nevinst_zvd84] .- nrgs[1],digits=4))")
            @assert nrgs[speccount_zvd84] - nrgs[1] < nrgs[speccount_zvd84+1] - nrgs[speccount_zvd84] "the lowest $(speccount_zvd84) levels at a = $(a) are not separated from the rest by a gap"
            [Vector{ComplexF64}(vecs[i]) for i in 1:speccount_zvd84]
        end
        starting_states_zvd84, target_states_zvd84 = manifolds_zvd84
        println("Quench (overlap of the endpoint manifolds): $(real(groundstate_manifold_fidelity(starting_states_zvd84,target_states_zvd84)))")
    end

    # the ZVD gradient on the run grid and the spacing it drives
    if false || if_all
        b_zvd84 = zvd_gradient_pulse(a_start_zvd84,a_end_zvd84,omega_zvd84,ramptime_zvd84,dt_zvd84)
        times_b_zvd84 = halfstep_times(ramptime_zvd84,dt_zvd84)
        a_zvd84 = a_start_zvd84 .+ magnetic_gradient_response(b_zvd84,dt_zvd84/2,omega_zvd84)[1]
        @assert minimum(a_zvd84) >= a_start_zvd84 - 1e-12 "a(t) dips below a_start, so dt_crit from the t = 0 Hamiltonian does not bound the run"
    end

    # time evolution with the instantaneous spectrum along the way
    if false || if_all
        # work on a copy: timeham writes the ramped values back into the dict
        hamilt_params_tevo_zvd84 = copy(hamilt_params_zvd84)
        hamilt_params_tevo_zvd84["scaling_type"] = "magnetic_gradient"
        time_running_args_zvd84 = (nev=nevinst_zvd84,output_level=1,if_instant_gs=true,if_instant_exact=if_exact_zvd84,if_reading=true,if_save_data=false,dataloc="tevo-daily-things-data/")
        tevo_params_zvd84 = Dict([ ("magnetic_gradient_time",(pulse_ramp,ramptime_zvd84,b_zvd84)),("trap_frequency",omega_zvd84),("tmax",ramptime_zvd84),("dt",dt_zvd84) ])
        tevo_data_zvd84,tevo_dict_zvd84,instdata_zvd84,_ = run_timeevo(starting_states_zvd84,tevo_params_zvd84,lattice_params_zvd84,hamilt_params_tevo_zvd84; time_running_args_zvd84...)

        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        final_manifold_zvd84 = [Vector{ComplexF64}(tevo_data_zvd84[1][i][:,end-1]) for i in 1:speccount_zvd84]
        fidelity_zvd84 = real(groundstate_manifold_fidelity(final_manifold_zvd84,target_states_zvd84))
        println("Fidelity with target manifold after the ZVD pulse: $(fidelity_zvd84)")
    end

    # spacing and gradient on top, instantaneous levels vs the transported energies of the two
    # evolved manifold states below
    if false || if_all
        fig, axs = subplots(2,1,sharex=true,figsize=(8,8),gridspec_kw=Dict("height_ratios"=>[1,2]))
        axs[1].plot(times_b_zvd84,a_zvd84,c="k")
        axs[1].set_ylabel("Spacing a(t)")
        axb = axs[1].twinx()
        axb.plot(times_b_zvd84,b_zvd84,c="b",lw=0.8)
        axb.set_ylabel("Gradient B(t)",color="b")
        axs[1].set_title("ZVD gradient pulse $(lx)x$(ly) N=$(n), a $(a_start_zvd84)→$(a_end_zvd84), w=$(omega_zvd84)\nfidelity = $(round(fidelity_zvd84,digits=6))")

        times_zvd84 = range(0.0,ramptime_zvd84,length=length(instdata_zvd84[2]["1"]))
        cols = ["b","g","r","m","c","y"]
        for i in 1:nevinst_zvd84
            axs[2].plot(times_zvd84,instdata_zvd84[2][string(i)],"-",c=cols[i],label="E$(i) instantaneous")
        end
        for (i,mk) in zip(1:speccount_zvd84,("x","+"))
            axs[2].plot(times_zvd84,tevo_data_zvd84[2][i][1:end-1],c="k",ls="none",marker=mk,ms=5,markevery=max(1,length(times_zvd84)÷40),label="transported manifold state $(i)")
        end
        axs[2].legend()
        axs[2].set_xlabel("Time")
        axs[2].set_ylabel("Energy")
        tight_layout()
    end

end





























"fin"


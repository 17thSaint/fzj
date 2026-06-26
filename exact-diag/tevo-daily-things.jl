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

if_all::Bool = true

# model parameters
if false || if_all
    lx,ly,n = 4,4,2
    
    intstren = 0.0
end

# define starting state
if false || if_all
    starting_config = [(1,1),(1,2)]

    psi_starting = position_state

        





















































"fin"
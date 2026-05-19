#####################################################
#=

This file contains tests of the time evolution code

Depends on:
    execute-ed.jl

=#
######################################################

include("execute-ed.jl")
include("time-evolution.jl")
include("control-functions.jl")
include("../other-funcs/basic-2d-plottings.jl")
include("plottings.jl")




#= testing the proper overlap with the manifold as a fidelity measure
if false
    lx,ly,n = 4,4,2
    intstren = 10.0
    num_samples = 100

    start_xi = 50.0
    pdict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("scaling_type","exp"),("corr_length",start_xi),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])

    all_states_gs1 = []
    all_states_gs2 = []

    for i in 1:num_samples
        states1,nrgs1,rhos1,filepath1,if_found1,lattice_params1,hamilt_params1 = run_normal_ed(pdict; output_level=0)
        push!(all_states_gs1,states1[1])
        push!(all_states_gs2,states1[2])
    end

    all_fidelities = zeros(Float64,num_samples,num_samples)
    for i1 in 1:num_samples
        for i2 in i1+1:num_samples
            fid_mat = zeros(ComplexF64,2,2)
            fid_mat[1,1] = adjoint(all_states_gs1[i1]) * all_states_gs1[i2]
            fid_mat[1,2] = adjoint(all_states_gs1[i1]) * all_states_gs2[i2]
            fid_mat[2,1] = adjoint(all_states_gs2[i1]) * all_states_gs1[i2]
            fid_mat[2,2] = adjoint(all_states_gs2[i1]) * all_states_gs2[i2]
            all_fidelities[i1,i2] = 0.5 * tr(adjoint(fid_mat)*fid_mat)
        end
    end

    fig = figure()
    imshow(all_fidelities,origin="lower",vmin=0,vmax=1)
    colorbar()
    xlabel("Sample index")
    ylabel("Sample index")
    title("Fidelity between different samples of the same state")
end=#

#= 4x4 N=2 ULR check
if false
    lx,ly,n = 4,4,2
    
    intstrens = range(0.0,10.0,length=20)

    for intstren in intstrens
        params_dict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])

        states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)

        if intstren == 0.0
            scatter(intstren,nrgs[1] - nrgs[1],c="b",label="E0")
            scatter(intstren,nrgs[2] - nrgs[1],c="g",label="E1")
            scatter(intstren,nrgs[3] - nrgs[1],c="k",label="E2")
            for i in 4:length(nrgs)
                scatter(intstren,nrgs[i] - nrgs[1],c="k")
            end
        else
            scatter(intstren,nrgs[1] - nrgs[1],c="b")
            scatter(intstren,nrgs[2] - nrgs[1],c="g")
            for i in 3:length(nrgs)
                scatter(intstren,nrgs[i] - nrgs[1],c="k")
            end
        end

        xlabel("Interaction strength")
        ylabel("Energy - E0")
        title("$(lx)x$(ly) N=$(n) Spectrum")
        legend()
    end
end=#

#= 4x4 N=2 twisting check
if false
    lx,ly,n = 4,4,2
    
    twist_angles = range(0.0,1.0,length=11)

    gaps02 = zeros(Float64,length(twist_angles),length(twist_angles))
    gaps12 = zeros(Float64,length(twist_angles),length(twist_angles))
    for (i1,tw1) in enumerate(twist_angles)
        for (i2,tw2) in enumerate(twist_angles)
            params_dict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("tw1",tw1),("tw2",tw2),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",0.0),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])

            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)

            gaps02[i1,i2] = nrgs[3] - nrgs[1]
            gaps12[i1,i2] = nrgs[3] - nrgs[2]
        end
    end

    fig = figure()
    imshow(gaps02,origin="lower",extent=(0,1,0,1),vmin=0)
    colorbar()
    xlabel(L"\theta_x/ 2\pi")
    ylabel(L"\theta_y/ 2\pi")
    title("$(lx)x$(ly) N=$(n) Gap E2-E0")

    fig = figure()
    imshow(gaps12,origin="lower",extent=(0,1,0,1),vmin=0)
    colorbar()
    xlabel(L"\theta_x/ 2\pi")
    ylabel(L"\theta_y/ 2\pi")
    title("$(lx)x$(ly) N=$(n) Gap E2-E1")

end=#

#= 4x4 N=2 dt benchmark
if false

    lx,ly,n = 4,4,2
    
    intstren = 0.0

    tmax_global = 10.0
    #dt_global = 0.001
    dataloc = get_folder_location("cluster-data/exact-diag/time-evo/dt-benchmark")

    start_tx = 1e-8

    params_dict_starting = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("tx",start_tx),("ty",1.0),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",5),("if_find_data",false),("if_save_data",false)])

    states_i,nrgs_i,rhos_i,filepath_i,if_found_i,lattice_params_i,hamilt_params_i = run_normal_ed(params_dict_starting; output_level=0)
    println("Found starting state")

    psi0_1 = states_i[1]
    psi0_2 = states_i[2]

    end_tx = 1.0

    params_dict_ending = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("tx",end_tx),("ty",1.0),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",5),("if_find_data",false),("if_save_data",false)])

    states_f,nrgs_f,rhos_f,filepath_f,if_found_f,lattice_params_f,hamilt_params_f = run_normal_ed(params_dict_ending; output_level=0)
    println("Found ending state")

    psif_1 = states_f[1]
    psif_2 = states_f[2]


    speccount = 2
    time_running_args = (nev=speccount,output_level=1,if_instant_gs=false,if_save_data=false,dataloc=dataloc,)

    tmax_global = 0.5
    ramptime = 0.1
    dts = 10 .^ range(-3,-1,length=21)
    for dt_global in dts
        tevo_params = Dict([ ("tx",(linear_ramp,start_tx,end_tx,ramptime)),("dt",dt_global),("tmax",tmax_global) ])
        tevo_gs,tevo_dict,intspec,saving_args = run_timeevo([psi0_1,psi0_2],tevo_params,lattice_params_i,hamilt_params_i; time_running_args...)

        final_gs_manifold = [tevo_gs[1][:,end-1],tevo_gs[2][:,end-1]]
        final_fidelity = groundstate_manifold_fidelity(states_f[1:speccount],final_gs_manifold)
        println("Final fidelity for dt=$(dt_global): $(final_fidelity)")
        scatter(dt_global,final_fidelity,c="b")
    end
    xscale("log")
    xlabel("Time step dt")
    ylabel("Fidelity with target manifold")
    title("Fidelity vs dt $(lx)x$(ly) N=$(n) U=$(intstren) ramp tx ramptime=$(ramptime)")

end=#





































"fin"
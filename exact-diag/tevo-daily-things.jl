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

if_all::Bool = true
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
if true || if_all
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

end#




    

        





















































"fin"
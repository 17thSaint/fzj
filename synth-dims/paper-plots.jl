#####################################################
#=

This file contains plotting functions for paper

Depends on:
    review-practice-codes/long-range-ttn.jl
    review-practice-codes/plottings.jl

=#
######################################################

include("../other-funcs/include-other-files.jl")

include_other_files(["synth-dims/long-range-ttn.jl","synth-dims/hatsugai-mbcn.jl","other-funcs/basic-2d-observables.jl"])
include_other_files(["review-practice-codes/plottings.jl","other-funcs/basic-2d-plottings.jl"])
include_other_files(["exact-diag/execute-ed.jl","exact-diag/observables.jl","exact-diag/plottings.jl"])
using JLD2,LaTeXStrings,LsqFit

########### TO-DO ############
#=

       Plots to be made
1. Hatsugai Chern number counting (or use old figure)
3. Fourier transform of density thermodynamic scaling to show flatness


       Data to be run
1. E2 and E3 for 16x8 N=8 ULR=(0.0 and 300.0)

=#
##############################



# plot of gap staying open for rho_1D = 0.5 8x4 N=4
function plot_opengap_8x4_ed()
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",1.0),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
    filter!(x -> !occursin("twist_angle",x), all_files)
    #display(all_files)

    for f in all_files
        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        all_nrgs = d["nrg"]
        
        params = get_params_dict_from_filename(f)
        ulr = params["interaction_strength"]

        scatter(ulr,0.0,c="r")
        scatter(ulr,all_nrgs[2] - all_nrgs[1],c="r")

        for i in 3:length(all_nrgs)
            scatter(ulr,all_nrgs[i] - all_nrgs[1],c="k")
        end
    end
    xlabel("Interaction Strength")
    ylabel("E - E0")
    title("Energy Spectrum $(lx)x$(ly) N=$(n)")
    ylim(-0.02,0.6)

end
#plot_opengap_8x4_ed()

# plot of gap staying open for rho_1D = 0.5 8x4 N=4
function plot_closedgap_8x4_ed()
    lx,ly,n = 4,8,4
    layers = Int(log(2,lx*ly))

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",1.0),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
    filter!(x -> !occursin("twist_angle",x), all_files)
    #display(all_files)

    for f in all_files
        params = get_params_dict_from_filename(f)
        ulr = params["interaction_strength"]

        if ulr > 10.0
            continue
        end

        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        all_nrgs = d["nrg"]

        scatter(ulr,0.0,c="r")
        scatter(ulr,all_nrgs[2] - all_nrgs[1],c="r")

        for i in 3:length(all_nrgs)
            scatter(ulr,all_nrgs[i] - all_nrgs[1],c="k")
        end
    end
    xlabel("Interaction Strength")
    ylabel("E - E0")
    title("Energy Spectrum $(lx)x$(ly) N=$(n)")
    #ylim(-0.02,0.6)

end
#plot_closedgap_8x4_ed()

# finite size scaling of the finite size gap
# waiting for 16x8 to finish running GS2
function plot_finitesize_gapscaling(ulr::Float64=0.0)
    fig = figure()

    # ED section
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("hopping_anisotropy",1.0),("interaction_strength",ulr),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2", output_level=0)
    filter!(x -> !occursin("twist_angle",x), all_files)
    display(all_files)

    used_files = []
    lxs = []
    for f in all_files
        params = get_params_dict_from_filename(f)

        if (params["N"] / params["Lx"] != 0.5) || (params["Lx"] != 2*params["Ly"])
            continue
        end

        append!(used_files,[f])

        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        all_nrgs = d["nrg"]

        scatter(params["Lx"],all_nrgs[2] - all_nrgs[1],c="b")
        append!(lxs,[params["Lx"]])
    end
    xlabel("Lx")
    ylabel("E1 - E0")
    title("Finite Gap Scaling for ULR=$ulr")
    yscale("log")


    #= TTN section
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("onsite_strength",ulr),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)

    for f in all_files_ttn
        params = get_params_dict_from_filename(f)
        Lx,Ly = get_lattice_dims_from_layers(params["layers"])


        if (params["particles"] / Lx != 0.5) || (Lx != 2*Ly) || Lx <= 10
            continue
        end

        append!(used_files,[f])

        d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

        e1 = m["observer"].nrg[end]
        e2 = m["observer_1"].nrg[end]
        scatter(Lx,e2 - e1,c="r")
        append!(lxs,[Lx])
    end
    =#

    xlim(0.0,1.1*(maximum(lxs)))

    return used_files
end
#kept_files = plot_finitesize_gapscaling(0.0)

# finite size scaling of the finite size gap with density perturbation at (1,1) strength 1e-4
function plot_finitesize_gapscaling_pinned(ulr::Float64=300.0)
    fig = figure()

    # ED section
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/pinned-scaling")
    pdict = Dict([("hopping_anisotropy",1.0),("if_pinning",true),("interaction_strength",ulr),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
    filter!(x -> !occursin("twist_angle",x), all_files)
    display(all_files)

    used_files = []
    lxs = []
    for f in all_files
        params = get_params_dict_from_filename(f)

        if (params["N"] / params["Lx"] != 0.5) || (params["Lx"] != 2*params["Ly"])
            continue
        end

        append!(used_files,[f])

        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        all_nrgs = d["nrg"]

        scatter(params["Lx"],all_nrgs[2] - all_nrgs[1],c="b")
        append!(lxs,[params["Lx"]])
    end
    xlabel("Lx")
    ylabel("E1 - E0")
    title("Finite Gap Scaling of Pinned State for ULR=$ulr")
    yscale("log")


    #= TTN section
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge/pinned-scaling")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("if_pinning",true),("onsite_strength",ulr),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)

    for f in all_files_ttn
        params = get_params_dict_from_filename(f)
        Lx,Ly = get_lattice_dims_from_layers(params["layers"])


        if (params["particles"] / Lx != 0.5) || (Lx != 2*Ly) || Lx <= 10
            continue
        end

        append!(used_files,[f])

        d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

        e1 = m["observer"].nrg[end]
        e2 = m["observer_1"].nrg[end]
        scatter(Lx,e2 - e1,c="r")
        append!(lxs,[Lx])
    end=#

    xlim(0.0,1.1*(maximum(lxs)))

    return used_files
end

# finite size scaling of the topological gap
# still need 16x8 to get E3
function plot_topo_gapscaling(ulr::Float64=0.0)

    # ED section
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("hopping_anisotropy",1.0),("interaction_strength",ulr),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
    filter!(x -> !occursin("twist_angle",x), all_files)
    #display(all_files)

    used_files = []
    allvals = []
    lxs = []
    for f in all_files
        params = get_params_dict_from_filename(f)

        if (params["N"] / params["Lx"] != 0.5) || (params["Lx"] != 2*params["Ly"])
            continue
        end

        append!(used_files,[f])

        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        all_nrgs = d["nrg"]

        scatter(params["Lx"],all_nrgs[3] - all_nrgs[1],c="b")
        
        append!(allvals,[all_nrgs[3] - all_nrgs[1]])
        append!(lxs,[params["Lx"]])
    end
    xlabel("Lx")
    ylabel("E2 - E0")
    title("Topological Gap Scaling")
    ylim(0.0,1.1*maximum(allvals))
    xlim(0.0,1.1*maximum(lxs))

    #= TTN section
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("onsite_strength",ulr),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)

    for f in all_files_ttn
        params = get_params_dict_from_filename(f)
        Lx,Ly = get_lattice_dims_from_layers(params["layers"])


        if (params["particles"] / Lx != 0.5) || (Lx != 2*Ly) || Lx <= 10
            continue
        end

        append!(used_files,[f])

        d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

        e1 = m["observer"].nrg[end]
        e2 = m["observer_1"].nrg[end]
        scatter(Lx,e2 - e1,c="r")
    end=#

    return used_files
end
#plot_topo_gapscaling(0.0)

# perimeter scaling of entanglement entropy for 16x8 N=8
function plot_ee_scaling()
    lx,ly,n = 16,8,8
    layers = Int(log(2,lx*ly))

    
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)

    perims = [8*3,4*4,4*3,2*4,2*3]
    all_ees = []
    for f in all_files_ttn

        params = get_params_dict_from_filename(f)
        d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

        ees = zeros(Float64,layers-2)
        entspecs = real.(m["entanglement_spectrum"])

        for k in 2:layers-1
            entspec = filter(x -> x != 0.0, entspecs[k,:])
            #display(entspec)

            ee = entanglement_entropy(entspec)
            ees[k-1] = ee
        end
        scatter(perims,ees,label="$(params["onsite_strength"])")
        append!(all_ees,[ees])
    end

    xlabel("Perimeter")
    ylabel("Entanglement Entropy")
    #title("Entanglement Entropy vs Perimeter $(lx)x$(ly) N=$(n)")
    #legend()
    ylim(-1.0,1.1*maximum(Iterators.flatten(all_ees)))
    xlim(0.0,1.1*maximum(perims))

    linmodel(x,p) = p[1] .* x .+ p[2]
    linfit = curve_fit(linmodel,perims,all_ees[1],[1.0,1.0])
    yintercept = linfit.param[2]
    title("Entanglement Entropy vs Perimeter $(lx)x$(ly) N=$(n)")
    xs = [0.0,1.2*maximum(perims)]
    plot(xs,linmodel(xs,linfit.param),label="Fit: "*L"$\gamma = $"*"$(round(-yintercept, digits=3))",c="b")
    legend()

end
#plot_ee_scaling()



































"fin"
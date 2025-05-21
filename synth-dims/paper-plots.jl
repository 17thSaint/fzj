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

    #xlim(0.0,1.1*(maximum(lxs)))

    return used_files
end
#kept_files = plot_finitesize_gapscaling(0.0)

# finite size scaling of the finite size gap with density perturbation at (1,1) strength 1e-4
function plot_finitesize_gapscaling_pinned(ulr::Float64=300.0)

    # ED section
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/pinned-scaling")
    pdict = Dict([("hopping_anisotropy",1.0),("if_pinning",true),("interaction_strength",ulr),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
    filter!(x -> !occursin("twist_angle",x), all_files)
    display(all_files)

    used_files = []
    lxs = Dict("0.01"=>[], "0.0001"=>[], "0.001"=>[])
    gaps = Dict("0.01"=>[], "0.0001"=>[], "0.001"=>[])
    for f in all_files
        params = get_params_dict_from_filename(f)

        if (params["N"] / params["Lx"] != 0.5) || (params["Lx"] != 2*params["Ly"])
            continue
        end

        append!(used_files,[f])

        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        all_nrgs = d["nrg"]

        if haskey(m,"pinning_strength")
            if haskey(gaps,string(m["pinning_strength"]))
                append!(gaps[string(m["pinning_strength"])],[all_nrgs[2] - all_nrgs[1]])
                append!(lxs[string(m["pinning_strength"])],[params["Lx"]])
            else
                gaps[string(m["pinning_strength"])] = [all_nrgs[2] - all_nrgs[1]]
                lxs[string(m["pinning_strength"])] = [params["Lx"]]
            end
        else
            append!(gaps["0.0001"],[all_nrgs[2] - all_nrgs[1]])
            append!(lxs["0.0001"],[params["Lx"]])
        end
        #append!(lxs[1],[params["Lx"]])
        #append!(gaps[1],[all_nrgs[2] - all_nrgs[1]])

        #scatter(params["Lx"],all_nrgs[2] - all_nrgs[1],c=col)

    end
    fig = figure()
    for (k,v) in gaps
        length(v) > 0 && scatter(lxs[k],v,label="$(parse(Float64,k)/10)")
    end
    xlabel("Lx")
    ylabel("E1 - E0")
    title("Finite Gap Scaling of Pinned State for ULR=$ulr")
    yscale("log")
    legend()

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

    #xlim(0.0,1.1*(maximum(lxs)))

    return used_files
end

# plot finite splitting scaling pinned and unpinned
function plot_finitesplitting_scaling(ulr::Float64=0.0)
    fig = figure()

    pinning_strength = ulr == 0.0 ? 0.1 : 0.0001

    # ED section: pinned
    dataloc_pin = get_folder_location("cluster-data/exact-diag/torus/new-gauge/pinned-scaling")
    pdict_pin = Dict([("hopping_anisotropy",1.0),("if_pinning",true),("interaction_strength",ulr),("if_periodic_x",true),("if_periodic_y",true)])
    all_files_pin = find_data_file(pdict_pin,"ed",dataloc_pin; file_type="jld2")
    display(all_files_pin)

    used_files = []
    lxs = Dict("0.1"=>[], "0.0001"=>[], "0.001"=>[],"0.0"=>[])
    gaps = Dict("0.1"=>[], "0.0001"=>[], "0.001"=>[],"0.0"=>[])
    for f in all_files_pin
        params = get_params_dict_from_filename(f)

        if (params["N"] / params["Lx"] != 0.5) || (params["Lx"] != 2*params["Ly"])
            continue
        end

        append!(used_files,[f])

        d,m = read_data_jld2(joinpath(dataloc_pin,f); output_level=0)
        all_nrgs = d["nrg"]

        if haskey(m,"pinning_strength") && m["pinning_strength"] == pinning_strength
            append!(gaps[string(pinning_strength)],[all_nrgs[2] - all_nrgs[1]])
            append!(lxs[string(pinning_strength)],[params["Lx"]])
        end
    end

    # ED section: unpinned
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("hopping_anisotropy",1.0),("interaction_strength",ulr),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
    display(all_files)

    for f in all_files
        params = get_params_dict_from_filename(f)

        if (params["N"] / params["Lx"] != 0.5) || (params["Lx"] != 2*params["Ly"])
            continue
        end

        append!(used_files,[f])

        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        all_nrgs = d["nrg"]

        append!(gaps["0.0"],[all_nrgs[2] - all_nrgs[1]])
        append!(lxs["0.0"],[params["Lx"]])
    end

    # TTN section: pinned
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge/pinned-scaling")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("onsite_strength",ulr),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)

    for f in all_files_ttn
        params = get_params_dict_from_filename(f)
        Lx,Ly = "Lx" in keys(params) ? (params["Lx"],params["Ly"]) : get_lattice_dims_from_layers(params["layers"])


        if (params["particles"] / Lx != 0.5) || (Lx != 2*Ly) || Lx <= 10
            continue
        end

        d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

        if haskey(m,"observer") && haskey(m,"observer_1")
            append!(used_files,[f])
            append!(gaps[string(pinning_strength)], [abs(m["observer_1"].nrg[end] - m["observer"].nrg[end])])
            append!(lxs[string(pinning_strength)], [Lx])
        end

    end

    # TTN section: unpinned
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("onsite_strength",ulr),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)

    for f in all_files_ttn
        params = get_params_dict_from_filename(f)
        Lx,Ly = "Lx" in keys(params) ? (params["Lx"],params["Ly"]) : get_lattice_dims_from_layers(params["layers"])

        # 16x8 is not converged yet
        Lx == 16 && continue

        if (params["particles"] / Lx != 0.5) || (Lx != 2*Ly) || Lx <= 10
            continue
        end

        d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

        if haskey(m,"observer") && haskey(m,"observer_1")
            append!(used_files,[f])
            append!(gaps["0.0"], [abs(m["observer_1"].nrg[end] - m["observer"].nrg[end])])
            append!(lxs["0.0"], [Lx])
        end

    end
    #



    for (k,v) in gaps
        local_label = k == "0.0" ? "Unpinned" : "$(parse(Float64,k)/10)"
        length(v) > 0 && scatter(lxs[k],v,label=local_label)
    end
    xlabel("Lx")
    ylabel("E1 - E0")
    title("Finite Splitting Scaling for ULR=$ulr")
    yscale("log")
    legend()

    return gaps,lxs
end
#rez_gaps,rez_lxs = plot_finitesplitting_scaling(0.0)

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

# correlations for 16x8 N=8
function plot_correlations(ulr::Float64=0.0)
    lx,ly,n = 16,8,8
    layers = Int(log(2,lx*ly))

    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("onsite_strength",ulr),("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)

    if length(all_files_ttn) > 1
        error("Multiple files found: $(all_files_ttn)")
    end
    
    f = all_files_ttn[1]

    params = get_params_dict_from_filename(f)
    d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

    rho = d["densmat"]
    corrs = physical_correlation(rho,lx,ly; if_plot=true,plot_title="$(lx)x$(ly) N=$n ULR=$ulr")
    #corrs = synthetic_correlation(rho,lx,ly; if_plot=true,plot_title="$(lx)x$(ly) N=$n ULR=$ulr")

    return corrs
end

# fourpt flatness scaling
# get fourpt data for ED 10x5
function plot_fourpt_flatness_scaling()

    # ED section
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict = Dict([("hopping_anisotropy",1.0),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2", output_level=0)
    filter!(x -> !occursin("twist_angle",x), all_files)
    display(all_files)

    ratios = Dict([("8",[]),("10",[]),("16",[])])
    intstrens = Dict([("8",[]),("10",[]),("16",[])])
    for f in all_files
        params = get_params_dict_from_filename(f)
        if params["Lx"] <= 6
            continue
        end

        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        
        !haskey(m,"fourpt_momentum") && continue

        fourpt_vals = m["fourpt_momentum"]
        #plot_four_point(fourpt_vals; plot_title="$(params["Lx"])x$(params["Ly"]) N=$(params["N"]) ULR=$ulr")

        restricted_fourpts = fourpt_vals[1,4:end-2]

        maxval = maximum(restricted_fourpts)
        minval = minimum(restricted_fourpts)
        maxminratio = minval / maxval
        append!(ratios[string(params["Lx"])],[maxminratio])
        append!(intstrens[string(params["Lx"])],[params["interaction_strength"]])
    end

    # TTN section
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("layers",7),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)

    for f in all_files_ttn
        params = get_params_dict_from_filename(f)
        Lx,Ly = get_lattice_dims_from_layers(params["layers"])

        d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

        fourpt_vals = m["fourpt_momentum"]
        #plot_four_point(fourpt_vals; plot_title="$(params["Lx"])x$(params["Ly"]) N=$(params["N"]) ULR=$ulr")

        restricted_fourpts = fourpt_vals[1,4:end-2]

        maxval = maximum(restricted_fourpts)
        minval = minimum(restricted_fourpts)
        maxminratio = minval / maxval
        append!(ratios[string(Lx)], [maxminratio])
        append!(intstrens[string(Lx)], [params["onsite_strength"]])
    end

    for (k,v) in ratios
        if length(v) > 0
            scatter(intstrens[k],v,label="$(k)")
        end
    end
    xlabel("Interaction Strength")
    ylabel("Min/Max Fourpt")
    legend()
    yscale("log")
    
    return
end

# fourpt momentum diagonal max as order parameter for k-DW transition
function plot_fourpt_kdw_transition()
    #=layers = 7

    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("layers",layers),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)

    orderparams_ttn = Float64[]
    intstrens_ttn = Float64[]
    for f in all_files_ttn
        d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

        if !("fourpt_momentum" in keys(m))
            continue
        end

        fourpt_vals = m["fourpt_momentum"]

        orderparam = maximum(diag(fourpt_vals))
        append!(orderparams_ttn,[orderparam])
        append!(intstrens_ttn,[m["onsite_strength"]])
    end
    scatter(intstrens_ttn,orderparams_ttn ./ orderparams_ttn[1],label="TTN N=8")
    xlabel("ULR")
    ylabel("Max Fourpt")
    title("Order Parameter for k-DW transition")=#

    laughlin_values = []
    for n in [3,4,5]
        lx,ly,n = Int(2*n),n,n

        orderparams_ed = Float64[]
        intstrens_ed = Float64[]

        dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",1.0),("if_periodic_x",true),("if_periodic_y",true)])
        all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
        for f in all_files
            d,m = read_data_jld2(joinpath(dataloc,f); output_level=1)

            if n == 5
                !haskey(m,"fourpt_momentum_diag") && continue
                fourpt_vals = m["fourpt_momentum_diag"]
                append!(orderparams_ed,[maximum(fourpt_vals)])
            else
                fourpt_vals = m["fourpt_momentum"]
                append!(orderparams_ed,[maximum(diag(fourpt_vals))])
            end

            append!(intstrens_ed,[m["U"][end]])

            if m["U"][end] == 0.0
                append!(laughlin_values,[orderparams_ed[end]])
            end

        end

        orderparams_ed = orderparams_ed ./ laughlin_values[n-2]

        scatter(intstrens_ed,orderparams_ed,label="N=$n")
        xlabel("Interaction Strength")
        ylabel("Max 4pt Diagonal")
        title("k-DW Order Parameter")
    end
    legend()

end

function plot_overlapslices_fourpt(Lx::Int64,ulr::Float64; kwargs...)
    plot_title::String = get(kwargs,:plot_title," ")

    if_ed::Bool = Lx <= 12
    if Lx > 10
        layers = Int(log(2,Lx*Lx/2))
        dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
        pdict_ttn = Dict([("hopping_anisotropy",1.0),("onsite_strength",ulr),("layers",layers),("if_periodic_phys",true),("if_periodic_synth",true)])
        all_files = find_data_file(pdict_ttn,"ttn",dataloc)
        display(all_files)
    else
        Ly,N = Int(Lx/2),Int(Lx/2)
        dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
        pdict = Dict([("Lx",Lx),("Ly",Ly),("N",N),("hopping_anisotropy",1.0),("interaction_strength",ulr),("if_periodic_x",true),("if_periodic_y",true)])
        all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2", output_level=0)
        filter!(x -> !occursin("twist_angle",x), all_files)
        display(all_files)
    end


    for f in all_files
        fig = figure()
        d,m = read_data(joinpath(dataloc,f); output_level=0)

        hasboth::Bool = false
        if haskey(m,"fourpt_momentum")
            if haskey(m,"fourpt_momentum_1")
                fourpt_vals = 0.5 .* (m["fourpt_momentum_1"] + m["fourpt_momentum"])
                hasboth = true
            else
                fourpt_vals = m["fourpt_momentum"]
            end
        else
            error("No fourpt data found in file $(f)")
        end

        for i in 1:Lx
            if hasboth || if_ed
                println("Working on slice $i")
                xs = collect(Int(-Lx/2+1):1:Int(Lx/2))
                reshaped_fourpts = circshift(fourpt_vals[i,:],Int(Lx/2 - i))
                plot(xs,reshaped_fourpts,"-p",label="$(i)")
            else
                if isodd(i)
                    println("Working on slice $i")
                    xs = collect(Int(-Lx/2+1):1:Int(Lx/2))
                    reshaped_fourpts = circshift(fourpt_vals[i,:],Int(Lx/2 - i))
                    plot(xs,reshaped_fourpts,"-p",label="$(i)")
                end
            end
        end
        xlabel("x")
        ylabel("Fourpt m,m+x")
        title("Fourpt Slices for ULR=$ulr")
        legend()
        title("Fourpt Slices $(Lx)x$(Int(Lx/2)) N=$(Int(Lx/2)) ULR=$ulr"*plot_title)
    end
end

# 4pt slices for ranging intstren to see change to k-DW
function plot_transition_slices(ulrs::Vector{Float64},Lx::Int64=16; kwargs...)

    if_ed::Bool = Lx <= 12
    intstren_string = if_ed ? "interacion_strength" : "onsite_strength"
    if Lx > 10
        layers = Int(log(2,Lx*Lx/2))
        dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
        pdict_ttn = Dict([("hopping_anisotropy",1.0),("layers",layers),("if_periodic_phys",true),("if_periodic_synth",true)])
        all_files = find_data_file(pdict_ttn,"ttn",dataloc)
        display(all_files)
    else
        Ly,N = Int(Lx/2),Int(Lx/2)
        dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
        pdict = Dict([("Lx",Lx),("Ly",Ly),("N",N),("hopping_anisotropy",1.0),("if_periodic_x",true),("if_periodic_y",true)])
        all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2", output_level=0)
        filter!(x -> !occursin("twist_angle",x), all_files)
        display(all_files)
    end
    which_slice = Int(Lx/2)-1

    for f in all_files
        fp = get_params_dict_from_filename(f)
        if fp[intstren_string] in ulrs
            d,m = read_data(joinpath(dataloc,f); output_level=0)

            !haskey(m,"fourpt_momentum") && error("No fourpt data found in file $(f)")
            if haskey(m,"fourpt_momentum_1")
                fourpt_vals = 0.5 .* (m["fourpt_momentum_1"] + m["fourpt_momentum"])
            else
                fourpt_vals = m["fourpt_momentum"]
            end

            xs = collect(Int(-Lx/2+1):1:Int(Lx/2))
            ys = circshift(fourpt_vals[which_slice,:],Int(Lx/2 - which_slice)) ./ sum(fourpt_vals)
            plot(xs,ys,"-p",label="$(fp[intstren_string])")

        end
    end
    xlabel("x")
    ylabel("Fourpt m,m+x")
    title("Fourpt Slices for $(Lx)x$(Int(Lx/2)) N=$(Int(Lx/2))")
    legend()
end






















"fin"
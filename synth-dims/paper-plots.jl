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
    cols = ["#82AC9F","#C73E1D","#36213E"]

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",1.0),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
    filter!(x -> !occursin("twist_angle",x), all_files)
    #display(all_files)

    g2s = []
    g3s = []
    ulrs_g2 = []
    ulrs_g3 = []
    for f in all_files
        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        all_nrgs = d["nrg"]
        
        params = get_params_dict_from_filename(f)
        ulr = params["interaction_strength"]

        if ulr > 2.0 && round(ulr,digits=0) != ulr
            continue
        end

        append!(g2s,[all_nrgs[2] - all_nrgs[1]])
        append!(ulrs_g2,[ulr])

        for i in 3:length(all_nrgs)
            append!(g3s,[all_nrgs[i] - all_nrgs[1]])
            append!(ulrs_g3,[ulr])
        end
    end
    #=fig = figure()
    scatter(ulrs_g3[1],g3s[1],c=cols[3],label="E2") # do this first to make E2 in legend on top
    scatter(ulrs_g2,g2s,c=cols[2],label="E1")
    scatter(ulrs_g3,g3s,c=cols[3]) # put them after so points are on top of E1
    scatter(ulrs_g2,zeros(length(ulrs_g2)),c=cols[1],label="E0")
    xlabel(L"U_{ir} / t")
    ylabel("E - E0")
    #title("Energy Spectrum $(lx)x$(ly) N=$(n)")
    xscale("log")
    legend()
    ylim([-0.02,0.65])=#

    # Split data for the two axis types
    log_transition_ulr = 2.1
    mask_linear_g2 = ulrs_g2 .<= log_transition_ulr
    mask_log_g2 = ulrs_g2 .> log_transition_ulr

    mask_linear_g3 = ulrs_g3 .<= log_transition_ulr
    mask_log_g3 = ulrs_g3 .> log_transition_ulr

    # Create subplots
    fig, (ax1, ax2) = subplots(1, 2, sharey=true, gridspec_kw=Dict("width_ratios" => [1, 2]))

    # Left: Linear scale
    ax1.scatter(ulrs_g3[mask_linear_g3][1], g3s[mask_linear_g3][1], c=cols[3],label="E2")
    ax1.scatter(ulrs_g2[mask_linear_g2], g2s[mask_linear_g2], c=cols[2], label="E1")
    ax1.scatter(ulrs_g3[mask_linear_g3], g3s[mask_linear_g3], c=cols[3])
    ax1.scatter(ulrs_g2[mask_linear_g2], zeros(sum(mask_linear_g2)), c=cols[1], label="E0")
    ax1.set_xlim(-0.2, log_transition_ulr)
    ax1.set_xscale("linear")
    ax1.set_ylabel("E - E₀")
    #ax1.set_xlabel(L"U_{ir} / t")

    # Right: Log scale
    ax2.scatter(ulrs_g2[mask_log_g2], g2s[mask_log_g2], c=cols[2])
    ax2.scatter(ulrs_g3[mask_log_g3], g3s[mask_log_g3], c=cols[3])
    ax2.scatter(ulrs_g2[mask_log_g2], zeros(sum(mask_log_g2)), c=cols[1])
    ax2.set_xlim(minimum(ulrs_g2[mask_log_g2]), 1.1*maximum(ulrs_g2[mask_log_g2]))
    ax2.set_xscale("log")
    ax2.set_xlabel(L"U_{ir} / t")
    # Synchronize y-limits
    ax1.set_ylim([-0.02, 0.65])

    # Formatting tweaks
    ax1.spines["right"].set_visible(false)
    ax2.spines["left"].set_visible(false)
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")

    # Optional: break marks
    d = .015
    ax1.plot([1-d, 1+d], [-d, +d], transform=ax1.transAxes, color="k", clip_on=false)
    ax1.plot([1-d, 1+d], [1-d, 1+d], transform=ax1.transAxes, color="k", clip_on=false)
    ax2.plot([-d, +d], [-d, +d], transform=ax2.transAxes, color="k", clip_on=false)
    ax2.plot([-d, +d], [1-d, 1+d], transform=ax2.transAxes, color="k", clip_on=false)

    ax2.set_title("8x4 N=4, "*L"\rho_{1D}=0.5")

    #ax1.legend()
    #ax2.legend()

end
#plot_opengap_8x4_ed()

# plot of gap staying open for rho_1D = 0.5 8x4 N=4
function plot_closedgap_8x4_ed()
    lx,ly,n = 4,8,4
    layers = Int(log(2,lx*ly))
    cols = ["#82AC9F","#C73E1D","#36213E"]

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",1.0),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
    filter!(x -> !occursin("twist_angle",x), all_files)
    #display(all_files)

    g2s = []
    g3s = []
    ulrs_g2 = []
    ulrs_g3 = []
    for f in all_files
        params = get_params_dict_from_filename(f)
        ulr = params["interaction_strength"]

        if ulr > 10.0
            continue
        end

        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        all_nrgs = d["nrg"]

        append!(g2s,[all_nrgs[2] - all_nrgs[1]])
        append!(ulrs_g2,[ulr])

        for i in 3:length(all_nrgs)
            append!(g3s,[all_nrgs[i] - all_nrgs[1]])
            append!(ulrs_g3,[ulr])
        end

        #scatter(ulr,0.0,c=cols[1])
        #scatter(ulr,all_nrgs[2] - all_nrgs[1],c=cols[2])

        #for i in 3:length(all_nrgs)
        #    scatter(ulr,all_nrgs[i] - all_nrgs[1],c=cols[3])
        #end
    end
    fig = figure()
    scatter(ulrs_g3[1],g3s[1],c=cols[3],label="E2") # do this first to make E2 in legend on top
    scatter(ulrs_g2,g2s,c=cols[2],label="E1")
    scatter(ulrs_g3,g3s,c=cols[3]) # put them after so points are on top of E1
    scatter(ulrs_g2,zeros(length(ulrs_g2)),c=cols[1],label="E0")
    xlabel(L"U_{ir} / t")
    ylabel("E - E0")
    title("$(lx)x$(ly) N=$(n), "*L"\rho_{1D}=1")
    legend()
    ylim([-0.02,0.65])


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
function plot_finitesplitting_scaling(ulr::Float64=0.0; kwargs...)
    if_plot = get(kwargs, :if_plot, true)

    cols = ["#82AC9F","#C73E1D","#36213E"]

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
        println("Working on file $f")
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
        Lx,Ly = "Lx" in keys(params) ? (params["Lx"],params["Ly"]) : get_tatami_lattice_dims(params["layers"])

        # 16x8 is not converged yet
        Lx == 16 && continue
        # 12x6 is not converged yet
        Lx == 12 && continue

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


    if if_plot
        fig = figure()
        for (k,v) in gaps
            local_label = k == "0.0" ? "Unpinned" : "$(parse(Float64,k)/10)"
            col = k == "0.0" ? cols[1] : cols[2]
            
            # temporary data manipulation for ulr = 300.0 Lx = 12
            if ulr == 300.0 && local_label == "1.0e-5"
                v[end] = 1e-5
                display(v)
            end  

            length(v) > 0 && scatter(lxs[k],v,label=local_label,c=col)
        end
        xlabel(L"L_x")
        ylabel("E1 - E0")
        title("Degeneracy Splitting "*L"U_{ir}"*"=$ulr")
        yscale("log")
        legend(loc="lower right")
    end

    return gaps,lxs
end
#rez_gaps,rez_lxs = plot_finitesplitting_scaling(300.0)

# plot finite splitting for ULR=0 and ULR=300 in same figure separate plots
function plot_paper_finitesplitting_scaling()
    gaps_laughlin,lxs_laughlin = plot_finitesplitting_scaling(0.0; if_plot=false)
    gaps_mblc,lxs_mblc = plot_finitesplitting_scaling(300.0; if_plot=false)

    cols = ["#82AC9F","#C73E1D","#36213E"]

    fig = figure()
    subplot(1,2,1)
    for (k,v) in gaps_laughlin
        local_label = k == "0.0" ? "Unpinned" : "$(parse(Float64,k)/10)"
        col = k == "0.0" ? cols[1] : cols[2]

        length(v) > 0 && scatter(lxs_laughlin[k],v,label=local_label,c=col)
    end
    xlabel(L"L_x")
    ylabel("E1 - E0")
    yscale("log")
    legend(loc="upper right")
    title("ULR=0.0")

    subplot(1,2,2)
    for (k,v) in gaps_mblc
        local_label = k == "0.0" ? "Unpinned" : "$(parse(Float64,k)/10)"
        col = k == "0.0" ? cols[1] : cols[2]
        
        # temporary data manipulation for ulr = 300.0 Lx = 12
        if local_label == "1.0e-5"
            v[end] = 1e-5
            display(v)
        end  

        length(v) > 0 && scatter(lxs_mblc[k],v,label=local_label,c=col)
    end
    xlabel(L"L_x")
    ylabel("E1 - E0")
    yscale("log")
    legend(loc="upper right")
    title("ULR=300.0")

    tight_layout()

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

    cols = ["#82AC9F","#C73E1D","#36213E"]
    
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)

    perims = [8*3,4*4,4*3,2*4,2*3]
    all_ees = []
    intstrens = [0.0,300.0]
    for f in all_files_ttn

        params = get_params_dict_from_filename(f)
        !(params["onsite_strength"] in intstrens) && continue
        if params["onsite_strength"] == 0.0
            col = cols[1]
        elseif params["onsite_strength"] == 300.0
            col = cols[2]
        end

        d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

        !("entanglement_spectrum" in keys(m)) && continue

        ees = zeros(Float64,layers-2)
        entspecs = real.(m["entanglement_spectrum"])

        for k in 2:layers-1
            entspec = filter(x -> x != 0.0, entspecs[k,:])
            #display(entspec)

            ee = entanglement_entropy(entspec)
            ees[k-1] = ee
        end
        scatter(perims,ees,label="$(params["onsite_strength"])",c=col)
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
    plot(xs,linmodel(xs,linfit.param),label="Fit: "*L"$\gamma = $"*"$(round(-yintercept, digits=3))",c=cols[3])
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

    ratios = Dict()
    intstrens = Dict()
    for f in all_files
        params = get_params_dict_from_filename(f)
        
        params["Lx"] <= 6 && continue
        params["Lx"] != 2*params["Ly"] && continue

        d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
        
        !haskey(m,"fourpt_momentum") && continue

        fourpt_vals = m["fourpt_momentum"]
        #plot_four_point(fourpt_vals; plot_title="$(params["Lx"])x$(params["Ly"]) N=$(params["N"]) ULR=$ulr")

        
        maxminratio = visibility_fourpt(m["fourpt_momentum"])

        if haskey(ratios,string(params["Lx"]))
            append!(ratios[string(params["Lx"])],[maxminratio])
            append!(intstrens[string(params["Lx"])],[m["U"][end]])
        else
            ratios[string(params["Lx"])] = [maxminratio]
            intstrens[string(params["Lx"])] = [m["U"][end]]
        end
    end

    # TTN section
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("hopping_anisotropy",1.0),("layers",7),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)

    for f in all_files_ttn
        params = get_params_dict_from_filename(f)
        if haskey(params,"Lx")
            Lx,Ly = params["Lx"],params["Ly"]
        else
            Lx,Ly = get_tatami_lattice_dims(params["layers"])
        end

        Lx != 2*Ly && continue

        d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)

        !haskey(m,"fourpt_momentum") && continue
        #plot_four_point(fourpt_vals; plot_title="$(params["Lx"])x$(params["Ly"]) N=$(params["N"]) ULR=$ulr")

        maxminratio = visibility_fourpt(m["fourpt_momentum"])

        if haskey(ratios,string(Lx))
            append!(ratios[string(Lx)],[maxminratio])
            append!(intstrens[string(Lx)],[m["onsite_strength"]])
        else
            ratios[string(Lx)] = [maxminratio]
            intstrens[string(Lx)] = [m["onsite_strength"]]
        end
    end

    # find the laughlin values
    laughlins = Dict()
    for (k,v) in intstrens
        where_zero = findfirst(x -> v[x] == 0.0, 1:length(v))
        println("Found zero for $k at index $where_zero")
        laughlins[k] = ratios[k][where_zero]
    end

    for (k,v) in ratios
        if length(v) > 0
            scatter(intstrens[k],v ./ laughlins[k],label="$(k)")
        end
    end
    
    ylabel("Visibility, min/max (normalized)")    
    xlabel("Interaction Strength")
    title("k-DW Transition Thermodynamic Finite Size Scaling")

    legend()
    #yscale("log")
    
    return ratios, intstrens
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
        pdict_ttn = Dict([("hopping_anisotropy",1.0),("onsite_strength",ulr),("Lx",Lx),("if_periodic_phys",true),("if_periodic_synth",true)])
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
    cols = ["#36213E","#82AC9F","#C73E1D"]

    if_ed::Bool = Lx <= 12
    intstren_string = if_ed ? "interaction_strength" : "onsite_strength"
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

        if haskey(fp,"Lx") && fp["Lx"] != Lx
            continue
        end
        if fp[intstren_string] in ulrs
            col = cols[findfirst(x -> fp[intstren_string] == ulrs[x],1:length(ulrs))]
            d,m = read_data(joinpath(dataloc,f); output_level=0)


            !haskey(m,"fourpt_momentum") && error("No fourpt data found in file $(f)")
            if haskey(m,"fourpt_momentum_1")
                fourpt_vals = 0.5 .* (m["fourpt_momentum_1"] + m["fourpt_momentum"])
            else
                fourpt_vals = m["fourpt_momentum"]
            end

            xs = collect(Int(-Lx/2+1):1:Int(Lx/2))
            ys = circshift(fourpt_vals[which_slice,:],Int(Lx/2 - which_slice)) ./ sum(fourpt_vals)
            plot(xs,ys,"-p",label="$(fp[intstren_string])",c=col)

        end
    end
    xlabel("x")
    ylabel(L"\hat{a}_{k}^{\dagger} \hat{a}_{k+x}^{\dagger} \hat{a}_{k+x} \hat{a}_{k}")
    #title("Fourpt Slices for $(Lx)x$(Int(Lx/2)) N=$(Int(Lx/2))")
    legend()
end
#plot_transition_slices([0.0,2.0,300.0],16; plot_title="")

# imshow of 4pt 16x8 comparing laughlin and ulr for poster
function plot_compare_fourpt(Lx::Int,Ly::Int,ulr::Float64)
    layers = Int(log(2,Lx*Ly))

    dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict = Dict([("hopping_anisotropy",1.0),("onsite_strength",ulr),("Lx",Lx),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files = find_data_file(pdict,"ttn",dataloc)

    mixed_fourpts = []
    for f in all_files
        d,m = read_data(joinpath(dataloc,f))

        intstren = m["onsite_strength"]
        if intstren != 0.0 && intstren != 300.0
            continue
        end

        plot_title = intstren == 0.0 ? "Laughlin" : L"U_{ir}"*"=$intstren"

        fourpt_vals = m["fourpt_momentum"]

        #fourpt_vals2 = m["fourpt_momentum_1"]
        fourpt_vals2 = zeros(Float64,size(fourpt_vals))
        for i in 1:size(fourpt_vals,1)-1
            for j in 1:size(fourpt_vals,2)-1
                fourpt_vals2[i,j] = fourpt_vals[i+1,j+1]
            end
        end
        fourpt_vals2[end,:] = vcat(fourpt_vals[1,:][2:end],fourpt_vals[1,:][1])
        fourpt_vals2[:,end] = vcat(fourpt_vals[:,1][2:end],fourpt_vals[:,1][1])

        mixed_fourpt = 0.25 .* (fourpt_vals .+ fourpt_vals2)

        fig = figure()
        imshow(mixed_fourpt; vmin=0.0, vmax=0.5, origin="lower")
        colorbar()
        xlabel("k")
        ylabel("k'")
        title(plot_title)

        append!(mixed_fourpts,[[mixed_fourpt]])
    end

    return mixed_fourpts
    
end

# energy spectrum where color of scatter point is the adiabatic condition matrix element with the groundstate manifold
function plot_adiabatic_spectrum(Lx::Int64)
    
    Ly,N = Int(Lx/2),Int(Lx/2)
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict = Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    display(all_files)

    f_ulrs = []
    ulrs = []
    allgaps = []

    for f in all_files
        fileparams = get_params_dict_from_filename(f)
        fileparams["interaction_strength"] > 10.0 && continue

        d,m = read_data(joinpath(dataloc,f); output_level=0)

        gaps = d["nrg"][3:end] .- d["nrg"][1]

        fuirs_0 = m["fuir_0"][2:end]
        fuirs_1 = m["fuir_1"]

        fuirs_avg = 0.5 .* (fuirs_0 .+ fuirs_1)

        append!(ulrs,[fileparams["interaction_strength"]])
        append!(f_ulrs,[[fuirs_avg]])
        append!(allgaps,[[gaps]])
    end


    clog = [log10.(f_ulrs[i][1]) for i in 1:length(f_ulrs)]

    cols = ["#82AC9F","#C73E1D","#36213E"]

    target = cols[3]
    whitetohex = matplotlib.colors.LinearSegmentedColormap.from_list(
        "white_to_hex", ["#ffffff", target]
    )

    for i in 1:length(ulrs)
        xs = ulrs[i] .* ones(length(allgaps[i][1]))
        ys = allgaps[i][1]
        scatter(xs,ys,c=clog[i],cmap=whitetohex)
    end
    colorbar().set_label(L"log_{10} (F_{U_{IR}})")
    xlabel("Interaction Strength, "*L"U_{LR}")
    ylabel("Energy Gap")
    ylim([-0.05,1.05])

    return ulrs,allgaps,clog
end

#= four point real space with new colors
if false
    lx,ly,n = 16,8,8
    intstren = 0.0
    layers = Int(log(2,lx*ly))
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus")
    pdict = Dict([("hopping_anisotropy",1.0),("onsite_strength",intstren),("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files_ttn = find_data_file(pdict,"ttn",dataloc_ttn)
    display(all_files_ttn)

    d,m = read_data(joinpath(dataloc_ttn,all_files_ttn[1]); output_level=0)
    display(keys(m))
end=#






















"fin"
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
using JLD2,LaTeXStrings,LsqFit,Dierckx

rc("font", family="serif", serif=["STIXGeneral", "Times New Roman", "Times"])
rc("mathtext", fontset="stix")

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

# plot of gap closing for rho_1D = 1.0 4x8 N=4
function plot_closedgap_8x4_ed()
    lx,ly,n = 4,8,4
    layers = Int(log(2,lx*ly))
    cols = ["#82AC9F","#C73E1D","#36213E"]

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",1.0),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
    filter!(x -> !occursin("twist_angle",x), all_files)
    #display(all_files)

    fs = 16

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
    scatter(ulrs_g3[1],g3s[1],c=cols[3],label="E2-9") # do this first to make E2 in legend on top
    scatter(ulrs_g2,g2s,c=cols[2],label="E1")
    scatter(ulrs_g3,g3s,c=cols[3]) # put them after so points are on top of E1
    scatter(ulrs_g2,zeros(length(ulrs_g2)),c=cols[1],label="E0")
    xlabel(L"U_{ir} / t",fontsize=fs)
    ylabel("E - E0",fontsize=fs)
    #title("$(lx)x$(ly) N=$(n), "*L"\rho_{1D}=1")
    legend()
    ylim([-0.02,0.65])


end

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

        length(v) > 0 && scatter(lxs_laughlin[k],v,label=local_label,c=col,marker="^")
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

        length(v) > 0 && scatter(lxs_mblc[k],v,label=local_label,c=col,marker="o")
    end
    xlabel(L"L_x")
    ylabel("E1 - E0")
    yscale("log")
    legend(loc="upper right")
    title("ULR=300.0")

    tight_layout()

end

# plot finite splitting for ULR=0 and ULR=300 in same figure separate plots
function plot_paper_finitesplitting_scaling_oneplot()
    gaps_laughlin,lxs_laughlin = plot_finitesplitting_scaling(0.0; if_plot=false)
    gaps_mblc,lxs_mblc = plot_finitesplitting_scaling(300.0; if_plot=false)

    cols = ["#82AC9F","#C73E1D","#36213E"]

    fig = figure()
    for (k,v) in gaps_laughlin
        local_label = k == "0.0" ? "Unpinned" : "$(parse(Float64,k)/10)"
        col = k == "0.0" ? cols[1] : cols[2]

        fc = k == "0.0" ? col : "none"

        length(v) > 0 && scatter(lxs_laughlin[k],v,label=local_label,marker="^",facecolors=fc,edgecolors=col)
    end
    xlabel(L"L_x")
    ylabel("E1 - E0")
    yscale("log")
    legend(loc="upper right")
    title("ULR=0.0")

    for (k,v) in gaps_mblc
        local_label = k == "0.0" ? "Unpinned" : "$(parse(Float64,k)/10)"
        col = k == "0.0" ? cols[1] : cols[2]
        
        # temporary data manipulation for ulr = 300.0 Lx = 12
        if local_label == "1.0e-5"
            v[end] = 1e-5
            display(v)
        end  

        fc = k == "0.0" ? col : "none"

        length(v) > 0 && scatter(lxs_mblc[k],v,label=local_label,marker="o",facecolors=fc,edgecolors=col)
    end
    xlabel(L"L_x")
    ylabel("E1 - E0")
    yscale("log")
    legend(loc="upper right")
    title("ULR=300.0")

    tight_layout()

end

# plot finite splitting for ULR=0 and ULR=300 in same figure separate plots
function plot_paper_finitesplitting_scaling_newplot()

    gaps_laughlin,lxs_laughlin = plot_finitesplitting_scaling(0.0; if_plot=false)
    gaps_mblc,lxs_mblc = plot_finitesplitting_scaling(300.0; if_plot=false)

    dataloc_ed = get_folder_location("cluster-data/exact-diag/torus/new-gauge/pinned-scaling")
    pdict_ed = Dict([("hopping_anisotropy",1.0),("if_pinning",true),("interaction_strength",300.0),("if_periodic_x",true),("if_periodic_y",true)])
    all_files_ed = find_data_file(pdict_ed,"ed",dataloc_ed; file_type="jld2")
    cols = ["b","g","r"]
    pins = [1e-3,1e-4,1e-5]
    datadict = Dict("0.01"=>[[],[]], "0.001"=>[[],[]], "0.0001"=>[[],[]])
    for f in all_files_ed
        d,m = read_data(joinpath(dataloc_ed,f); output_level=0)

        splitting = d["nrg"][2] - d["nrg"][1]
        Lx = m["Lx"]
        pinstren = m["pinning_strength"]

        append!(datadict[string(pinstren)][1],[Lx])
        append!(datadict[string(pinstren)][2],[splitting])

    end
    gaps_mblc["0.001"] = datadict["0.001"][2]
    lxs_mblc["0.001"] = datadict["0.001"][1]

    dataloc_ed_laughlin = get_folder_location("cluster-data/exact-diag/torus/new-gauge/pinned-scaling")
    pdict_ed_laughlin = Dict([("hopping_anisotropy",1.0),("if_pinning",true),("interaction_strength",0.0),("if_periodic_x",true),("if_periodic_y",true)])
    all_files_ed_laughlin = find_data_file(pdict_ed_laughlin,"ed",dataloc_ed_laughlin; file_type="jld2")
    datadict_laughlin = Dict("0.0001"=>[[],[]],"0.001"=>[[],[]],"0.01"=>[[],[]],"0.1"=>[[],[]])
    for f in all_files_ed_laughlin
        d,m = read_data(joinpath(dataloc_ed_laughlin,f); output_level=0)

        splitting = d["nrg"][2] - d["nrg"][1]
        Lx = m["Lx"]

        !haskey(m,"pinning_strength") && continue

        append!(datadict_laughlin[string(m["pinning_strength"])][1],[Lx])
        append!(datadict_laughlin[string(m["pinning_strength"])][2],[splitting])

    end
    gaps_laughlin["0.001"] = datadict_laughlin["0.001"][2]
    lxs_laughlin["0.001"] = datadict_laughlin["0.001"][1]

    cols = ["#82AC9F","#C73E1D","#36213E"]

    fig, axs = subplots(1,2; figsize=(6,4))
    
    axs[1].scatter(lxs_laughlin["0.0"],gaps_laughlin["0.0"],label=L"U_{ir}=0",c=cols[2],marker="^")
    axs[1].scatter(lxs_mblc["0.0"],gaps_mblc["0.0"],label=L"U_{ir}=300",c=cols[1],marker="o")
    axs[1].set_xlabel(L"L_x",fontsize=16)
    axs[1].set_ylabel(L"\Delta_{12}",fontsize=16)
    axs[1].set_yscale("log")
    axs[1].legend(loc="center right",fontsize=12)
    axs[1].set_title("Without Pinning",fontsize=16)

    axs[1].tick_params(axis="both", which="major", labelsize=14)

    ticks = [6,7,8,9,10,11,12,13,14]
    labels = ["6","","8","","10","","12","","14"]
    axs[1].set_xticks(ticks,labels)

    for (k,v) in gaps_laughlin
        local_label = k == "0.0" ? "Unpinned" : "$(parse(Float64,k)/10)"
        
        if k == "0.0"
            continue
        end


        if length(v) > 0
            if local_label == "0.0001"
                local_label = "1.0e-4"
                axs[2].scatter(lxs_laughlin[k],v,label=L"\delta"*"=$(local_label)",marker="^",facecolors="none",edgecolors=cols[2])
            elseif k == "0.1"
                axs[2].scatter(lxs_laughlin[k],v,label=L"\delta"*"=$(local_label)",c=cols[2],marker="^")
            end
        end
    end
    for (k,v) in gaps_mblc
        local_label = k == "0.0" ? "Unpinned" : "$(parse(Float64,k)/10)"
        
        if k == "0.0"
            continue
        end

        # temporary data manipulation for ulr = 300.0 Lx = 12
        if local_label == "1.0e-5"
            v[end] = 1e-5
            display(v)
        end

        if length(v) > 0
            if local_label == "0.0001"
                local_label = "1.0e-4"
                axs[2].scatter(lxs_mblc[k],v,label=L"\delta"*"=$(local_label)",marker="o",facecolors="none",edgecolors=cols[1])
            elseif k == "0.0001"
                axs[2].scatter(lxs_mblc[k],v,label=L"\delta"*"=$(local_label)",c=cols[1],marker="o")
            end
        end
    end
    axs[2].legend(loc="upper right",fontsize=12)
    axs[2].set_xlabel(L"L_x", fontsize=16)
    #axs[2].set_ylabel("GS Splitting",fontsize=16)
    axs[2].set_yscale("log")
    axs[2].set_title("With Pinning",fontsize=16)

    axs[2].tick_params(axis="both", which="major", labelsize=14)
    axs[2].set_xticks(ticks,labels)

    axs[1].tick_params(top=true, right=true)
    axs[1].tick_params(labeltop=false, labelright=false)
    axs[2].tick_params(top=true, right=true)
    axs[2].tick_params(labeltop=false, labelright=false)

    yticks_pin = [ii*10.0^(-jj) for ii in 1:9 for jj in 1:6]
    ylabels_pin = ["" for ii in 1:9 for jj in 1:6]
    ylabels_pin[1] = L"10^{-1}"
    ylabels_pin[2] = L"10^{-2}"
    ylabels_pin[3] = L"10^{-3}"
    ylabels_pin[4] = L"10^{-4}"
    ylabels_pin[5] = L"10^{-5}"
    axs[2].set_yticks(yticks_pin, ylabels_pin)

    axs[2].set_ylim(5e-6,3e-1)

    
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
            col = cols[2]
            mm = "^"
        elseif params["onsite_strength"] == 300.0
            col = cols[1]
            mm = "o"
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
        scatter(perims,ees,label="$(params["onsite_strength"])",c=col,marker=mm)
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

# four point slices with both 0 and 300 and correct markers/colors
function plot_four_point_slices_both()
    cols = ["#82AC9F","#C73E1D","#36213E"]

    d_300,m_300 = read_data("../cluster-data/synth-dims/torus/new-gauge/ttn-if_periodic_phys-true-onsite_strength-300.0-lr-7-particles-8-alpha-0.125-if_periodic_synth-true-layers-7-hopping_anisotropy-1.0.h5");
    fourpt_1_300 = m_300["fourpt_momentum"]
    
    d_0,m_0 = read_data("../cluster-data/synth-dims/torus/new-gauge/ttn-if_periodic_phys-true-onsite_strength-0.0-lr-7-particles-8-alpha-0.125-if_periodic_synth-true-layers-7-hopping_anisotropy-1.0.h5");
    fourpt_1_0 = m_0["fourpt_momentum"]
    fourpt_2_0 = m_0["fourpt_momentum_1"]
    mixed_fourpt = 0.5 .* (fourpt_1_0 .+ fourpt_2_0)
    
    xs = collect(-7:7)

    ys_0 = circshift(mixed_fourpt[1,:],7)[1:end-1]
    
    ys_300 = circshift(fourpt_1_300[7,:],1)[1:end-1] ./ 2


    fig = figure(figsize=(6,4))
    plot(xs,ys_0,c=cols[2],label=L"U_{ir}=0",marker="^")
    plot(xs,ys_300,c=cols[1],label=L"U_{ir}=300",marker="o")
    legend(loc="upper center",fontsize=14)
    xlabel(L"k_y", fontsize=18)
    #ylabel(L"\langle \hat{a}_{k_y}^{\dagger} \hat{a}_{k_{y}^{'}}^{\dagger} \hat{a}_{k_{y}^{'}} \hat{a}_{k_{y}} \rangle / \langle \hat{n}_{k_{y}^{'}} \rangle \langle \hat{n}_{k_{y}} \rangle",fontsize=14)
    ylabel(L"C^{(4)} (k_y, k_y'=0)",fontsize=16)
    tick_params(axis="both", which="major", labelsize=14)

    xtickvals = collect(-7:7)
    xticklabels = ["","-6","","-4","","-2","","0","","2","","4","","6",""]
    xticks(xtickvals,xticklabels)

    tick_params(top=true, right=true)
    tick_params(labeltop=false, labelright=false)

    fig.tight_layout()
end

# energy spectrum where color of scatter point is the adiabatic condition matrix element with the groundstate manifold
function plot_ULR_adiabatic_spectrum(Lx::Int64)
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

    fs = 15


    clog = [log10.(f_ulrs[i][1]) for i in 1:length(f_ulrs)]

    cols = ["#82AC9F","#C73E1D","#36213E","#B47EB3"]

    target = cols[2]
    starting = cols[1]
    whitetohex = matplotlib.colors.LinearSegmentedColormap.from_list(
        "c1_to_c2", ["#ffffff", target]
    )

    for i in 1:length(ulrs)
        xs = ulrs[i] .* ones(length(allgaps[i][1]))
        ys = allgaps[i][1]
        scatter(xs,ys,c=clog[i],cmap=whitetohex)
    end
    colorbar().set_label(L"log_{10} (F_{U_{i}})",fontsize=fs)
    xlabel("Interaction Strength, "*L"U_{i}",fontsize=fs)
    ylabel("Energy Gap",fontsize=fs)
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

# TEE as a function of ULR if there are enough data points
function plot_tee_vs_ulr()
    cols = ["#82AC9F","#C73E1D","#36213E"]

    lx,ly,n = 16,8,8
    layers = 7
    dataloc1 = get_folder_location("cluster-data/synth-dims/torus")
    pdict = Dict([("hopping_anisotropy",1.0),("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files1 = find_data_file(pdict,"ttn",dataloc1)
    
    dataloc2 = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    all_files2 = find_data_file(pdict,"ttn",dataloc2)

    all_locs = (dataloc2*"/") .* all_files2
    
    linmodel(x,p) = p[1] .* x .+ p[2]
    perims = [8*3,4*4,4*3,2*4,2*3]

    yints = []
    sigmas = []
    ulrs = []
    for (idx,f) in enumerate(all_locs)
        d,m = read_data(f; output_level=0)

        if haskey(m,"entanglement_spectrum")
            println("Processing file $idx / $(length(all_locs))")
            params = get_params_dict_from_filename(f)
            intstren = haskey(m,"onsite_strength") ? m["onsite_strength"] : params["onsite_strength"]

            ees = zeros(Float64,layers-2)
            entspecs = m["entanglement_spectrum"]

            for k in 2:layers-1
                entspec = filter(x -> x != 0.0, entspecs[:,k-1])
                #display(entspec)

                ee = entanglement_entropy(entspec)
                ees[k-1] = ee
            end

            linfit = curve_fit(linmodel,perims,ees,[1.0,1.0])
            sigma = stderror(linfit)[2]
            yintercept = -linfit.param[2]
            append!(yints,[yintercept])
            append!(sigmas,[sigma])
            append!(ulrs,[intstren])

            fig = figure()
            scatter(perims,ees)
            title("ULR=$intstren")
            xlim([0,1.1*maximum(perims)])
            ylim([-1,1.1*maximum(ees)])
        end
    end

    fig = figure()
    errorbar(ulrs,yints,yerr=sigmas,fmt="o",c=cols[3])
    plot(range(0,1000,length=10),0.5 .* ones(10),"--",c=cols[2],label="Laughlin: "*L"\gamma = 1/2")
    legend(fontsize=12)

    fs = 16

    xscale("symlog"; linthresh=2.0, linscale=1.0)
    ymin,ymax = -0.1,1.1
    ylim([ymin, ymax])
    xmin,xmax = -0.25, 1200.0
    xlim([xmin,xmax])
    xlabel(L"U_{\mathrm{i}} / t", fontsize=fs)
    ylabel(L"\gamma", fontsize=fs)

    ticks = [0,0.5,1,1.5,2,10,100,1000]
    labels = ["0","0.5","1","1.5","2",L"10^1",L"10^2",L"10^3"]
    xticks(ticks,labels)

    tight_layout()

end


function get_chern_number_plot(intstren::Float64; if_spline::Bool=false, cutoff::Float64=0.27)
    lx,ly,n = 8,4,4
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    display(all_files)

    d,m = read_data(joinpath(dataloc,all_files[1]); output_level=0)

    tw1s = m["tw1s"]
    tw2s = m["tw2s"]
    omegas = m["omegas"]
    lambda1 = m["lambda1s"]


    which_to_keep = []
    for (idx,tw1) in enumerate(tw1s)
        if round(tw1,digits=1) == tw1 && round(tw2s[idx],digits=1) == tw2s[idx]
            append!(which_to_keep,[idx])
        end
    end

    squaresize = Int(sqrt(length(which_to_keep)))
    new_tw1s = unique(tw1s[which_to_keep]) 
    new_tw2s = unique(tw2s[which_to_keep]) .- 0.5
    new_wrongshape_omegas = reshape(omegas[which_to_keep],squaresize,squaresize)
    new_omegas = zeros(ComplexF64,squaresize,squaresize)
    new_omegas[:,1:5] = new_wrongshape_omegas[:,7:end]
    new_omegas[:,6:end] = new_wrongshape_omegas[:,1:6]
    rez = plot_omega(new_tw1s,new_tw2s,new_omegas; plot_title=L"U_{ir}"*"=$intstren",if_plot=false)

    new_wrongshape_lambda1 = reshape(lambda1[which_to_keep],squaresize,squaresize)
    new_lambda1 = zeros(ComplexF64,squaresize,squaresize)
    new_lambda1[:,1:5] = new_wrongshape_lambda1[:,7:end]
    new_lambda1[:,6:end] = new_wrongshape_lambda1[:,1:6]
    plotgamma1 = plot_gamma(new_tw1s,new_tw2s,new_lambda1,1; plot_title="ULR=$intstren", if_plot=false)
    spline1 = get_spline_outline(plotgamma1,cutoff,range(0,1,length=11),range(0,1,length=11))

    chern_count,chern_imshow,chern_xs,chern_ys = count_chern_number(new_tw1s,new_tw2s,rez; if_plot=false)

    return chern_imshow,chern_xs,chern_ys,spline1
end

function get_spline_outline(A::AbstractMatrix,
                            cutoff::Real,
                            xcoords::AbstractVector,
                            ycoords::AbstractVector)

    fig, ax = subplots()

    # Show image
    ax.imshow(A; interpolation="nearest",
              extent=(minimum(xcoords), maximum(xcoords),
                      minimum(ycoords), maximum(ycoords)),
              origin="lower")

    # Contour at the actual cutoff level
    cs = ax.contour(xcoords, ycoords, A; levels=[cutoff], colors="none")

    paths = cs[:collections][1][:get_paths]()

    all_xs = Vector{Vector{Float64}}()
    all_ys = Vector{Vector{Float64}}()

    for path in paths
        verts = Array(path[:vertices])   # N×2
        if size(verts, 1) < 10
            continue
        end

        xs = verts[:, 1]
        ys = verts[:, 2]

        # Smooth the outline with splines
        n  = length(xs)
        t  = range(0, 1; length=n)

        sx = Spline1D(t, xs; k=3, s=0.0)
        sy = Spline1D(t, ys; k=3, s=0.0)

        tt  = range(0, 1; length=500)
        xsp = sx.(tt)
        ysp = sy.(tt)

        push!(all_xs, xsp)
        push!(all_ys, ysp)
    end

    close(fig)
    return all_xs, all_ys
end

# plot Matteo Hatsugai step-by-step End Matter
function plot_hatsugai_stepbystep()
    lx,ly,n = 8,4,4
    intstren = 1000.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    display(all_files)

    d,m = read_data(joinpath(dataloc,all_files[1]); output_level=0)

    tw1s = m["tw1s"]
    tw2s = m["tw2s"]
    omegas = m["omegas"]
    lambda1 = m["lambda1s"]
    lambda2 = m["lambda2s"]

    which_to_keep = []
    for (idx,tw1) in enumerate(tw1s)
        if round(tw1,digits=1) == tw1 && round(tw2s[idx],digits=1) == tw2s[idx]
            append!(which_to_keep,[idx])
        end
    end

    squaresize = Int(sqrt(length(which_to_keep)))
    new_tw1s = unique(tw1s[which_to_keep])
    new_tw2s = unique(tw2s[which_to_keep]) .- 0.5
    new_wrongshape_omegas = reshape(omegas[which_to_keep],squaresize,squaresize)
    new_omegas = zeros(ComplexF64,squaresize,squaresize)
    new_omegas[:,1:5] = new_wrongshape_omegas[:,7:end]
    new_omegas[:,6:end] = new_wrongshape_omegas[:,1:6]

    new_wrongshape_lambda1 = reshape(lambda1[which_to_keep],squaresize,squaresize)
    new_lambda1 = zeros(ComplexF64,squaresize,squaresize)
    new_lambda1[:,1:5] = new_wrongshape_lambda1[:,7:end]
    new_lambda1[:,6:end] = new_wrongshape_lambda1[:,1:6]
    new_wrongshape_lambda2 = reshape(lambda2[which_to_keep],squaresize,squaresize)
    new_lambda2 = zeros(ComplexF64,squaresize,squaresize)
    new_lambda2[:,1:5] = new_wrongshape_lambda2[:,7:end]
    new_lambda2[:,6:end] = new_wrongshape_lambda2[:,1:6]

    fig, axs = subplots(1,4; figsize=(13,3))
    cutoffval = 0.27

    xticks = [0,0.25,0.5,0.75,1.0]
    xlabels = ["0.0","0.25","0.5","0.75","1.0"]

    yticks = [-0.5,-0.25,0,0.25,0.5]
    ylabels = ["-0.5","-0.25","0.0","0.25","0.5"]

    fs = 14

    plotgamma2 = plot_gamma(new_tw1s,new_tw2s,new_lambda2,2; plot_title="ULR=$intstren", if_plot=false)
    im1 = axs[1].imshow(plotgamma2; extent=(minimum(new_tw1s), maximum(new_tw1s), minimum(new_tw2s), maximum(new_tw2s)), vmin=0.0, vmax=1.0, origin="lower")
    axs[1].set_ylabel(L"\theta_y / 2\pi",fontsize=fs)
    axs[1].set_xlabel(L"\theta_x / 2\pi",fontsize=fs)
    fig.colorbar(im1, ax=axs[1])
    axs[1].set_title(L"\vert \Lambda_{\Phi} \vert",fontsize=fs+2)
    spline2 = get_spline_outline(plotgamma2,cutoffval,new_tw1s,new_tw2s)
    for i in 1:length(spline2[1])
        axs[1].plot(spline2[1][i],spline2[2][i],"-m",linewidth=2)
    end
    axs[1].set_xticks(xticks,xlabels)
    axs[1].set_yticks(yticks,ylabels)

    plotgamma1 = plot_gamma(new_tw1s,new_tw2s,new_lambda1,1; plot_title="ULR=$intstren", if_plot=false)
    im2 = axs[2].imshow(plotgamma1; extent=(minimum(new_tw1s), maximum(new_tw1s), minimum(new_tw2s), maximum(new_tw2s)), vmin=0.0, vmax=1.0, origin="lower")
    axs[2].set_ylabel(L"\theta_y / 2\pi",fontsize=fs)
    axs[2].set_xlabel(L"\theta_x / 2\pi",fontsize=fs)
    fig.colorbar(im2, ax=axs[2])
    axs[2].set_title(L"\vert \Lambda_{\Phi'} \vert",fontsize=fs+2)
    spline1 = get_spline_outline(plotgamma1,cutoffval,new_tw1s,new_tw2s)
    for i in 1:length(spline1[1])
        axs[2].plot(spline1[1][i],spline1[2][i],"-k",linewidth=2)
    end
    axs[2].set_xticks(xticks,xlabels)
    axs[2].set_yticks(yticks,ylabels)

    plotomega = plot_omega(new_tw1s,new_tw2s,new_omegas; plot_title=L"U_{ir}"*"=$intstren", if_plot=false, if_count_chern=false)
    im3 = axs[3].imshow(plotomega; extent=(minimum(new_tw1s), maximum(new_tw1s), minimum(new_tw2s), maximum(new_tw2s)), vmin=0, vmax=2*pi, cmap="hsv", origin="lower")
    xs = transpose(repeat(new_tw1s,1,length(new_tw2s)))
    ys = reverse(repeat(new_tw2s,1,length(new_tw1s)),dims=1)
    us = cos.(plotomega)
    vs = sin.(plotomega)
    axs[3].quiver(xs, ys, us, vs)
    axs[3].set_ylabel(L"\theta_y / 2\pi",fontsize=fs)
    axs[3].set_xlabel(L"\theta_x / 2\pi",fontsize=fs)
    fig.colorbar(im3, ax=axs[3])
    axs[3].set_title(L"\Omega_{\Phi \rightarrow \Phi'}",fontsize=fs+2)
    axs[3].set_xticks(xticks,xlabels)
    axs[3].set_yticks(yticks,ylabels)


    aa,plotcherns,bb,cc = count_chern_number(new_tw1s,new_tw2s,plotomega; if_plot=false)
    im4 = axs[4].imshow(plotcherns ./ (2*pi); extent=(minimum(new_tw1s), maximum(new_tw1s), minimum(new_tw2s), maximum(new_tw2s)), vmin=-1, vmax=1, cmap="bwr", origin="lower")
    axs[4].set_ylabel(L"\theta_y / 2\pi",fontsize=fs)
    axs[4].set_xlabel(L"\theta_x / 2\pi",fontsize=fs)
    fig.colorbar(im4, ax=axs[4])
    axs[4].set_title(L"\frac{1}{2\pi} \sum_{\diamond} \Delta \Omega",fontsize=fs+2)
    for i in 1:length(spline1[1])
        axs[4].plot(spline1[1][i],spline1[2][i],"-k",linewidth=2)
    end
    for i in 1:length(spline2[1])
        axs[4].plot(spline2[1][i],spline2[2][i],"-m",linewidth=2)
    end
    axs[4].set_xticks(xticks,xlabels)
    axs[4].set_yticks(yticks,ylabels)

    tight_layout()

end

# plot gap spectrum with Chern insets
function plot_opengap_spectrum_with_chern()
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


    fig = figure()

    fs = 16

    x1_1000 = 290.0
    y_1000 = 0.108
    x0_1000 = 60.0
    plot([1000.0,x1_1000],[0.0,y_1000],c="gray",linewidth=1)
    plot([1000.0,x0_1000],[0.0,y_1000],c="gray",linewidth=1)

    x1_0 = 2.4
    y_0 = 0.108
    x0_0 = 0.92
    plot([0.0,x1_0],[0.0,y_0],c="gray",linewidth=1)
    plot([0.0,x0_0],[0.0,y_0],c="gray",linewidth=1)

    scatter(ulrs_g3,g3s,c=cols[3],label="E2-9")
    scatter(ulrs_g2,g2s,c=cols[2],label="E1")
    scatter(ulrs_g2,zeros(length(ulrs_g2)),c=cols[1],label="E0")

    xscale("symlog"; linthresh=2.0, linscale=1.0)
    ymin,ymax = -0.02,0.5
    ylim([ymin, ymax])
    xmin,xmax = -0.25, 1200.0
    xlim([xmin,xmax])
    xlabel(L"U_{\mathrm{i}} / t", fontsize = fs)
    ylabel(L"E-E_0", fontsize = fs)
    legend(loc="upper left")

    ticks = [0,0.5,1,1.5,2,10,100,1000]
    labels = ["0","0.5","1","1.5","2",L"10^1",L"10^2",L"10^3"]
    xticks(ticks,labels)

    #tight_layout()

    height = 0.085
    renorm_height = height / (ymax - ymin)

    # get chern number calcs
    chern_imshow_1000,chern_xs_1000,chern_ys_1000,spline1 = get_chern_number_plot(1000.0; if_spline=true,cutoff=0.28)
    chern_imshow_0 = chern_imshow_1000

    ax_inset1 = PyPlot.axes([0.25, 0.3, renorm_height, renorm_height])
    im1 = ax_inset1.imshow(chern_imshow_0 ./ (2*pi), extent=(0,1,0,1), cmap="bwr", aspect="auto", vmin=-1, vmax=1, origin="lower")
    fig.colorbar(im1, ax=ax_inset1)
    title("Winding Defects", fontsize=8)
    for i in 1:length(spline1[1])
        ax_inset1.plot(spline1[1][i],spline1[2][i],"-k",linewidth=2)
    end

    ax_inset2 = PyPlot.axes([0.65, 0.3, renorm_height, renorm_height])
    im2 = ax_inset2.imshow(chern_imshow_1000 ./ (2*pi), extent=(0,1,0,1), cmap="bwr", aspect="auto", vmin=-1, vmax=1, origin="lower")
    fig.colorbar(im2, ax=ax_inset2)
    title("Winding Defects", fontsize=8)
    for i in 1:length(spline1[1])
        ax_inset2.plot(spline1[1][i],spline1[2][i],"-k",linewidth=2)
    end

end

# combines spectrum with Chern inset with TEE vs ULR
function make_topomarkers_paperplot()

    # top figure of spectrum with Chern insets

    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",1.0),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
    filter!(x -> !occursin("twist_angle",x), all_files)

    g2s = []
    g3s = []
    ulrs_g2 = []
    ulrs_g3 = []
    g1s = []
    ulrs_g1 = []
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

        #append!(g1s,[all_nrgs[1]])
        #append!(ulrs_g1,[ulr])

        for i in 3:length(all_nrgs)
            append!(g3s,[all_nrgs[i] - all_nrgs[1]])
            append!(ulrs_g3,[ulr])
        end
    end#
    

    fig = figure(figsize=(6, 7), constrained_layout=true)
    gs  = fig.add_gridspec(2, 1, height_ratios=[2, 1])  # top:bottom = 2:1
    ax_top = fig.add_subplot(gs[0, 0])
    ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)
    axs = [ax_bot, ax_top]

    mm = "_"
    ss = 40


    cols = ["#82AC9F","#C73E1D","#36213E"]
    fs = 20

    x1_1000 = 340.0
    y_1000 = 0.108
    x0_1000 = 75.0
    axs[1].plot([1000.0,x1_1000],[0.0,y_1000],c="gray",linewidth=1)
    axs[1].plot([1000.0,x0_1000],[0.0,y_1000],c="gray",linewidth=1)

    x1_0 = 2.05
    y_0 = 0.108
    x0_0 = 0.83
    axs[1].plot([0.0,x1_0],[0.0,y_0],c="gray",linewidth=1)
    axs[1].plot([0.0,x0_0],[0.0,y_0],c="gray",linewidth=1)

    axs[1].scatter(ulrs_g3,g3s,c=cols[3],marker=mm,s=ss)
    axs[1].scatter(ulrs_g2,g2s,c=cols[2],marker=mm,s=ss)
    axs[1].scatter(ulrs_g2,zeros(length(ulrs_g2)),c=cols[1],marker=mm,s=ss)

    axs[1].set_xscale("symlog"; linthresh=2.0, linscale=1.0)
    axs[2].set_xscale("symlog"; linthresh=2.0, linscale=1.0)
    ymin,ymax = -0.02,0.5
    axs[1].set_ylim([ymin, ymax])
    xmin,xmax = -0.25, 1200.0
    axs[1].set_xlim([xmin,xmax])
    #axs[1].set_xlabel(L"U_{\mathrm{i}} / t", fontsize = fs)
    axs[1].set_ylabel(L"E-E_0", fontsize = fs)

    axs[2].set_xlabel(L"U_{\mathrm{i}} / t", fontsize=fs)
    axs[2].set_ylabel(L"\gamma", fontsize=fs)

    ticks = [0,0.5,1,1.5,2,3,4,5,6,7,8,9,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000]
    labels = ["0","0.5","1","1.5","2","","","","","","","",L"10^1","","","","","","","","",L"10^2","","","","","","","","",L"10^3"]
    axs[1].set_xticks(ticks,labels)
    axs[2].set_xticks(ticks,labels)

    gammaticks = collect(0:8) ./ 8
    gammalabels = ["0.0","","0.25","","0.5","","0.75","","1.0"]
    axs[2].set_yticks(gammaticks,gammalabels)

    axs[2].tick_params(axis="both", which="major", labelsize=16)
    axs[1].tick_params(axis="both", which="major", labelsize=16)

    # add ticks to top and left
    axs[1].tick_params(top=true, right=true)
    axs[1].tick_params(labeltop=false, labelright=false)

    # add ticks to top and left
    axs[2].tick_params(top=true, right=true)
    axs[2].tick_params(labeltop=false, labelright=false)

    axs[2].plot(range(0,1000,length=10),0.5 .* ones(10),"--",c=cols[2],label="Laughlin")

    fig.tight_layout()

    height = 0.085
    renorm_height = height / (ymax - ymin)

    # get chern number calcs
    chern_imshow_1000,chern_xs_1000,chern_ys_1000,spline1 = get_chern_number_plot(1000.0; if_spline=true,cutoff=0.28)
    chern_imshow_0 = chern_imshow_1000

    ax_inset1 = axs[1].inset_axes([0.15, 0.25, renorm_height, renorm_height*(8/6)])
    im1 = ax_inset1.imshow(chern_imshow_0 ./ (2*pi), extent=(0,1,0,1), cmap="bwr", aspect="auto", vmin=-1, vmax=1, origin="lower")
    inset1_colorbar = fig.colorbar(im1; ax=ax_inset1, anchor=(3,0.5))
    ax_inset1.set_title("Winding Defects", fontsize=12)
    for i in 1:length(spline1[1])
        ax_inset1.plot(spline1[1][i],spline1[2][i],"-k",linewidth=2)
    end

    ax_inset2 = axs[1].inset_axes([0.7, 0.25, renorm_height, renorm_height*(8/6)])
    im2 = ax_inset2.imshow(chern_imshow_1000 ./ (2*pi), extent=(0,1,0,1), cmap="bwr", aspect="auto", vmin=-1, vmax=1, origin="lower")
    inset2_colorbar = fig.colorbar(im2; ax=ax_inset2, anchor=(3,0.5))
    ax_inset2.set_title("Winding Defects", fontsize=12)
    for i in 1:length(spline1[1])
        ax_inset2.plot(spline1[1][i],spline1[2][i],"-k",linewidth=2)
    end
    
    inset1_colorbar.ax.tick_params(labelsize=12)
    inset2_colorbar.ax.tick_params(labelsize=12)

    insetticks = [0,0.5,1.0]
    insetlabels = ["0","0.5","1"]
    ax_inset1.set_yticks(insetticks,insetlabels)
    ax_inset1.set_xticks(insetticks,insetlabels)
    ax_inset2.set_yticks(insetticks,insetlabels)
    ax_inset2.set_xticks(insetticks,insetlabels)

    ax_inset2.tick_params(axis="both", which="major", labelsize=12)
    ax_inset1.tick_params(axis="both", which="major", labelsize=12)

    ax_inset1.tick_params(top=true, right=true)
    ax_inset1.tick_params(labeltop=false, labelright=false)

    ax_inset2.tick_params(top=true, right=true)
    ax_inset2.tick_params(labeltop=false, labelright=false)

    ax_inset2.set_xlabel(L"\theta_x", fontsize=12)
    ax_inset2.set_ylabel(L"\theta_y", fontsize=12)
    ax_inset1.set_xlabel(L"\theta_x", fontsize=12)
    ax_inset1.set_ylabel(L"\theta_y", fontsize=12)


    ##### Second figure of TEE vs ULR
    #
    lx,ly,n = 16,8,8
    layers = 7
    dataloc1 = get_folder_location("cluster-data/synth-dims/torus")
    pdict = Dict([("hopping_anisotropy",1.0),("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files1 = find_data_file(pdict,"ttn",dataloc1)
    
    dataloc2 = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    all_files2 = find_data_file(pdict,"ttn",dataloc2)

    all_locs = (dataloc2*"/") .* all_files2
    
    linmodel(x,p) = p[1] .* x .+ p[2]
    perims = [8*3,4*4,4*3,2*4,2*3]

    yints = []
    sigmas = []
    ulrs = []
    for (idx,f) in enumerate(all_locs)
        d,m = read_data(f; output_level=0)

        if haskey(m,"entanglement_spectrum")
            println("Processing file $idx / $(length(all_locs))")
            params = get_params_dict_from_filename(f)
            instren = haskey(m,"onsite_strength") ? m["onsite_strength"] : params["onsite_strength"]

            ees = zeros(Float64,layers-2)
            entspecs = real.(m["entanglement_spectrum"])

            for k in 2:layers-1
                entspec = filter(x -> x != 0.0, entspecs[:,k-1])
                #display(entspec)

                ee = entanglement_entropy(entspec)
                ees[k-1] = ee
            end
            linfit = curve_fit(linmodel,perims,ees,[1.0,1.0])
            sigma = stderror(linfit)[2]
            yintercept = -linfit.param[2]
            append!(yints,[yintercept])
            append!(sigmas,[sigma])
            append!(ulrs,[instren])
        end
    end#

    axs[2].errorbar(ulrs,yints,yerr=sigmas,fmt="o",c=cols[3])

    ymin,ymax = 0.0,1.0
    axs[2].set_ylim([ymin, ymax])
    xmin,xmax = -0.25, 1200.0
    axs[2].set_xlim([xmin,xmax])

    axs[2].legend(loc="upper left", fontsize=14)
end

# energy spectrum under twisting plot
function plot_spectrum_twisting3D(lx::Int,n::Int,intstren::Float64)
    ly = lx == n ? n*2 : n
    dataloc = intstren == 0.0 ? get_folder_location("cluster-data/exact-diag/torus") : get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",1.0),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true)])
    all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2") 
    display(all_files)

    d,m = read_data(joinpath(dataloc,all_files[1]); output_level=0)

    nrgs0 = m["twisted_nrgs_0"]
    nrgs1 = m["twisted_nrgs_1"]
    nrgs2 = m["twisted_nrgs_2"]

    twists = range(0,1,length=size(nrgs0,1))

    cols = ["#82AC9F","#C73E1D","#36213E"]

    nx = length(twists)
    ny = length(twists)
    xx = repeat(reshape(twists, 1, nx), ny, 1)
    yy = repeat(reshape(twists, ny, 1), 1, nx)

    fig = figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(xx,yy,nrgs0,color=cols[1])
    ax.plot_surface(xx,yy,nrgs1,color=cols[2])
    ax.plot_surface(xx,yy,nrgs2,color=cols[3])

    xlabel(L"\theta_x",fontsize=16)
    ylabel(L"\theta_y",fontsize=16)
    zlabel("Energy",fontsize=16)
    title("$(lx)x$(ly) N=$(n) ULR=$(intstren)",fontsize=16)

end

# nrg spectra under twisting 3D for 8x4 0, 8x4 1000, 4x8 1000
function plot_spectrum_twisting3D_allthree()
    fig = figure()
    ax = [fig.add_subplot(1, 3, i, projection="3d") for i in 1:3]
    idx = 1
    for (lx,n,intstren) in [(8,4,0.0),(8,4,1000.0),(4,4,1000.0)]
        ly = lx == n ? n*2 : n
        dataloc = intstren == 0.0 ? get_folder_location("cluster-data/exact-diag/torus") : get_folder_location("cluster-data/exact-diag/torus/new-gauge")
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",1.0),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true)])
        all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2") 
        display(all_files)

        d,m = read_data(joinpath(dataloc,all_files[1]); output_level=0)

        nrgs0 = m["twisted_nrgs_0"]
        nrgs1 = m["twisted_nrgs_1"]
        nrgs2 = m["twisted_nrgs_2"]

        twists = range(0,1,length=size(nrgs0,1))

        cols = ["#82AC9F","#C73E1D","#36213E"]

        nx = length(twists)
        ny = length(twists)
        xx = repeat(reshape(twists, 1, nx), ny, 1)
        yy = repeat(reshape(twists, ny, 1), 1, nx)

        s1 = ax[idx].plot_surface(xx,yy,nrgs0,color=cols[1])
        s2 = ax[idx].plot_surface(xx,yy,nrgs1,color=cols[2])
        s3 = ax[idx].plot_surface(xx,yy,nrgs2,color=cols[3])

        ax[idx].set_title("$(lx)x$(ly), N=$(n), "*L"U_{\mathrm{i}}="*"$(intstren)",fontsize=12)
        
        if idx == 1
            ax[idx].set_xlabel(L"\theta_x",fontsize=16)
            ax[idx].set_ylabel(L"\theta_y",fontsize=16)
            ax[idx].set_zlabel("Energy",fontsize=16)
            #fig.legend(handles=[s1, s2, s3],labels=["E0", "E1", "E2"],loc="upper center")
        end
        idx += 1
    end
end

# pair density distribution vs ULR for 16x8
function plot_realspace_pair_dist_stacked()
    lx,ly,n = 16,8,8
    dataloc = get_folder_location("cluster-data/synth-dims/torus")
    pdict = Dict([("hopping_anisotropy",1.0),("layers",7),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    display(all_files)

    idx = 1
    
    fig, axs = subplots(3,1; sharex = true,sharey = true,figsize = (4, 6),constrained_layout = true)

    for f in all_files
        d,m = read_data(joinpath(dataloc,f); output_level=0)
        
        if haskey(m,"densitydensity")
            ulr = get_params_dict_from_filename(f)["onsite_strength"]
            dd = m["densitydensity"]
            occs = transpose(get_occupancy(d["densmat"]; if_plot=false))
            pairdist = pairdistribution(dd,occs; if_plot=false)
            
            thisim = axs[idx].imshow(pairdist, origin="lower",vmin=0.0,vmax=1.6,extent=[1,16,1,8])
            axs[idx].set_title(L"U_{\mathrm{i}}="*"$ulr",fontsize=10)
            if idx == 3
                axs[idx].set_xlabel("Physical",fontsize=10)
                axs[idx].set_ylabel("Synthetic",fontsize=10)
                fig.colorbar(thisim, ax=axs[idx])
            end

            global idx += 1
        end

    end
        

end

# particle entanglement spectrum
function plot_particle_entanglement_spectrum_paperplot()
    #using DelimitedFiles,PyCall

    cols = ["#C73E1D","#82AC9F","#36213E"]
    markers = "_"

    dataloc = get_folder_location("synth-dims/local-paperstuff/N=5_Nx=10_Ny=5_p=1_q=5/ParticleEntanglementSpectrum_NA=3")
    all_files = reverse(readdir(dataloc))

    fig,axs = subplots(1,2; figsize=(6,4))

    gs_counts = zeros(Int, length(all_files))
    
    for (ii,f) in enumerate(all_files)
        data = readdlm(joinpath(dataloc, f), '\t')

        ulr = split(split(f,"_")[end],"=")[end]

        idx   = data[:, 1]
        col2  = data[:, 2]
        col3  = data[:, 3]
        value = data[:, 4]

        gs_count = length(filter(x -> x < 6, value))

        xs = col3 .+ 5 .* col2
        ys = value

        gs_counts[ii] = gs_count

        axs[ii].scatter(xs, ys, c=cols[ii], marker=markers, s=60, label=L"U_{\mathrm{i}}"*"= $ulr")
        axs[ii].set_xlabel(L"K", fontsize=14)
        ii == 1 && (axs[ii].set_ylabel("PES", fontsize=14))
        axs[ii].set_ylim(1.0,15.0)
        axs[ii].legend(loc="upper right", fontsize=14)
    end

    xtickvals = collect(0:9)
    xticklabels = ["0","","2","","4","","6","","8",""]
    axs[1].set_xticks(xtickvals, xticklabels)
    axs[2].set_xticks(xtickvals, xticklabels)

    axs[1].tick_params(axis="both", which="major", labelsize=14)
    axs[2].tick_params(axis="both", which="major", labelsize=14)

    axs[1].tick_params(top=true, right=true)
    axs[1].tick_params(labeltop=false, labelright=false)

    axs[2].tick_params(top=true, right=true)
    axs[2].tick_params(labeltop=false, labelright=false)

    # label ulr = 0 ground state count
    axs[1].annotate(L"\mathcal{N}"*" = $(gs_counts[1])",
         xy=(4.5, 6.5), xycoords="data",
         fontsize=14,
         ha="center", va="center")

    # label ulr = 10 ground state count
    axs[2].annotate(L"\mathcal{N}"*" = $(gs_counts[2])",
         xy=(4.5, 4.3), xycoords="data",
         fontsize=14,
         ha="center", va="center")

    
    # ellipses which enclose the ground state manifold for counting
    patches = pyimport("matplotlib.patches")
    ellipse_u0 = patches.Ellipse(
        (4.5, 4.0),    # (x_center, y_center)
        12.0,           # width (x-extent)
        3.0;           # height (y-extent)
        angle = 0,    # rotation in degrees (counterclockwise)
        fill = false,
        linewidth = 2,
        edgecolor = cols[3],
        zorder = 10
    )
    ellipse_u10 = patches.Ellipse(
        (4.5, 3.0),    # (x_center, y_center)
        12.0,           # width (x-extent)
        1.0;           # height (y-extent)
        angle = 0,    # rotation in degrees (counterclockwise)
        fill = false,
        linewidth = 2,
        edgecolor = cols[3],
        zorder = 10
    )

    axs[1].add_patch(ellipse_u0)
    axs[2].add_patch(ellipse_u10)

    tight_layout()

end





########## Felix Palm ULR Length Plots ##########

# Felix phase diagram compare 16x8 to 8x4
#=function plot_finitesizescaling_ulr_phasetransition()
if false
    intstren = 4.0

    lx,ly,n = 12,8,6
    dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge/")#ulr-length")
    pdict = Dict([("Lx",lx),("Ly",ly),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc) 
    display(all_files)

    xis = range(0.0,8.0,length=11)
    flatnesses = zeros(Float64,length(xis))
    for f in all_files

        d,m = read_data(joinpath(dataloc,f); output_level=0)

        xi = m["corr_length"]

        xi_index = findfirst(x -> xis[x] == xi,1:11)

        if haskey(m,"fourpt_momentum")
            fourpt1 = m["fourpt_momentum"]

            #=fig = figure()
            fourpt1_includedsubset = zeros(Float64,lx,lx)
            for i in 3:lx-3
                for j in 1:length(diag(fourpt1,i))
                    fourpt1_includedsubset[j,j+i] = diag(fourpt1,i)[j]
                    fourpt1_includedsubset[j+i,j] = diag(fourpt1,-i)[j]
                end 
            end
            imshow(fourpt1_includedsubset,extent=(1,lx,1,lx),origin="lower",vmin=0.0)
            colorbar()
            title("Four-Point Correlator $(lx)x$(ly) N=$(n) ULR Length = $(xi)")=#

            subset_fourpt = vcat([diag(fourpt1,-i) for i in 3:lx-3]...,[diag(fourpt1,i) for i in 3:lx-3]...)
            flatness = minimum(subset_fourpt) / maximum(subset_fourpt)

            flatnesses[xi_index] = flatness
        end

    end

    #=dataloc_laughlin = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_laughlin16x8 = Dict([("layers",7),("particles",8),("hopping_anisotropy",1.0),("if_periodic_phys",true),("if_periodic_synth",true),("onsite_strength",0.0)])
    all_files_laughlin16x8 = find_data_file(pdict_laughlin16x8,"ttn",dataloc_laughlin)
    d,m = read_data(joinpath(dataloc_laughlin,all_files_laughlin16x8[1]); output_level=0)
    fourpt_laughlin16x8 = m["fourpt_momentum"]
    subset_fourpt_laughlin16x8 = vcat([diag(fourpt_laughlin16x8,i) for i in 3:lx-3]...)
    flatness_laughlin16x8 = minimum(subset_fourpt_laughlin16x8) / maximum(subset_fourpt_laughlin16x8)
    flatnesses[1] = flatness_laughlin16x8=#

    # normalize from laughlin tao-thouless
    #flatnesses ./= flatnesses[1]

    fig = figure()
    plot(xis,flatnesses ./ flatnesses[1],"-o",c="b",label="$(lx)x$(ly)")


    #=lx,ly,n = 8,8,4
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/ulr-length")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    
    xis = range(0.0,ly,length=11)
    #intstrens = range(0.0,4.0,length=11)
    flatnesses = ones(Float64,length(xis))
    gaps = ones(Float64,length(xis))
    splittings = ones(Float64,length(xis))
    for f in all_files

        d,m = read_data(joinpath(dataloc,f); output_level=0)

        xi = m["corr_length"]

        xi_index = findfirst(x -> xis[x] == xi,1:11)
        #intstren_index = findfirst(x -> intstrens[x] == m["U"][1],1:11)

        if !haskey(m,"fourpt_momentum_1")
            continue
        end

        fourpt1 = m["fourpt_momentum"]
        fourpt2 = m["fourpt_momentum_1"]

        fourpt_mixed = 0.5 * (fourpt1 .+ fourpt2)

        subset_fourpt = vcat([diag(fourpt_mixed,i) for i in 3:lx-3]...)
        flatness = minimum(subset_fourpt) / maximum(subset_fourpt)

        flatnesses[xi_index] = flatness

    end

    # normalize from laughlin tao-thouless
    flatnesses ./= flatnesses[1]

    plot(xis ./ ly, flatnesses,"-o",c="r",label="$(lx)x$(ly)")
    ylim([-0.05,1.1])=#

    legend()
    xlabel("Interaction Length")
    ylabel("Normalized Fourpt Flatness")
end=#






"fin"
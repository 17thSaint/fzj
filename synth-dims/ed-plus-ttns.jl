#####################################################
#=

This file contains data comparisons for both TTNs and ED

Depends on:
    synth-dims/long-range-ttn.jl
    review-practice-codes/observables.jl
    review-practice-codes/plottings.jl

    exact-diag/execute-ed.jl
    exact-diag/observables.jl
    exact-diag/plottings.jl

=#
######################################################

using JLD2
include("../other-funcs/include-other-files.jl")
include_other_files(["synth-dims/long-range-ttn.jl","review-practice-codes/observables.jl","other-funcs/basic-2d-observables.jl","exact-diag/execute-ed.jl","exact-diag/observables.jl"])
#include_other_files(["synth-dims/oneD-effective-LR.jl","synth-dims/plottings-oneD.jl"])
include_other_files(["other-funcs/basic-2d-plottings.jl","review-practice-codes/plottings.jl","exact-diag/plottings.jl"])

function plot_nrg_vs_intstren_fromdata_ttn(layers::Int64,which_strens::Union{String,Vector{Float64}}="all"; kwargs...)
    hanis = get(kwargs, :hopping_anisotropy, 1.0)
    if_gap = get(kwargs, :if_gap, true)
    levels_count = get(kwargs, :levels_count, 3)
    if_plot = get(kwargs, :if_plot, true)
    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    lx,ly = get_lattice_dims_from_layers(layers)
    n = get(kwargs, :particles, lx)
    pdict = Dict([("layers",layers),("particles",n),("hopping_anisotropy",hanis),("if_periodic_phys",true),("if_periodic_synth",true)])

    if_check = which_strens == "all" ? false : true

    all_files = find_data_file(pdict,"ttn",dataloc)

    nrgs::Dict{String,Vector{Float64}} = Dict([("1",Float64[])])
    for i in 2:levels_count
        nrgs[string(i)] = []
    end
    intstrens = []
    for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        if_done,all_checks = check_nrg_convergence(m,false)
        if if_check && !(m["onsite_strength"] in which_strens)
            continue
        end
        if !any(collect(values(all_checks)))
            continue
        end
        append!(intstrens,[m["onsite_strength"]])
        for i in 1:levels_count
            keyname = i == 1 ? "observer" : "observer_$(i-1)"
            if keyname in keys(m) && all_checks[string(i-1)]
                append!(nrgs[string(i)],[m[keyname].nrg[end]])
            else
                append!(nrgs[string(i)],[0.0])
            end
        end
    end

    if if_plot
        shift_nrg = if_gap ? nrgs["1"] : zeros(length(nrgs["1"]))
        cols = ["b","r","k"]
        fig = figure()
        for i in 1:levels_count
            scatter(intstrens,nrgs[string(i)] .- shift_nrg,c=cols[i],label="E$(i-1)")
        end
        xlabel("Interaction Strength")
        ylabel("Energy Difference")
        legend()
        if_gap && ylim([-0.1,0.7])

        title("Energy Spectrum for "*L"\rho_{1D}=1.0"*" TTNs $(lx)x$ly lattice")
    end

    return intstrens,nrgs
end

function plot_nrg_vs_intstren_fromdata_ed(lx::Int64,ly::Int64,which_strens::Union{String,Vector{Float64}}="all"; kwargs...)
    hanis = get(kwargs, :hopping_anisotropy, 1.0)
    if_gap = get(kwargs, :if_gap, true)
    levels_count = get(kwargs, :levels_count, 5)
    if_plot = get(kwargs, :if_plot, true)
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    n = get(kwargs, :particles, lx)
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",hanis),("if_periodic_x",true),("if_periodic_y",true)])

    if_check = which_strens == "all" ? false : true

    all_files = find_data_file(pdict,"ed",dataloc)

    display(all_files)

    nrgs::Dict{String,Vector{Float64}} = Dict([("1",Float64[])])
    for i in 2:levels_count
        nrgs[string(i)] = []
    end
    intstrens = []
    for f in all_files
        filename_dict = get_params_dict_from_filename(f)
        if haskey(filename_dict,"twist_angle1") || haskey(filename_dict,"twist_angle2")
            continue
        end
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        if if_check && !(m["U"][end] in which_strens)
            continue
        end
        append!(intstrens,[m["U"][end]])
        for i in 1:levels_count
            if length(d["nrg"]) < i
                append!(nrgs[string(i)],[0.0])
            else
                append!(nrgs[string(i)],[d["nrg"][i]])
            end
        end
    end

    if if_plot
        shift_nrg = if_gap ? nrgs["1"] : zeros(length(nrgs["1"]))
        cols = ["b","r","k","k","k"]
        fig = figure()
        for i in 1:levels_count
            scatter(intstrens,nrgs[string(i)] .- shift_nrg,c=cols[i],label="E$(i-1)")
        end
        xlabel("Interaction Strength")
        ylabel("Energy Difference")
        legend()
        if_gap && ylim([-0.1,0.7])
        xlim([-0.1,2.1])

        title("Energy Spectrum for "*L"\rho_{1D}=1.0"*" ED $(lx)x$ly lattice")
    end
    
    return intstrens,nrgs
end


#= look at finite size scaling of commensurate filling interaction strength spectrum
if false
    eight_xs, eight_ys = plot_nrg_vs_intstren_fromdata_ttn(6; particles=8, if_gap=true,if_plot=false)
    four_xs, four_ys = plot_nrg_vs_intstren_fromdata_ed(4,8; particles=4, if_gap=true,if_plot=false)
    three_xs, three_ys = plot_nrg_vs_intstren_fromdata_ed(3,8; particles=3, if_gap=true,if_plot=false)
    five_xs, five_ys = plot_nrg_vs_intstren_fromdata_ed(5,8; particles=5, if_gap=true,if_plot=false)

    eight_23gap = eight_ys["3"] .- eight_ys["2"]
    four_23gap = four_ys["3"] .- four_ys["2"]
    three_23gap = three_ys["3"] .- three_ys["2"]
    five_23gap = five_ys["3"] .- five_ys["2"]

    cols = get_cmap("plasma")
    divisor = 7

    fig = figure()
    plot(three_xs,three_23gap,"-p",c=cols(4/divisor),label="ED 3x8")
    plot(four_xs,four_23gap,"-p",c=cols(3/divisor),label="ED 4x8")
    plot(five_xs,five_23gap,"-p",c=cols(2/divisor),label="ED 5x8")
    plot(eight_xs,eight_23gap,"-p",c=cols(1/divisor),label="TTN 8x8")
    xlim([-0.1,2.1])
    ylim([-0.1,0.7])
    legend()
    xlabel("Interaction Strength")
    ylabel("E3 - E2")
end

# check the NRGs of synth-rectangle TTNs vs ED
if false
    ttn_intstrens,ttn_nrgs = plot_nrg_vs_intstren_fromdata_ttn(5; particles=4, if_plot=false)
    ed_intstrens,ed_nrgs = plot_nrg_vs_intstren_fromdata_ed(4,8; particles=4, if_plot=false)
    
    for level in 1:length(keys(ttn_nrgs))
        fig = figure()
        plot(ed_intstrens,ed_nrgs[string(level)],c="r",label="ED")
        scatter(ttn_intstrens,ttn_nrgs[string(level)],c="b",label="TTN")
        xlabel("Interaction Strength")
        ylabel("Energy Difference")
        legend()
        title("Energy Spectrum for "*L"\rho_{1D}=1.0"*" E$(level-1)")
    end
end

# build phase diagram for ULR vs rho1D at finite and infinite ULR from flatness
function build_phase_diagram_ulr_rho1d_flatness()
    configs = [(8,3,3),(8,4,4),(4,6,3),(3,8,3),(8,5,5),(8,7,7)]

    ulrs::Vector{Float64} = Float64[]
    flatnesses::Vector{Float64} = Float64[]
    oneDrhos::Vector{Float64} = Float64[]

    for (idx,config) in enumerate(configs)
        lx,ly,n = config
        append!(flatnesses,[twist_flatness_1deff(lx,ly,n; if_plot_spectrum=true,plot_title=" $(lx)x$ly N=$n")])
        append!(ulrs,[400.0])
        append!(oneDrhos,[n/lx])
    end

    #=for (lx,ly,n) in configs
        if (lx,ly,n) == (8,7,7)
            continue
        end
        local_strens,local_flats = twist_flatness_ed(lx,ly,n; if_plot=true, max_intstren=500.0)
        append!(ulrs,local_strens)
        append!(flatnesses,local_flats)
        append!(oneDrhos,ones(Float64,length(local_strens)) .* (n / lx))
    end=#

    bin_count = 100
    data_dict = bin_values(flatnesses,bin_count)
    bv = [data_dict[val] for val in flatnesses]
    min_nrgs2, max_nrgs2 = 0.0,1.0
    normalized_bv = [(val - minimum(bv)) / (maximum(bv) - minimum(bv)) * (max_nrgs2 - min_nrgs2) + min_nrgs2 for val in bv]

    #=fig = figure()
    scatter(oneDrhos, ulrs, c=normalized_bv, cmap="viridis")
    colorbar()
    xlabel(L"\rho_{1D}")
    ylabel("ULR")
    title("Flatness Phase Diagram")
    ylim([-0.1,4.1])=#

    plot_phasediag_ulrrho1d_flatness(oneDrhos,ulrs,bv,500.0)

end

# build_phase_diagram_ulr_rho1d for ft dd
if false
    configs = [(8,3,3),(4,6,3),(3,8,3),(8,4,4),(8,5,5)]

    ulrs::Vector{Float64} = Float64[]
    ftdds::Vector{Float64} = Float64[]
    oneDrhos::Vector{Float64} = Float64[]

    for (lx,ly,n) in configs
        local_strens,local_ftdd = findall_ft_dd(lx,ly,n; if_plot=false)
        append!(ulrs,local_strens)
        append!(ftdds,1 ./ local_ftdd)
        append!(oneDrhos,ones(Float64,length(local_strens)) .* (n / lx))
    end

    #=for (idx,config) in enumerate(configs)
        lx,ly,n = config
        append!(ftdds,[findall_ft_dd_1deff(lx,ly,n)])
        append!(ulrs,[900.0])
        append!(oneDrhos,[n/lx])
    end=#

    bin_count = 100
    data_dict = bin_values(ftdds,bin_count)
    bv = [data_dict[val] for val in ftdds]
    min_nrgs2, max_nrgs2 = minimum(ftdds),maximum(ftdds)
    normalized_bv = [(val - minimum(bv)) / (maximum(bv) - minimum(bv)) * (max_nrgs2 - min_nrgs2) + min_nrgs2 for val in bv]

    #=fig = figure()
    scatter(oneDrhos, ulrs, c=normalized_bv, cmap="plasma")
    colorbar()
    xlabel(L"\rho_{1D}")
    ylabel("ULR")
    title("FT-DD Ratio Phase Diagram")
    yscale("log")=#

    plot_phasediag_ulrrho1d_flatness(oneDrhos,ulrs,normalized_bv,1000.0)
    title("Inverse FT DD Max Phase Diagram")
end

# check perturbative 1Deff against ED # it seems this only works for rho1D = 1/2 and 1.0
if false
    lx,ly,n = 8,3,3

    intstrens = vcat(range(1000.0,10000.0,length=21),range(0.0,100.0,length=21))
    cols = get_colors(3)
    for (idx,intstren) in enumerate(intstrens)
        if intstren > 100.0
            pdict_mps = Dict([("Lphys",lx),("Lsynth",ly),("particles",n),("interaction_strength",intstren),("if_remapping",false),("es_count",2),("nrgtol",1e-6),("mdim",200),("if_periodic_phys",true),("if_periodic_synth",true),("filling",0.5),("if_find_data",false),("if_save_data",false)])
            psis,rhos,nrgs,mparas,if_found = run_normal_1deffmps(pdict_mps)
            result = ft_densitydensity_correlation(pi/2,psis[1])
        else
            pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
            psis,nrgs,rhos,filepath_ed,if_exists_ed,latparas,hamparas = run_normal_ed(pdict_ed)
            result = ft_densitydensity_correlation(pi/2,psis[1],latparas)
        end

        scatter(intstren,result,c="b")
        xscale("log")
        xlabel("Interaction Strength")
        ylabel("FT-DD pi/2")

    end
end

# checking ED vs full infinite limit 1Deff # need larger intstren for (8,5,5)
if false
    cols = get_colors(3)
    configs = [(8,3,3),(8,4,4),(4,6,3),(3,8,3),(8,5,5)]
    for (lx,ly,n) in configs
        #lx,ly,n = 8,3,3
        pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
        all_files_ed = find_data_file(pdict_ed,"ed",get_folder_location("cluster-data/exact-diag/torus"))
        filter!(x -> !occursin("twist",x),all_files_ed)
        fig = figure()
        intstrens_ed = Float64[]
        for f in all_files_ed
            d,m = read_data_jld2(get_folder_location("cluster-data/exact-diag/torus") * "/" * f; output_level=0)
            nrgs_ed = d["nrg"]
            append!(intstrens_ed,[m["U"][end]])
            for i in 1:3
                scatter(m["U"][end],nrgs_ed[i],c=cols[i])
            end
        end

        pdict_mps = Dict([("Lphys",lx),("Lsynth",ly),("nbosons",n),("if_periodic_phys",true),("if_periodic_synth",true)])
        all_files = find_data_file(pdict_mps,"mps",get_folder_location("cluster-data/synth-dims/excited-states"))
        d,m = read_data_jld2(get_folder_location("cluster-data/synth-dims/excited-states") * "/" * all_files[1]; output_level=0)
        nrgs_mps = zeros(Float64,3)
        for i in 1:3
            obs_string = i == 1 ? "observer" : "observer_$(i-1)"
            nrgs_mps[i] = m[obs_string].energies[end]
        end
        sort!(nrgs_mps)
        for i in 1:3
            plot([0.0,maximum(intstrens_ed)],[nrgs_mps[i],nrgs_mps[i]],c=cols[i],label="E$i")
        end
        legend()

        xlabel("Interaction Strength")
        ylabel("Energy")
        title("ED (scatter) vs 1Deff (line) for $(lx)x$ly N=$n")
    end

    
    
end

# energy spectrum scaling with Lx/Ly = 2.0
if false
    configs = [(6,3,3),(8,4,4),(10,5,5)]
    for config in configs
        lx,ly,n = config
        if lx == 10
            filename = "../cluster-data/exact-diag/torus/ed-Ly-5-interaction_strength-0.0-Lx-10-if_pinning-true-N-5-alpha-0.2-if_periodic_x-true-if_periodic_y-true-hopping_anisotropy-1.0.jld2"
            d,m = read_data_jld2(filename; output_level=0)
            nrgs = d["nrg"]
        else
            pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",0.0),("filling",0.5),("lr","all"),("if_find_data",true),("if_save_data",false),("nev",10)])
            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)
        end

        scatter(lx,nrgs[1] - nrgs[1],c="b")
        scatter(lx,nrgs[2] - nrgs[1],c="g")
        scatter(lx,nrgs[3] - nrgs[1],c="r")

        xlabel("Lx")
        ylabel("Energy Gap")
        title("Energy Spectrum Finite Size Scaling")
    end

    ttn_configs = [(16,8,8)]
    for config in ttn_configs
        lx,ly,n = config
        layers = get_layers_from_latticesize(lx,ly)
        pdict = Dict([("layers",layers),("particles",n),("hopping_anisotropy",1.0),("if_periodic_phys",true),("if_periodic_synth",true),("onsite_strength",0.0)])
        dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
        all_files = find_data_file(pdict,"ttn",dataloc)
        d,m = read_data_jld2(dataloc * "/" * all_files[1]; output_level=0)
        nrgs::Dict{String,Float64} = Dict([("1",m["observer"].nrg[end])])

        for i in 2:3
            keyname = "observer_$(i-1)"
            if keyname in keys(m)
                nrgs[string(i)] = m[keyname].nrg[end]
            else
                nrgs[string(i)] = 0.0
            end
        end

        scatter(lx,nrgs["1"] - nrgs["1"],c="b")
        scatter(lx,nrgs["2"] - nrgs["1"],c="g")
        scatter(lx,nrgs["3"] - nrgs["1"],c="r")

        xlabel("Lx")
        ylabel("Energy Gap")
        title("Energy Spectrum Finite Size Scaling")
    end
end

# check the FT-D scaling with Lx/Ly = 2.0
if false
    limit_vals = Float64[]
    limit_sizes = Int64[]

    # ED section
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    ks = range(0,2*pi,length=50)
    configs = [(10,5,5),(8,4,4)]#(6,4,3),(8,4,4),(10,4,5)]
    for (lx,ly,n) in configs
        if_pinning = true
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
        all_files = find_data_file(pdict,"ed",dataloc; output_level=0)

        filter!(x -> !occursin("twist_angle1",x),all_files)
        filter!(x -> !occursin("mk",x),all_files)
        filter!(x -> occursin("if_pinning",x),all_files)

        display(all_files)

        cdw_sfs = Float64[]
        intstrens = Float64[]
        if_plot = false
        for f in all_files

            filename_dict = get_params_dict_from_filename(f)
            if lx == 10 && filename_dict["interaction_strength"] == 0.8
                continue
            end
            d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

            if lx == 8 && round(m["U"][end],digits=0) != m["U"][end]
                continue
            end

            append!(intstrens,m["U"][end])

            latparas = get_lattice_params_from_metadata(m)
            occs = get_occupancy(d["state"][1],latparas; if_plot=if_plot,plot_title="$(lx)x$(ly) n=$n ULR=$(m["U"][end])")
            ftd = ft_density([pi,0],occs)
            append!(cdw_sfs,[abs(ftd)])

            if m["U"][end] == 1000.0
                push!(limit_vals,abs(ftd))
                push!(limit_sizes,lx)
            end

        end

        #fig = figure()
        scatter(intstrens,cdw_sfs,label="$(lx)x$(ly)")
        xlabel("Interaction Strength")
        ylabel("CDW SF at k=(pi,0)")
        title("CDW SF at "*L"\rho_{1D}=1/2")
        xscale("log")
        legend()
    end

    # TTN section
    dataloc = get_folder_location("cluster-data/synth-dims/torus")
    configs = [(7,8),(7,6)]
    for (layers,particles) in configs
        pdict = Dict([("layers",layers),("particles",particles),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
        all_files = find_data_file(pdict,"ttn",dataloc)
        display(all_files)

        intstrens = Float64[]
        ftvals = ComplexF64[]
        for f in all_files
            d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

            #= checking convergence of energies
            all_convs = [false,false,false]
            for i in 1:3
                observerkey = i == 1 ? "observer" : "observer_$(i-1)"
                if observerkey in keys(m)
                    all_convs[i] = abs(m[observerkey].nrg[end] - m[observerkey].nrg[end-1]) < m["nrgtol"]
                end
            end
            println("For $(m["onsite_strength"]) all convs are ",all_convs)=#

            lx,ly = 2*particles,particles
            occs = get_occupancy(d["densmat"]; plot_title="Intstren = $(m["onsite_strength"])",if_plot=false)[1:lx,1:ly]
            ftval = ft_density([pi,0.0],occs)
            push!(intstrens,m["onsite_strength"])
            push!(ftvals,ftval)

            if m["onsite_strength"] == 1000.0
                push!(limit_vals,abs(ftval))
                push!(limit_sizes,particles*2)
            end
        end
        scatter(intstrens,abs.(ftvals),label="$(2*particles)x$particles")
        legend()
    end


    fig = figure()
    scatter(limit_sizes,limit_vals)
    xlabel("Lx")
    ylabel("CDW SF at k=(pi,0)")
    title("CDW SF at "*L"\rho_{1D}=1/2"*", ULR=1000.0")
    yscale("log")

end

# finite size scaling of FT-DD at ULR = 0.0
if false
    intstren = 100.0

    # ED version with 8x4 and 10x5
    configs = [(8,4,4),(10,5,5)]
    for config in configs
        lx,ly,n = config
        dataloc = get_folder_location("cluster-data/exact-diag/torus")
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
        all_files = find_data_file(pdict,"ed",dataloc; output_level=0)
        filter!(x -> !occursin("twist_angle1",x),all_files)
        filter!(x -> !occursin("mk",x),all_files)
        filter!(x -> !occursin("pinning",x),all_files)
        f = all_files[1]
        d,m = read_data_jld2(dataloc * "/" * f; output_level=1)
        dds = m["densitydensity"]
        ft_dds_stripe = abs(ft_densitydensity([pi,0],dds))
        scatter(lx,ft_dds_stripe,c="r")
    end

    # TTN version with 16x8 and 12x6
    configs = [8,6]
    for n in configs
        dataloc = get_folder_location("cluster-data/synth-dims/torus")
        pdict = Dict([("layers",7),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
        all_files = find_data_file(pdict,"ttn",dataloc)
        f = all_files[1]
        d,m = read_data_jld2(dataloc * "/" * f; output_level=1)
        dds = m["densitydensity"]
        ft_dds_stripe = abs(ft_densitydensity([pi,0],dds))
        scatter(n*2,ft_dds_stripe,c="b")
    end

    xlabel("Lx")
    ylabel("FT-DD at k=(" * L"\pi" * ",0)")
    title("Finite Size Scaling of FT-DD at ULR=$intstren")
    yscale("log")
end

# testing 2pt MPO FT constructor against known density matrix FT
if false
    lx,ly,n = 4,4,2
    params_dict = Dict([("hopping_anisotropy",1.0),("if_check_fluxes",false),("es_count",0),("expander_fraction",0.5),("particles",n),("layers",Int(log(2,lx*ly))),("mdim",100),("if_save_data",false),("alpha",0.0),("onsite_strength",0.0),("lr",0),("if_periodic_phys",true),("if_periodic_synth",true)])
    #psi, hamilthere, obs, rho, rt = run_synth_dims_generic(params_dict)

    ks = range(0,1.0,length=16)
    vals = zeros(Float64,length(ks),length(ks))
    for (idx,kx) in enumerate(ks)
        for (idx2,ky) in enumerate(ks)
            println("Working on kx = $kx, ky = $ky")
            rho = two_point_mpo(psi; momentum = [kx,ky])
            vals[idx2,idx] = real(calculate_mpo_expectation(psi,rho))
        end
    end
    fig = figure()
    imshow(vals,origin="lower",extent=[0,1,0,1])
    colorbar()
    xlabel("kx")
    ylabel("ky")

    rho_real = density_matrix(psi)
    dm_vals = zeros(Float64,length(ks),length(ks))
    for (idx,kx) in enumerate(ks)
        for (idx2,ky) in enumerate(ks)
            println("Working on kx = $kx, ky = $ky")
            dm_vals[idx2,idx] = real(ft_density_matrix(rho_real,[kx,ky],lx,ly))
        end
    end
    fig = figure()
    imshow(dm_vals,origin="lower",extent=[0,1,0,1])
    colorbar()
    xlabel("kx")
    ylabel("ky")

end

# test zigzag curve
if false
    lx,ly = 8,4
    zc = zigzag_curve(lx,ly)
    display(zc)
end=#

# do 4pt momentum MPO
if false
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))
    intstren = 0.0

    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    pdict = Dict([("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    #pdict = Dict([("particles",n),("expander_fraction",100),("flux_direction","synth"),("layers",layers),("mdim",200),("if_save_data",false),("if_find_data",false),("filling",0.5),("onsite_strength",intstren),("lr",0),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    filter!(x -> !occursin("if_synth_rectangle",x),all_files)
    display(all_files)

    #for f in all_files
        f = all_files[1]
        d,m = read_data(dataloc * "/" * f; output_level=0)
        #!haskey(d,"ttn") && continue
        #haskey(m,"fourpt_momentum") && continue
        psi = d["ttn"]
        #TTN.move_ortho!(psi,(layers,1))
        #psi2 = d["ttn_1"]
        #get_occupancy(psi; if_plot=true,plot_title="$(m["onsite_strength"])")
        #all_results = run_synth_dims_generic(pdict)
        #psi = all_results[1]
        lat = TTN.physical_lattice(psi.net)

        mapss = zigzag_curve(lx,ly)

        ks = [n/ly for n in 1:lx]
        mp = [0.0,4/ly]
        twopt_vals = zeros(Float64,length(ks))
        real_twopt_vals = zeros(Float64,length(ks))
        for (idx,ky) in enumerate(ks)
            twopt_vals[idx] = two_point(psi,mp,[0.0,ky])
            scatter(ky *ly,twopt_vals[idx],c="b")
            #real_twopt_vals[idx] = abs(two_point_real(psi, mp, [0.0,ky]))
            #scatter(ky *ly,real_twopt_vals[idx],c="g")
        end

        #display(vals)
        #fig = figure()
        #scatter(ks .* ly,fourpt_vals,c="b",label="Real")
        xlabel("Momentum k = m / Ly, m' = $(Int(mp[2]*ly))")
        ylabel("Two Point Momentum")
        title("Two Point Momentum for $(lx)x$(ly) N=$n")


    #end
end

# do 2pt momentum MPO
if false
    lx,ly,n = 4,4,2
    layers = Int(log(2,lx*ly))
    intstren = 0.0

    
    #=dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    pdict = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    filter!(x -> !occursin("if_synth_rectangle",x),all_files)
    display(all_files)
    f = all_files[1]
    d,m = read_data(dataloc * "/" * f; output_level=0)
    psi = d["ttn"]=#
    
    pdict = Dict([("particles",n),("es_count",1),("if_pinning",true),("if_check_fluxes",false),("max_occ",1),("expander_fraction",100),("flux_direction","synth"),("layers",layers),("mdim",200),("if_save_data",false),("if_find_data",false),("filling",0.5),("onsite_strength",1000000.0),("lr",0),("if_periodic_phys",true),("if_periodic_synth",true)])
    #all_results = run_synth_dims_generic(pdict)
    psi1 = all_results[1][1]
    psi2 = all_results[1][2]

    ks = [2/ly]#[n/ly for n in 1:lx]
    mp = [0.0,2/ly]

    twopt_vals = two_point([psi1,psi2],mp,mp)

    #=twopt_vals = zeros(Float64,length(ks))
    real_twopt_vals = zeros(Float64,length(ks))
    for (idx,ky) in enumerate(ks)
        #fourpt_vals[idx] = four_point(psi,mp,[0.0,ky])
        real_twopt_vals[idx] = abs(two_point_real(psi, mp, [0.0,ky]))
    end

    #scatter(ks .* ly,fourpt_vals,c="b",label="MPO")
    scatter(ks .* ly,real_twopt_vals,c="g",label="Real")
    legend()
    xlabel("Momentum k = m / Ly, m' = $(Int(mp[2]*ly))")
    ylabel("Two Point Momentum")
    title("Two Point Momentum for $(lx)x$(ly) N=$n")=#

    rho_ttn1 = all_results[end-1][1]
    rho_ttn2 = all_results[end-1][2]

end

# plot 4pt momentum
if false
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))
    #intstren = 100.0

    dataloc = [get_folder_location("cluster-data/synth-dims/excited-states"),get_folder_location("cluster-data/synth-dims/torus")]
    pdict = Dict([("layers",layers),("particles",n),("if_periodic_phys",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    filter!(x -> !occursin("if_synth_rectangle",x),all_files)
    display(all_files)

    #f = all_files[1]
    ulrs = [0.0,0.5,1.0,10.0,300.0]
    fig1 = figure()
    for f in all_files
        d,m = read_data(f; output_level=0)

        #get_occupancy(d["densmat"]; plot_title="ULR=$(m["onsite_strength"])")

        #!haskey(m,"twopt_momentum") && continue
        !haskey(m,"fourpt_momentum") && continue
        !haskey(d,"densmat") && continue

        ulr = get_params_dict_from_filename(split(f,"/")[end])["onsite_strength"]
        !(ulr in ulrs) && continue

        #twopt_vals = m["twopt_momentum"]

        #=fig = figure()
        imshow(twopt_vals,origin="lower",extent=[0,lx,0,lx],vmin=0.0)
        colorbar()
        xlabel("kx")
        ylabel("k2")
        title("2-Point Momentum at ULR=$(m["onsite_strength"])")=#
        #twopt_vals = two_point_densmat(d["densmat"],lx,ly)
        #display(twopt_vals)
        #scatter(m["onsite_strength"],mean(twopt_vals),c="b")
        fourpt_vals = m["fourpt_momentum"]
        #scatter(m["onsite_strength"],fourpt_vals[5,5],c="b")
        #xlabel("Interaction Strength")
        #ylabel("Fourpt at (4,4)")

        #normalized_4pt = normalize_four_point(fourpt_vals,twopt_vals)

        plot(0:size(fourpt_vals,1)-1,fourpt_vals[:,Int(lx/2)+1],"-p",label="$ulr")
        title(L"$\langle \hat{a}_{k}^{\dagger} \hat{a}_{k'}^{\dagger} \hat{a}_{k'} \hat{a}_k \rangle / \langle \hat{n}_{k} \rangle \langle \hat{n}_{k'} \rangle$"*" for $(lx)x$(ly) N=$n Range ULR")
        xlabel("Momentum k = n / Lx, n' = $(Int(lx/2))")
        ylabel("Four Point Momentum")
        legend()

        #=plot(0:size(fourpt_vals,1)-1,diag(normalized_4pt),"-p",label="$ulr")
        title(L"$\langle \hat{a}_{k}^{\dagger} \hat{a}_{k}^{\dagger} \hat{a}_{k} \hat{a}_k \rangle / \langle \hat{n}_{k} \rangle \langle \hat{n}_{k} \rangle$"*" for $(lx)x$(ly) N=$n Range ULR")
        xlabel("Momentum k = n / Lx")
        ylabel("Four Point Momentum")
        legend()=#

        #=scatter(m["onsite_strength"],sum(fourpt_vals),c="b")
        xlabel("Interaction Strength")
        ylabel("Sum of 4pt")
        title("Sum of 4pt on $(lx)x$ly N=$n")=#

        #plot(0:lx,fourpt_vals[:,5] ./ (lx*ly)^4,"-p",label="$(m["onsite_strength"])")
        #xlabel("Momentum")
        #ylabel("Four Point Momentum")

        #=fig = figure()
        imshow(fourpt_vals,origin="lower",extent=[0,lx,0,lx],vmin=0.0)
        colorbar()
        xlabel("kx")
        ylabel("k2")
        title("4-Point Momentum at ULR=$(m["onsite_strength"])")=#
    end
        

end

# testing 4pt momentum on known density wave
if false
    lx,ly,n = 4,4,2
    layers = Int(log(2,lx*ly))

    stren = 0.0
    params_dict = Dict([("hopping_anisotropy",1.0),("es_count",1),("cutoff",1e-10),("particles",n),("layers",layers),("mdim",200),("expander_fraction",100),("if_save_data",false),("filling",0.5),("onsite_strength",stren),("lr","all"),("if_periodic_phys",true),("if_periodic_synth",true)])
    psi, hamilthere, obs, rho, rt = run_synth_dims_generic(params_dict)

    #lat = TTN.physical_lattice(psi[1].net)
    #mapss = zigzag_curve(lx,ly)

    fig = figure()
    ks = [n/lx for n in 0:lx]
    fourpt_vals1 = zeros(Float64,length(ks))
    fourpt_vals2 = zeros(Float64,length(ks))
    for (idx,kx) in enumerate(ks)
        println("Working on kx = $kx")

        fourpt_val = real.(four_point([psi[1],psi[2]], [kx,0], [kx,0]))
        fourpt_vals1[idx] = fourpt_val[1]
        fourpt_vals2[idx] = fourpt_val[2]

        #twopt_val = real.(two_point([psi[1],psi[2]], [kx,0], [kx,0]))
        #twopt_vals1[idx] = twopt_val[1]
        #twopt_vals2[idx] = twopt_val[2]

        if idx == 1
            scatter(kx*lx,fourpt_val[1],c="b",label="TTN MPO 1")
            scatter(kx*lx,fourpt_val[2],c="g",label="TTN MPO 2")
        else
            scatter(kx*lx,fourpt_val[1],c="b")
            scatter(kx*lx,fourpt_val[2],c="g")
        end
        xlabel("Momentum k, k' = k")
        ylabel("Four Point Momentum")
    end

end

# do multi-state 4pt momentum MPO
if false
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))
    intstren = 0.0

    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    pdict = Dict([("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    filter!(x -> !occursin("if_synth_rectangle",x),all_files)
    display(all_files)

    f = all_files[1]

    d,m = read_data(joinpath(dataloc,f); output_level=0)

    @assert haskey(d,"ttn")
    @assert haskey(d,"ttn_1")

    psi1 = d["ttn"]
    psi2 = d["ttn_1"]

    lat = TTN.physical_lattice(psi1.net)
    mapss = zigzag_curve(lx,ly)
    ks = [n/lx for n in 0:lx]
    for (idx,kx) in enumerate(ks)
        fourpt = four_point_mpo(psi1; momentum1 = [kx,0], momentum2 = [0.5,0], mapping = mapss)
        fourpt_wrapped = easy_mpowrapper(fourpt, lat; mapping=mapss)
        mat::Matrix{Float64} = zeros(Float64,2,2)
        mat[1,1] = real(calculate_mpo_expectation(psi1, psi1, fourpt_wrapped)) / (lx*ly)^2
        mat[1,2] = real(calculate_mpo_expectation(psi1, psi2, fourpt_wrapped)) / (lx*ly)^2
        mat[2,1] = real(calculate_mpo_expectation(psi2, psi1, fourpt_wrapped)) / (lx*ly)^2
        mat[2,2] = real(calculate_mpo_expectation(psi2, psi2, fourpt_wrapped)) / (lx*ly)^2

        if idx == 1
            scatter(idx-1,mat[1,1],c="r",label="11")
            scatter(idx-1,mat[1,2],c="g",label="12")
            scatter(idx-1,mat[2,1],c="b",label="21")
            scatter(idx-1,mat[2,2],c="k",label="22")
            title(L"$\langle \hat{a}_{k}^{\dagger} \hat{a}_{k'}^{\dagger} \hat{a}_{k'} \hat{a}_k \rangle$"*" for $(lx)x$(ly) N=$n ULR=$(m["onsite_strength"])")
            xlabel("Momentum k = n / Lx, n' = $(Int(lx/2))")
            ylabel("Four Point Momentum")
            legend()
        else
            scatter(idx-1,mat[1,1],c="r")
            scatter(idx-1,mat[1,2],c="g")
            scatter(idx-1,mat[2,1],c="b")
            scatter(idx-1,mat[2,2],c="k")
        end
    end

end

# test 4pt momentum with ED
if false
    lx,ly,n = 8,4,4
    intstren = 300.0
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("flux_direction","y"),("if_check_fluxes",false),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",false)])
    states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)
    #=lattice_params,hamilt_params,running_args = get_normal_model_params_ed(pdict)
    full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=running_args.basis_dataloc)   
    lattice_params["full_basis"] = full_basis=#


    #fig = figure()
    ks = [n/ly for n in 1:lx]
    ed_fourpt_vals = zeros(Float64,length(ks))
    mp = [0.0,4/ly]
    for (idx,ky) in enumerate(ks)
        ed_fourpt_vals[idx] = abs(ft_fourpt_alberto(states[1],mp,[0.0,ky],lattice_params))
        scatter(ky*ly,ed_fourpt_vals[idx],c="b")
    end
    xlabel("Momentum k = n / Ly, m' = $(Int(mp[2]*ly))")
    ylabel("Four Point Momentum")
    title("ED 4pt Momentum $(lx)x$(ly) N=$n ULR=$intstren")
end

# test 2pt momentum with ED
if false
    println("Starting ED")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_pinning",true),("flux_direction","y"),("if_check_fluxes",false),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",false)])
    states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)

    #=ed_twopt_vals = zeros(Float64,length(ks))
    for (idx,ky) in enumerate(ks)
        ed_twopt_vals[idx] = abs(ft_twopt_alberto(states[1],mp,[0.0,ky],lattice_params))
    end
    scatter(ks .* ly,ed_twopt_vals,c="r",label="ED")
    legend()=#

    #ed_twopt_vals = ft_twopt_alberto(states[1:2],mp,mp,lattice_params)

    #rho_ed1 = density_matrix(states[1],lattice_params)
    #rho_ed2 = density_matrix(states[2],lattice_params)
end

# compare 4pt momentum max_occs MPO
if false
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))

    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    pdict = Dict([("layers",layers),("particles",n),("onsite_strength",0.0),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    display(all_files)

    ho_dataloc = get_folder_location("cluster-data/synth-dims/torus")
    higherocc_filename = "ttn-if_periodic_phys-true-onsite_strength-1000.0-lr-0-particles-4-if_gpu-true-max_occ-2-alpha-0.25-if_periodic_synth-true-layers-5-hopping_anisotropy-1.0.h5"
    working_files = [joinpath(dataloc,all_files[1]),joinpath(ho_dataloc,higherocc_filename)]

    ks = [n/ly for n in 1:lx]
    mp = [0.0,4/ly]
    col = "b"
    for f in working_files
        d,m = read_data(f)
        psi = d["ttn"]
        fourpt_vals = zeros(Float64,length(ks))
        for (idx,ky) in enumerate(ks)
            fourpt_vals[idx] = four_point(psi,mp,[0.0,ky])
        end
        scatter(ks .* ly,fourpt_vals,c="b",label="MaxOcc=$(m["max_occ"])")
        xlabel("Momentum k = n / Lx, m' = $(Int(mp[2]*ly))")
        ylabel("Four Point Momentum")
        title("TTN 4pt Momentum $(lx)x$(ly) N=$n range Maximum Occupancy")
        global col = "r"
    end
end

# compare ED TTN ground states
if false
    s1c = (2,1)
    s2c = (2,1)
    s1l = linear_index(s1c,lx,ly)
    s2l = linear_index(s2c,lx,ly)

    ed_hoppingcorr = abs(hopping_probability(states[1],s1c,s2c,lattice_params))
    println("ED Hopping Probability = ",ed_hoppingcorr)

    ttn_hoppingcorr = abs(TTN.correlation(psi,"Adag","A",s1c,s2c))
    println("TTN Hopping Probability = ",ttn_hoppingcorr)

    s1_revlin = lx*ly - s1l + 1
    s2_revlin = lx*ly - s2l + 1
    s1_revcorr = coordinate(s1_revlin,lx,ly)
    s2_revcorr = coordinate(s2_revlin,lx,ly)
    ttn_revhopcorr = abs(TTN.correlation(psi,"Adag","A",s1_revcorr,s2_revcorr))
    println("TTN Reverse Hopping Probability = ",ttn_revhopcorr)
end

# test sum rule using ED
if false
    lx,ly,n = 4,4,2
    intstren = 0.0
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_pinning",true),("flux_direction","y"),("if_check_fluxes",false),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",false)])
    states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)

    sumval = 0.0*im
    ks = [n/lx for n in 1:lx]
    for (idx,kx) in enumerate(ks)
        for (idx2,ky) in enumerate(ks)
            global sumval += ftfull_twopt(states[1],[kx,ky],[kx,ky],lattice_params)
        end
    end

    println("Final Sum value is ",sumval)
end

# test sum rule using TTNs
if true
    lx,ly,n = 4,4,2
    layers = Int(log(2,lx*ly))
    intstren = 0.0
    pdict = Dict([("particles",n),("if_check_fluxes",false),("max_occ",1),("expander_fraction",100),("flux_direction","synth"),("layers",layers),("mdim",200),("if_save_data",false),("if_find_data",false),("filling",0.5),("onsite_strength",1000000.0),("lr",0),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_results = run_synth_dims_generic(pdict)

    sumval = 0.0*im
    ks = [n/lx for n in 1:lx]
    for (idx,kx) in enumerate(ks)
        for (idx2,ky) in enumerate(ks)
            global sumval += two_point(all_results[1],[kx,ky],[kx,ky])
        end
    end

    println("Final Sum value is ",sumval)
end




























"fin"
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

include("../other-funcs/include-other-files.jl")
include_other_files(["synth-dims/long-range-ttn.jl","review-practice-codes/observables.jl","exact-diag/execute-ed.jl","exact-diag/observables.jl"])
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


# look at finite size scaling of commensurate filling interaction strength spectrum
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

# testing MPO construction of 4point momentum correlator
if true
    lx,ly,n = 2,2,2
    params_dict = Dict([("hopping_anisotropy",1.0),("if_check_fluxes",false),("es_count",0),("expander_fraction",0.5),("particles",n),("layers",Int(log(2,lx*ly))),("mdim",100),("if_save_data",false),("alpha",0.0),("onsite_strength",0.0),("lr",0),("if_periodic_phys",true),("if_periodic_synth",true)])
    psi, hamilthere, obs, rho, rt = run_synth_dims_generic(params_dict)



    #creat = projected_op_mpo(psi,"Adag"; if_wrap=false)#,what_given_inds=["pull","pull"])
    #anh = projected_op_mpo(psi, "A"; if_wrap = false)
    
    creat = single_point_mpo(psi, "Adag"; if_wrap=false)
    annih = single_point_mpo(psi, "A"; if_wrap=false)
    

    #rho = two_point_mpo_reverse(psi; if_wrap=true)
    #val = calculate_mpo_expectation(psi,creat)
    #println("The value is $val")
end































"fin"
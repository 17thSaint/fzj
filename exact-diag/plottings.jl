#####################################################
#=

This file contains plotting functions for ED results

Depends on:
    other-funcs/basic-2d-stuff.jl
    two-dimensions.jl
    observables.jl
    hatsugai-mbcn.jl

=#
######################################################

using LaTeXStrings,PyPlot

function get_colors(nev::Int)
    cols = ["b","g","r","m","c"]
    if nev > length(cols)
        cols = repeat(cols,ceil(Int,nev/length(cols)))
    end
    return cols
end

# 2D phase diagram of hopping anisotropy vs interaction strength and energy gap
# needs testing, not sure about the commented out stuff
function nrg_gap_phase_diagram(lx::Int64=4,ly::Int64=8,N::Int64=4)

    cols = ["b","r","g","m","c"]
    if 10 > length(cols)
        cols = repeat(cols,ceil(Int,10/length(cols)))
    end

    nrgs1 = []
    nrgs2 = []
    nrgs3 = []
    nrgs4 = []
    xs = []
    ys = []
    
    params_dict = Dict([("Lx",lx),("Ly",ly),("N",N),("if_periodic_x",true),("if_periodic_y",true)])
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)

        for f in files
            data,metadata = read_data_jld2(f,dataloc; output_level=0)
            nrgs = data["nrg"]
            intstren = metadata["U"][end]
            anis = metadata["hopping_anisotropy"]

            if anis > 1.3
                continue
            end

            append!(xs,[anis])
            append!(ys,[intstren])
            #append!(nrgs1,[nrgs[2] - nrgs[1]])
            append!(nrgs2,[nrgs[3] - nrgs[2]])
            #append!(nrgs3,[nrgs[4] - nrgs[1]])
            #append!(nrgs4,[nrgs[5] - nrgs[1]])

        end
    
    bin_count = 100
    data_dict = bin_values(nrgs2,bin_count)
    bv = [data_dict[val] for val in nrgs2]
    min_nrgs2, max_nrgs2 = minimum(nrgs2), maximum(nrgs2)
    normalized_bv = [(val - minimum(bv)) / (maximum(bv) - minimum(bv)) * (max_nrgs2 - min_nrgs2) + min_nrgs2 for val in bv]

    #= now try to find linear fit to maximum of nrgs2 at each hopping anis row
    hanises = unique(xs)
    #max_intstrens = zeros(Float64,length(hanises))
    min_intstrens = zeros(Float64,length(hanises))
    for (j,hanis) in enumerate(hanises)
        indices = findall(x -> x == hanis, xs)
        if length(indices) == 1
            #println("Only one data point for Anis=$hanis")
            continue
        end
        #max_index = findfirst(x -> x == maximum(nrgs2[i] for i in indices), nrgs2)
        #max_intstrens[j] = ys[max_index]
        relevant_nrgs = [nrgs2[i] for i in indices]
        min_index = findfirst(i -> isapprox(relevant_nrgs[i],0.0,atol=1e-6),1:length(relevant_nrgs))
        println("Found Level Crossing for Anis=$hanis at index $min_index")
        println("It has value $(relevant_nrgs[min_index])")
        println("This is intstrength of $(ys[indices[min_index]]) \n")
        min_intstrens[j] = ys[indices[min_index]]
    end

    linfit(x,p) = p[1] .* x .+ p[2]
    fit = curve_fit(linfit,min_intstrens,hanises,[1.0,0.0])
    fit_xs = range(minimum(ys),2.0,length=10)
    plot(fit_xs,linfit(fit_xs,fit.param),c="r",label="m=$(round(fit.param[1],digits=3))")=#

    fig = figure()
    scatter(ys, xs, c=normalized_bv, cmap="viridis")
    colorbar()
    
    ylim(minimum(xs)-0.05,0.05+maximum(xs))

    ylabel("Hopping Anisotropy")
    xlabel("Interaction Strength")
    legend()
    title("Energy Gap btw 2nd and 3rd 4x8 N=4")

end

# plots time scaling with Hilbert space dimension for ED using functional form instead of full matrix
function timescaling_functional_ed()
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    params_dict = Dict([("if_periodic_x",true)])
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    for f in files
        data,metadata = read_data_jld2(f,dataloc)
        if_func = metadata["if_function"]
        lattice_params = get_lattice_params_from_metadata(metadata)
        dimHilb = binomial(lattice_params["Lx"]*lattice_params["Ly"],lattice_params["N"])

        col = if_func ? "r" : "b"
        scatter(dimHilb,metadata["runtime"],c=col)
    end
    xlabel("Hilbert Space Dimension")
    ylabel("Runtime")
    xscale("log")
    yscale("log")
end

# pretty plot of N=4 4x8 and 8x4 on torus for poster
function pretty_spectrum_closedopen(which_one::String="closed",anis::Float64=1.0)
    if which_one == "closed"
        params_dict = Dict([("Lx",4),("Ly",8),("N",4),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",anis)])
    elseif which_one == "open"
        params_dict = Dict([("Lx",8),("Ly",4),("N",4),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",anis)])
    else
        error("Invalid which_one")
    end
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    cols = ["b","r","k","k","k","k","k","k","k","k"]

    nrgs::Dict{String,Vector} = Dict([("intstrens",[])])
    etop::Int = params_dict["Lx"] == 8 ? 10 : 5
    for i in 1:etop
        nrgs["E$(i-1)"] = []
    end
    
    for f in files
        d,m = read_data_jld2(f,dataloc; output_level=0)
        if m["twist_angle"] != [0.0,0.0]
            continue
        end
        intstren = m["U"][end]
        if intstren > 2.0 && anis == 1.0
            continue
        end
        append!(nrgs["intstrens"],intstren)
        for i in 1:etop
            append!(nrgs["E$(i-1)"],d["nrg"][i] - d["nrg"][1])
        end
    end
    
    for i in 1:etop
        if i > 3
            scatter(nrgs["intstrens"],nrgs["E$(i-1)"],c=cols[i])
        elseif i == 3
            scatter(nrgs["intstrens"],nrgs["E$(i-1)"],c=cols[i],label="E2-$(etop-1)")
        else
            scatter(nrgs["intstrens"],nrgs["E$(i-1)"],c=cols[i],label="E$(i-1)")
        end
    end

    if which_one == "closed" && anis == 1.0
        plot([1.0,1.0],[-0.1,0.5],c="k",linestyle="--")
    end

    legend()
    ylim([-0.01,0.4])
    xlabel("Long Range Interaction Strength V/t")
    ylabel("E - E0")
    anis == 1.0 ? title("Energy Spectrum N=4 $(params_dict["Lx"])x$(params_dict["Ly"])") : title("Energy Spectrum N=4 $(params_dict["Lx"])x$(params_dict["Ly"]) Hopping Anisotropy=$(anis)")
end

# quick plot energies of N=4 4x8 and 8x4 on torus
function quick_spectrum_closedopen(which_one::String="closed")
    if which_one == "closed"
        params_dict = Dict([("Lx",4),("Ly",8),("N",4),("if_periodic_x",true),("if_periodic_y",true)])
    elseif which_one == "open"
        params_dict = Dict([("Lx",4),("Ly",8),("N",4),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    else
        error("Invalid which_one")
    end
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)


    cols = ["b","g","r","m","c"]
    if 10 > length(cols)
        cols = repeat(cols,ceil(Int,10/length(cols)))
    end

    
    for f in files
        d,m = read_data_jld2(f,dataloc)
        intstren = m["U"][end]
        if intstren > 2.0
            continue
        end
        for i in 1:length(d["nrg"])
            scatter(intstren,d["nrg"][i] - d["nrg"][1],c=cols[i])
        end
    end
    xlabel("Interaction Strength")
    ylabel("Energy")

end

#=function plot_spectrum(xxs::Vector,nrgs::Vector,idx::Int,nev::Int,xstring::AbstractString="x",if_diff::Bool=true; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    if_labels = get(kwargs,:if_labels,false)

    cols = ["b","g","r","m","c"]
    if nev > length(cols)
        cols = repeat(cols,ceil(Int,nev/length(cols)))
    end

    for i in 1:nev
        change = abs(xxs[1] - xxs[2])
        xval = xxs[idx]
        shift = (i - nev/2) * ((0.1*change)/(nev/2))
        idx == 1 && if_labels ? scatter(xval + shift,nrgs[i] - if_diff*nrgs[1],c=cols[i],label="E$(i-1)") : scatter(xval + shift,nrgs[i] - if_diff*nrgs[1],c=cols[i])
    end
    xlabel(xstring)
    ystring = if_diff ? "NRG - E0" : "NRG"
    ylabel(ystring)
    title("Energy Spectrum"*plot_title)
    idx == 1 && if_labels ? legend() : nothing

    return
end
plot_spectrum(xxs::StepRangeLen,nrgs::Vector,idx::Int,nev::Int,xstring::AbstractString="x",if_diff::Bool=true; kwargs...) = plot_spectrum(collect(xxs),nrgs,idx,nev,xstring,if_diff; kwargs...)=#

#=function plot_omega_ed(theta_xs::Vector{Float64},theta_ys::Vector{Float64},omegas::Matrix{ComplexF64}; kwargs...)
    plot_title::String = get(kwargs,:plot_title,"")

    omegas_phase::Matrix{Float64} = zeros(Float64,length(theta_xs),length(theta_ys))
    for i in 1:length(theta_xs)
        for j in 1:length(theta_ys)
            omegas_phase[i,j] = angle(omegas[i,j]) + pi
        end
    end
    omegas_phase[1,1] = 0.0
    reverse!(omegas_phase, dims=1)

    fig = figure()
    imshow(omegas_phase; cmap="hsv", extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)], vmax=2*pi, vmin=0.0)
    colorbar()
    reverse!(omegas_phase,dims=1)
    diag_shift = 0.0 * maximum(theta_xs) / length(theta_xs)
    xs = transpose(repeat(theta_xs,1,length(theta_ys))) .+ diag_shift
    ys = repeat(theta_ys,1,length(theta_xs)) .+ diag_shift
    us = cos.(omegas_phase)
    vs = sin.(omegas_phase)
    quiver(xs, ys, us, vs)
    title("Phase of "*L"\Omega"*plot_title)
    xlabel(L"\theta_x / 2\pi")
    ylabel(L"\theta_y / 2\pi")
    xlim([minimum(theta_xs),maximum(theta_xs)])
    ylim([minimum(theta_ys),maximum(theta_ys)])

end
plot_omega_ed(theta_xs::StepRangeLen,theta_ys::StepRangeLen,omegas::Matrix{ComplexF64}; kwargs...) = plot_omega(collect(theta_xs),collect(theta_ys),omegas; kwargs...)=#

function plot3d_flow_with_levelcrossing_3x7n3()
    tw2 = 0.0
    amplification_factor = 12.0
    nev = 6

    cols = ["b","g","r","m","c"]
    if nev > length(cols)
        cols = repeat(cols,ceil(Int,nev/length(cols)))
    end

    intstrens = range(0.0,2.0,length=10)
    tws = range(0.0,1.0,length=30)
    these_nrgs = zeros(Float64,nev,length(tws))
    for (idx2,intstren) in enumerate(intstrens)
        where_fqh = 1
        if intstren > 1.5
            where_fqh += 2
        end
        if intstren > 2
            where_fqh += 2
        end
        for (idx,tw1) in enumerate(tws)
            params_dict = Dict([("Lx",3),("Ly",7),("N",3),("interaction_strength",intstren),("nev",nev),("if_save_data",false),("if_find_data",false),("if_periodic_x",true),("if_periodic_y",true),("tw1",tw1),("tw2",tw2)])
            states,nrgs,rhos,filepath,if_found,latpara,hamiltpara = run_normal_ed(params_dict; output_level=1)
            these_nrgs[:,idx] = nrgs
        end

        # this next section will amplify the flow by the given factor
        plotting_nrgs = these_nrgs
        max_fqh_nrg = maximum(plotting_nrgs[where_fqh,2:end-1])
        min_fqh_nrg = minimum(plotting_nrgs[where_fqh,2:end-1])
        middle_fqh_nrg = (max_fqh_nrg + min_fqh_nrg)/2
        plotting_nrgs[where_fqh:where_fqh+1,:] = (plotting_nrgs[where_fqh:where_fqh+1,:] .- middle_fqh_nrg) .* amplification_factor .+ middle_fqh_nrg

        botshift = where_fqh == 1 ? middle_fqh_nrg : minimum(plotting_nrgs)
        for i in 1:nev
            scatter3D(tws,intstren .* ones(length(tws)),plotting_nrgs[i,:] .- botshift,c=cols[i])
        end
        xlabel("Theta_x / 2pi")
        ylabel("Interaction Strength")
        zlabel("Energy")
        #zlim([0.0,0.1])
        title("Spectrum for Lx=3 Ly=7 N=3 ULR=$intstren")
    end

end

# make nice figures for presentation for spectral flow twisting
function pretty_plot_twisting(points_count::Int,if_gap::Bool=true,if_diff::Bool=false)
    lx,ly,n = 6,5,3
    tw2 = 0.0

    cols = get_colors(3)

    tws = range(0.0,1.0,length=points_count)
    for (idx,tw1) in enumerate(tws)
        params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",0.0),("nev",3),("if_save_data",false),("if_find_data",false),("if_periodic_x",true),("if_periodic_y",true),("tw1",tw1),("tw2",tw2)])
        states,nrgs,rhos,filepath,if_found,latpara,hamiltpara = run_normal_ed(params_dict; output_level=1)
        if if_gap
            plot_spectrum(tws,nrgs,idx,3,L"\theta_{x} / 2 \pi ",if_diff; if_labels=true)
        else
            plot_spectrum(tws,nrgs,idx,2,L"\theta_{x} / 2 \pi ",if_diff; if_labels=true)
        end
    end

end

# plots the spectrum as a function of interaction strength for all available data at zero twist angle
function plot_intstren_spectrum(lx::Int64,ly::Int64,n::Int64; kwargs...)
    hanis = get(kwargs,:hanis,1.0)    
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0)

    intstrens = []
    nrgs = Dict()
    for i in 1:20
        nrgs[string(i)] = []
    end
    for f in all_files
        filename_dict = get_params_dict_from_filename(f)
        if haskey(filename_dict,"twist_angle1")
            continue
        end
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        for i in 1:20
            if i > length(d["nrg"])
                append!(nrgs[string(i)],1000.0)
            else
                append!(nrgs[string(i)],d["nrg"][i])
            end
        end
        append!(intstrens,filename_dict["interaction_strength"])
    end
    plot_fullspectrum(intstrens,nrgs,"Interaction Strength",true; plot_title=" $(lx)x$(ly) N=$(n)")

end

# plot phase diagram ULR vs rho1D using flatness
function make_phasediag_ulrrho1d_flatness(configs=[(3,8,3),(8,3,3),(4,6,3),(8,4,4)]; kwargs...)
    max_intstren::Float64 = get(kwargs,:max_intstren,2.0)
    if_plot::Bool = get(kwargs,:if_plot,true)

    ulrs::Vector{Float64} = Float64[]
    flatnesses::Vector{Float64} = Float64[]
    oneDrhos::Vector{Float64} = Float64[]

    for (lx,ly,n) in configs
        local_strens,local_flats = twist_flatness_ed(lx,ly,n; if_plot=false, max_intstren=max_intstren)
        append!(ulrs,local_strens)
        append!(flatnesses,local_flats)
        append!(oneDrhos,ones(Float64,length(local_strens)) .* (n / lx))
    end

    bin_count = get(kwargs,:bin_count,100)
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
    ylim([-0.05,max_intstren + 0.05])=#

    if_plot ? plot_phasediag_ulrrho1d_flatness(oneDrhos,ulrs,normalized_bv,max_intstren) : nothing

    return oneDrhos,ulrs,normalized_bv,max_intstren
end

function plot_phasediag_ulrrho1d_flatness(rho1Ds::Vector{Float64},ulrs::Vector{Float64},normalized_bv,max_intstren::Float64)
    
    linear_range = (-0.05,2.5)
    log_range = (2.5,max_intstren)
    # Create a figure with two subplot axes
    fig = figure(figsize=(8,6))
    
    # Create two axes with specified height ratios
    gs = plt.GridSpec(2, 1, height_ratios=[1, 1], hspace=0.05)
    ax1 = fig.add_subplot(gs[0])  # upper subplot
    ax2 = fig.add_subplot(gs[1])  # lower subplot
    
    # Plot data on both axes
    ax1.scatter(rho1Ds, ulrs, c=normalized_bv, cmap="viridis")
    ax2.scatter(rho1Ds, ulrs, c=normalized_bv, cmap="viridis")
    ax2.plot([minimum(rho1Ds),maximum(rho1Ds)],[maximum(ulrs),maximum(ulrs)],c="r",linestyle="--",label=L"\infty")
    
    # Set the ranges for each axis
    ax1.set_ylim(linear_range)  # log scale portion
    ax2.set_ylim(log_range)  # linear scale portion
    
    # Set scale for each axis
    ax1.set_yscale("linear")
    ax2.set_yscale("log")
    
    # Hide the appropriate tick marks
    ax2.spines["bottom"].set_visible(false)
    ax1.spines["top"].set_visible(false)
    ax2.xaxis.set_visible(false)

    ylabel("ULR")
    xlabel(L"\rho_{1D}")
    title("Flatness Phase Diagram")
    legend()

end

function plot_hatsugai_fromsaveddata(lx::Int64,ly::Int64,N::Int64; kwargs...)
    intstren::Float64 = get(kwargs,:intstren,0.0)
    hanis::Float64 = get(kwargs,:hanis,1.0)
    if_pinning::Bool = get(kwargs,:if_pinning,false)

    dataloc::String = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",N),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis),("interaction_strength",intstren),("if_pinning",if_pinning)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0)

    tw1s::Vector{Float64} = Float64[]
    tw2s::Vector{Float64} = Float64[]
    omegas::Vector{ComplexF64} = ComplexF64[]
    lambda1s::Vector{ComplexF64} = ComplexF64[]
    lambda2s::Vector{ComplexF64} = ComplexF64[]
    for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        if haskey(m,"omega")
            append!(tw1s,m["twist_angle"][1])
            append!(tw2s,m["twist_angle"][2])
            append!(omegas,m["omega"])
            append!(lambda1s,m["gamma1"])
            append!(lambda2s,m["gamma2"])
        else
            println("No omega data found")
        end
    end

    square_size = sqrt(length(tw1s))
    unique!(tw1s)
    unique!(tw2s)
    if isinteger(square_size)
        square_size = Int(square_size)
        plot_omega(tw1s,tw2s,reshape(omegas,square_size,square_size))
        plot_gamma(tw1s,tw2s,reshape(lambda1s,square_size,square_size),1)
        plot_gamma(tw1s,tw2s,reshape(lambda2s,square_size,square_size),2)
    else
        error("Not a square number of data points, size of data is $(length(tw1s))")
    end

end





























"fin"
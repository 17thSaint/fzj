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

using PyPlot,LaTeXStrings

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

function plot_gamma_fromsaveddata_scatter(lx::Int,ly::Int,N::Int,which_gamma::Int; kwargs...)
    intstren = get(kwargs,:interaction_strength,0.0)
    hanis = get(kwargs,:hopping_anisotropy,1.0)
    ref_twists = get(kwargs,:ref_twists,[0.33,0.67])
    bin_count = get(kwargs,:bin_count,100)

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    params_dict = Dict([("Lx",lx),("Ly",ly),("N",N),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("hopping_anisotropy",hanis)])
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    theta_xs::Vector{Float64} = []
    theta_ys::Vector{Float64} = []
    gammas_mag::Vector{Float64} = []
    #gammas_phase::Vector{Float64} = []
    for f in files
        filename_dict = get_params_dict_from_filename(f)
        try
            if filename_dict["twist_angle1"] in ref_twists || filename_dict["twist_angle2"] in ref_twists
                continue
            end
        catch
            continue
        end

        data,metadata = read_data_jld2(f,dataloc;output_level=0)
        append!(theta_xs,filename_dict["twist_angle1"])
        append!(theta_ys,filename_dict["twist_angle2"])
        append!(gammas_mag,abs(metadata["gamma$(which_gamma)"]))
        #append!(gammas_phase,angle(metadata["gamma$(which_gamma)"]))
    end
    
    data_dict_mag = bin_values(gammas_mag,bin_count)
    gmag = [data_dict_mag[val] for val in gammas_mag]
    min_gammas_mag, max_gammas_mag = minimum(gammas_mag), maximum(gammas_mag)
    normalized_gmag = [(val - minimum(gmag)) / (maximum(gmag) - minimum(gmag)) * (max_gammas_mag - min_gammas_mag) + min_gammas_mag for val in gmag]

    fig = figure()
    scatter(theta_xs, theta_ys, c=normalized_gmag, cmap="viridis")
    colorbar()
    xlabel("Theta_x / 2pi")
    ylabel("Theta_y / 2pi")
    title("Magnitude of Gamma$(which_gamma)")
end

function count_chern_number_scatter(theta_xs::Vector{Float64},theta_ys::Vector{Float64},omega_phases::Vector{Float64}; kwargs...)
    plot_title::String = get(kwargs,:plot_title,"")
    
    append!(theta_xs,0.0)
    append!(theta_ys,0.0)
    append!(omega_phases,0.0)
    omega_phase_mat = reshape(omega_phases,Int(sqrt(length(theta_xs))),Int(sqrt(length(theta_xs))))

    angle_diffs_x::Matrix{Float64} = zeros(Float64,size(omega_phase_mat))
    angle_diffs_y::Matrix{Float64} = zeros(Float64,size(omega_phase_mat))
    for i in 2:size(omega_phase_mat,1) - 1
        for j in 2:size(omega_phase_mat,2) - 1
            angle_diffs_x[i,j] += (omega_phase_mat[i+1,j] - omega_phase_mat[i-1,j]) / (2*pi)
            angle_diffs_y[i,j] += (omega_phase_mat[i,j+1] - omega_phase_mat[i,j-1]) / (2*pi)
        end
    end
    
    fig = figure()
    imshow(angle_diffs_x; cmap="jet")
    colorbar()
    title("Omega X Angle Difference"*plot_title)

    fig = figure()
    imshow(angle_diffs_y; cmap="jet")
    colorbar()
    title("Omega Y Angle Difference"*plot_title)
end

function plot_omega_fromsaveddata_scatter(lx::Int,ly::Int,N::Int; kwargs...)
    intstren = get(kwargs,:interaction_strength,0.0)
    hanis = get(kwargs,:hopping_anisotropy,1.0)
    ref_twists = get(kwargs,:ref_twists,[0.33,0.67])
    bin_count = get(kwargs,:bin_count,100)
    plot_title = get(kwargs,:plot_title,"")
    if_angle_diff = get(kwargs,:if_angle_diff,true)
    if_omega_mag = get(kwargs,:if_omega_mag,false)

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    params_dict = Dict([("Lx",lx),("Ly",ly),("N",N),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("hopping_anisotropy",hanis)])
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    theta_xs::Vector{Float64} = []
    theta_ys::Vector{Float64} = []
    omegas_mag::Vector{Float64} = []
    omegas_phase::Vector{Float64} = []
    for f in files
        filename_dict = get_params_dict_from_filename(f)
        try
            if filename_dict["twist_angle1"] in ref_twists || filename_dict["twist_angle2"] in ref_twists
                continue
            end
        catch
            continue
        end
        data,metadata = read_data_jld2(f,dataloc;output_level=0)
        append!(theta_xs,filename_dict["twist_angle1"])
        append!(theta_ys,filename_dict["twist_angle2"])
        append!(omegas_mag,abs(metadata["omega"]))
        append!(omegas_phase,angle(metadata["omega"]))
    end

    if if_omega_mag
        data_dict_mag = bin_values(omegas_mag,bin_count)
        omag = [data_dict_mag[val] for val in omegas_mag]
        min_omegas_mag, max_omegas_mag = minimum(omegas_mag), maximum(omegas_mag)
        normalized_omag = [(val - minimum(omag)) / (maximum(omag) - minimum(omag)) * (max_omegas_mag - min_omegas_mag) + min_omegas_mag for val in omag]

        fig = figure()
        scatter(theta_xs, theta_ys, c=normalized_omag, cmap="viridis")
        colorbar()
        title("Magnitude of Omega"*plot_title)
        xlabel("Theta_x / 2pi")
        ylabel("Theta_y / 2pi")
    end

    data_dict_angle = bin_values(omegas_phase,bin_count)
    oang = [data_dict_angle[val] for val in omegas_phase]
    min_omegas_phase, max_omegas_phase = minimum(omegas_phase), maximum(omegas_phase)
    normalized_oang = [(val - minimum(oang)) / (maximum(oang) - minimum(oang)) * (max_omegas_phase - min_omegas_phase) + min_omegas_phase for val in oang]

    fig = figure()
    scatter(theta_xs, theta_ys, c=normalized_oang, cmap="jet")
    colorbar()
    title("Phase of Omega"*plot_title)
    xlabel("Theta_x / 2pi")
    ylabel("Theta_y / 2pi")

    if_angle_diff && count_chern_number(theta_xs,theta_ys,omegas_phase)
end

function plot_gamma_fromsaveddata(lx::Int,ly::Int,N::Int,which_gamma::Int; kwargs...)
    intstren = get(kwargs,:interaction_strength,0.0)
    hanis = get(kwargs,:hopping_anisotropy,1.0)
    ref_twists = get(kwargs,:ref_twists,[0.33,0.67])
    plot_title = get(kwargs,:plot_title,"")

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    params_dict = Dict([("Lx",lx),("Ly",ly),("N",N),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("hopping_anisotropy",hanis)])
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    length(files) == 0 && error("No data found for given parameters")

    theta_xs::Vector{Float64} = []
    theta_ys::Vector{Float64} = []
    for f in files
        filename_dict = get_params_dict_from_filename(f)
        try
            if filename_dict["twist_angle1"] in ref_twists || filename_dict["twist_angle2"] in ref_twists
                continue
            end
        catch
            continue
        end
        if !(filename_dict["twist_angle1"] in theta_xs)
            append!(theta_xs,filename_dict["twist_angle1"])
            sort!(theta_xs)
        end
        if !(filename_dict["twist_angle2"] in theta_ys)
            append!(theta_ys,filename_dict["twist_angle2"])
            sort!(theta_ys; rev=true)
        end
    end

    gammas::Matrix{Float64} = zeros(Float64,length(theta_xs),length(theta_ys))
    for f in files
        filename_dict = get_params_dict_from_filename(f)
        try
            if filename_dict["twist_angle1"] in ref_twists || filename_dict["twist_angle2"] in ref_twists
                continue
            end
        catch
            continue
        end
        data,metadata = read_data_jld2(f,dataloc;output_level=0)
        i = findfirst(x -> x == filename_dict["twist_angle1"],theta_xs)
        j = findfirst(x -> x == filename_dict["twist_angle2"],theta_ys)
        gammas[i,j] = abs(metadata["gamma$(which_gamma)"])
    end
   
    fig = figure()
    imshow(transpose(gammas); extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)], vmin=0.0, vmax=1.0)
    colorbar()
    xlabel(L"\theta_x / 2\pi")
    ylabel(L"\theta_y / 2\pi")
    which_gamma == 1 ? title("Magnitude of "*L"\Lambda_1" *plot_title) : title("Magnitude of "*L"\Lambda_2" *plot_title)
end

function plot_omega_fromsaveddata(lx::Int,ly::Int,N::Int; kwargs...)
    intstren = get(kwargs,:interaction_strength,0.0)
    hanis = get(kwargs,:hopping_anisotropy,1.0)
    ref_twists = get(kwargs,:ref_twists,[0.33,0.67])
    bin_count = get(kwargs,:bin_count,100)
    plot_title = get(kwargs,:plot_title,"")
    if_angle_diff = get(kwargs,:if_angle_diff,true)
    if_omega_mag = get(kwargs,:if_omega_mag,false)
    if_save_fig = get(kwargs,:if_save_fig,false)

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    params_dict = Dict([("Lx",lx),("Ly",ly),("N",N),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("hopping_anisotropy",hanis)])
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    length(files) == 0 && error("No data found for given parameters")

    theta_xs::Vector{Float64} = []
    theta_ys::Vector{Float64} = []
    for f in files
        filename_dict = get_params_dict_from_filename(f)
        try
            if filename_dict["twist_angle1"] in ref_twists || filename_dict["twist_angle2"] in ref_twists
                continue
            end
        catch
            continue
        end
        if !(filename_dict["twist_angle1"] in theta_xs)
            append!(theta_xs,filename_dict["twist_angle1"])
            sort!(theta_xs)
        end
        if !(filename_dict["twist_angle2"] in theta_ys)
            append!(theta_ys,filename_dict["twist_angle2"])
            sort!(theta_ys; rev=true)
        end
    end

    omegas_phase::Matrix{Float64} = zeros(Float64,length(theta_xs),length(theta_ys))
    for f in files
        filename_dict = get_params_dict_from_filename(f)
        try
            if filename_dict["twist_angle1"] in ref_twists || filename_dict["twist_angle2"] in ref_twists
                continue
            end
        catch
            continue
        end
        data,metadata = read_data_jld2(f,dataloc;output_level=0)
        i = findfirst(x -> x == filename_dict["twist_angle1"],theta_xs)
        j = findfirst(x -> x == filename_dict["twist_angle2"],theta_ys)
        omegas_phase[i,j] = angle(metadata["omega"]) + pi
    end

    fig = figure()
    imshow(transpose(omegas_phase); cmap="hsv", extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)])
    colorbar()
    diag_shift = 0.5 * maximum(theta_xs) / length(theta_xs)
    xs = transpose(repeat(theta_xs,1,length(theta_ys))) .+ diag_shift
    ys = repeat(theta_ys,1,length(theta_xs)) .+ diag_shift
    us = cos.(transpose(omegas_phase))
    vs = sin.(transpose(omegas_phase))
    quiver(xs, ys, us, vs)
    title("Phase of "*L"\Omega"*plot_title)
    xlabel(L"\theta_x / 2\pi")
    ylabel(L"\theta_y / 2\pi")
    xlim([minimum(theta_xs),maximum(theta_xs)])
    ylim([minimum(theta_ys),maximum(theta_ys)])

    if_save_fig && fig.savefig("local-figs/omegaphase-mbcn1p0-6x5-n3-ulr0p0-arrowsheat-morepixels.png",dpi=300)

    #if_angle_diff ? count_chern_number(theta_xs,theta_ys,omegas_phase) : nothing
end

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

# plotting the minimum amount the GSs go into the excited states as a function of interaction strength
function plot_twistflatness_vs_intstren_ed(lx,ly,n; kwargs...)
    intstrens = get(kwargs,:intstrens, range(0.0,2.0,length=11))
    hanis = get(kwargs,:hanis,1.0)
    if_plot = get(kwargs,:if_plot,true)
    if_plot_flatness = get(kwargs,:if_plot_flatness,false)
    plot_title = get(kwargs,:plot_title,"")
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    
    flatnesses = zeros(Float64,length(intstrens))

    for (idx,intstren) in enumerate(intstrens)
        params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis)])
        all_files = find_data_file(params_dict,"ed",dataloc; output_level=0)

        tw1s = Float64[]
        tw2s = Float64[]
        nrgs = Dict([("1",Float64[]),("2",Float64[]),("3",Float64[])])
        for f in all_files
            d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
            append!(tw1s,m["twist_angle"][1])
            append!(tw2s,m["twist_angle"][2])
            for i in 1:3
                append!(nrgs[string(i)],d["nrg"][i])
            end
        end

        all_flatnesses = (nrgs["2"] .- nrgs["1"]) ./ (nrgs["3"] .- nrgs["1"])
        flatnesses[idx] = maximum(all_flatnesses)

        if if_plot_flatness
            fig = figure()
            scatter3D(tw1s,tw2s,all_flatnesses,c="b")
            #scatter3D(tw1s,tw2s,nrgs["1"],c="b")
            #scatter3D(tw1s,tw2s,nrgs["2"],c="g")
            #scatter3D(tw1s,tw2s,nrgs["3"],c="r")
            xlabel(L"\theta_x / 2\pi")
            ylabel(L"\theta_y / 2\pi")
            title("Flatness for $(lx)x$(ly) N=$n ULR=$intstren ")
        end
    end

    if if_plot
        fig = figure()
        scatter(intstrens,flatnesses,c="b")
        xlabel("Interaction Strength")
        ylabel("Minimum Twist Flatness")
        ylim(-0.02,1.05)
        title("Flatness for $(lx)x$(ly) N=$n "*plot_title)
    end

    return intstrens,flatnesses
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





























"fin"
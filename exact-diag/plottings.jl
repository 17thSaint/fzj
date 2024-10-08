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
function pretty_spectrum_closedopen(which_one::String="closed")
    if which_one == "closed"
        params_dict = Dict([("Lx",4),("Ly",8),("N",4),("if_periodic_x",true),("if_periodic_y",true)])
    elseif which_one == "open"
        params_dict = Dict([("Lx",8),("Ly",4),("N",4),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
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
        d,m = read_data_jld2(f,dataloc)
        intstren = m["U"][end]
        if intstren > 2.0
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
    legend()
    ylim([-0.01,0.4])
    xlabel("Long Range Interaction Strength V/t")
    ylabel("E - E0")
    title("Energy Spectrum N=4 $(params_dict["Lx"])x$(params_dict["Ly"])")
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

function plot_spectrum(xxs::Vector,nrgs::Vector,idx::Int,nev::Int,xstring::String="x",if_diff::Bool=true; kwargs...)
    plot_title = get(kwargs,:plot_title,"")

    cols = ["b","g","r","m","c"]
    if nev > length(cols)
        cols = repeat(cols,ceil(Int,nev/length(cols)))
    end

    for i in 1:nev
        change = abs(xxs[1] - xxs[2])
        xval = xxs[idx]
        shift = (i - nev/2) * ((0.1*change)/(nev/2))
        scatter(xval + shift,nrgs[i] - if_diff*nrgs[1],c=cols[i])
    end
    xlabel(xstring)
    ystring = if_diff ? "NRG - E0" : "NRG"
    ylabel(ystring)
    title("Energy Spectrum"*plot_title)

    return
end

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

function count_chern_number(theta_xs::Vector{Float64},theta_ys::Vector{Float64},omega_phases::Matrix{Float64}; kwargs...)
    plot_title::String = get(kwargs,:plot_title,"")

    angle_diffs::Matrix{Float64} = zeros(Float64,size(omega_phases))
    #angle_diffs_y::Matrix{Float64} = zeros(Float64,size(omega_phases))
    for i in 2:size(omega_phases,1) - 1
        for j in 2:size(omega_phases,2) - 1
            for xshift in [1,-1]
                for yshift in [1,-1]
                    angle_diffs[i,j] += (omega_phases[i+xshift,j+yshift] - omega_phases[i,j]) / (2*pi)
                end
            end
        end
    end
    
    fig = figure()
    imshow(transpose(angle_diffs); cmap="jet", extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)])
    colorbar()
    title("Omega Angle Difference"*plot_title)
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
    imshow(transpose(gammas); extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)])
    colorbar()
    xlabel("Theta_x / 2pi")
    ylabel("Theta_y / 2pi")
    title("Magnitude of Gamma$(which_gamma)"*plot_title)
end

function plot_omega_fromsaveddata(lx::Int,ly::Int,N::Int; kwargs...)
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
    title("Phase of Omega"*plot_title)
    xlabel("Theta_x / 2pi")
    ylabel("Theta_y / 2pi")
    xlim([minimum(theta_xs),maximum(theta_xs)])
    ylim([minimum(theta_ys),maximum(theta_ys)])

    #if_angle_diff ? count_chern_number(theta_xs,theta_ys,omegas_phase) : nothing
end

function plot_omega(theta_xs::Vector{Float64},theta_ys::Vector{Float64},omegas::Matrix{ComplexF64}; kwargs...)
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
    imshow(omegas_phase; cmap="hsv", extent=[minimum(theta_xs),maximum(theta_xs),minimum(theta_ys),maximum(theta_ys)])
    colorbar()
    reverse!(omegas_phase,dims=1)
    diag_shift = 0.0 * maximum(theta_xs) / length(theta_xs)
    xs = transpose(repeat(theta_xs,1,length(theta_ys))) .+ diag_shift
    ys = repeat(theta_ys,1,length(theta_xs)) .+ diag_shift
    us = cos.(omegas_phase)
    vs = sin.(omegas_phase)
    quiver(xs, ys, us, vs)
    title("Phase of Omega"*plot_title)
    xlabel("Theta_x / 2pi")
    ylabel("Theta_y / 2pi")
    xlim([minimum(theta_xs),maximum(theta_xs)])
    ylim([minimum(theta_ys),maximum(theta_ys)])

end































"fin"
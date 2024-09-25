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

include("execute-ed.jl")



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

function plot_spectrum(xxs::Vector,nrgs::Vector,nev::Int,xstring::String="x",if_diff::Bool=true; kwargs...)
    plot_title = get(kwargs,:title,"")

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
    title(plot_title)

    return
end


































"fin"
#####################################################
#=

This file contains plotting functions for TTN results

Depends on:
    other-funcs/basic-2d-stuff.jl
    review-practice-codes/ttn.jl
    review-practice-codes/observables.jl

=#
######################################################

using PyPlot,LaTeXStrings

include("../other-funcs/include-other-files.jl")

include_other_files(["other-funcs/basic-2d-stuff.jl","review-practice-codes/ttn.jl","review-practice-codes/observables.jl","other-funcs/basic-2d-plottings.jl"])

#### Generic Plotting

#=function plot_occupancy(exp_occ; kwargs...)
    fix_colorbar = get(kwargs, :fix_colorbar, false)
	data_dict = get(kwargs, :data_dict, nothing)
	if !isnothing(data_dict)
		exp_occ = data_dict["vals"]
	end
	fig = figure()
	fix_colorbar ? imshow(transpose(exp_occ),vmin=0.0,vmax=maximum(exp_occ)) : imshow(transpose(exp_occ))
	colorbar()
	plot_title = get(kwargs, :plot_title, "")
	title_string = "Occupancy, " * plot_title
	title(title_string)
	ylabel("Synthetic")
	xlabel("Physical")
	
	if get(kwargs, :if_save_fig, false)
		location = get(kwargs, :location, pwd())
		fig_name = get(kwargs, :name, "occs")
		fig_name = check_plot_label(fig_name,"occs")
		save_figure(fig_name; kwargs...)
	end
	return
end=#

# needs testing
function plot_occupancy_3d(exp_occ; kwargs...)
	m,n = size(exp_occ)

	x = collect(1:m)  # x-coordinates (1 to m)
	y = collect(1:n)  # y-coordinates (1 to n)

	# Create the corresponding z-values from the matrix
	z = exp_occ

	# Now, we need to turn x, y, z into vectors for scatter3D
	# Use the function `repeat` and `vec` to flatten them
	x_flat = repeat(x, n)     # Repeat x n times
	y_flat = repeat(y', m)    # Repeat each y m times
	z_flat = vec(z)           # Flatten the matrix into a vector

	# Create the 3D scatter plot
	fig = figure()
	scatter3D(x_flat, y_flat, z_flat)
	return
end

function plot_synthetic_current(currents::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    which_lines = get(kwargs,:which_lines,1:size(currents)[2])
    fig = figure()
    for i in which_lines
        plot(1:size(currents)[1],currents[:,i],"-p",label="$i")
    end
    xlabel("Physical Site")
    ylabel("Current")
    title("Synthetic Current "*plot_title)
    legend()
    return nothing
end

function plot_physical_current(currents::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    which_lines = get(kwargs,:which_lines,1:size(currents)[1])
    fig = figure()
    for i in which_lines
        plot(1:size(currents)[2],currents[i,:],"-p",label="$i")
    end
    xlabel("Synthetic Site")
    ylabel("Current")
    title("Physical Current "*plot_title)
    legend()
    return nothing
end

function plot_synthetic_correlation(syn_corrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    which_lines = get(kwargs,:which_lines,1:size(syn_corrs)[1])
    fig = figure()
    for i in which_lines
        plot(1:size(syn_corrs)[2],syn_corrs[i,:],"-p",label="$i")
    end
    xlabel("Synthetic Distance")
    ylabel("Correlation")
    title("Synthetic Correlation "*plot_title)
    legend()
    yscale("log")
    return nothing
end

function plot_physical_correlation(phys_corrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    which_lines = get(kwargs,:which_lines,1:size(phys_corrs)[2])
    fig = figure()
    for i in which_lines
        plot(1:size(phys_corrs)[1],phys_corrs[:,i],"-p",label="$i")
    end
    xlabel("Physical Distance")
    ylabel("Correlation")
    title("Physical Correlation "*plot_title)
    legend()
    yscale("log")
    return nothing
end




# make pretty plots for OBC FQH droplet
function prettyplot_obc_occupancy(layers::Int64; kwargs...)
    Lx,Ly = get_lattice_dims_from_layers(layers)
    intstren = get(kwargs, :onsite_strength, 0.0)
    particles = get(kwargs, :particles, Int(sqrt(2^layers)/2))
    hanis = get(kwargs, :hopping_anisotropy, 1.0)

    dataloc = get_folder_location("cluster-data/synth-dims/obc")
    params_dict = Dict([("layers",layers),("onsite_strength",intstren),("particles",particles),("hopping_anisotropy",hanis)])

    filename = find_data_file(params_dict,"ttn",dataloc)[1]
    d,m = read_data_jld2(filename,dataloc; output_level=0)

    occs = get_occupancy(d["densmat"]; if_plot=false)

    fig = figure()
    imshow(transpose(occs),vmin=0.0,vmax=maximum(occs),extent=[0.5,Lx+0.5,0.5,Ly+0.5])
    colorbar()
    title("Onsite Density")
    ylabel("Synthetic")
    xlabel("Physical")

    for i in [1,2,4,8]
        plot(1:Lx,i .* ones(Lx))
    end
end

function prettyplot_obc_correlations(layers::Int64; kwargs...)
    intstren = get(kwargs, :onsite_strength, 0.0)
    particles = get(kwargs, :particles, Int(sqrt(2^layers)/2))
    hanis = get(kwargs, :hopping_anisotropy, 1.0)

    Lx,Ly = Int(sqrt(2^layers)),Int(sqrt(2^layers))

    dataloc = get_folder_location("cluster-data/synth-dims/obc")
    params_dict = Dict([("layers",layers),("onsite_strength",intstren),("particles",particles),("hopping_anisotropy",hanis)])

    filename = find_data_file(params_dict,"ttn",dataloc)[1]
    d,m = read_data_jld2(filename,dataloc; output_level=0)

    phys_corrs = physical_correlation(d["densmat"],Lx,Ly; if_plot=true,which_lines=[1,2,4,8])
    title("Correlations")
    xlabel("Distance (j - 1)")
    ylabel(L"\vert \langle a_{1,s}^{\dagger} a_{j,s} \rangle \vert")
    #syn_corrs = synthetic_correlation(d["densmat"],Lx,Ly; if_plot=true)
end

function prettyplot_obc_currents(layers::Int64; kwargs...)
    intstren = get(kwargs, :onsite_strength, 0.0)
    particles = get(kwargs, :particles, Int(sqrt(2^layers)/2))
    hanis = get(kwargs, :hopping_anisotropy, 1.0)

    Lx,Ly = Int(sqrt(2^layers)),Int(sqrt(2^layers))

    dataloc = get_folder_location("cluster-data/synth-dims/obc")
    params_dict = Dict([("layers",layers),("onsite_strength",intstren),("particles",particles),("hopping_anisotropy",hanis)])

    filename = find_data_file(params_dict,"ttn",dataloc)[1]
    d,m = read_data_jld2(filename,dataloc; output_level=0)

    syn_currents = synthetic_current(d["densmat"],Lx,Ly; if_plot=true, which_lines=[1,2,4,8])
    title("Currents")
    xlabel("Site (j)")
    ylabel(L"i\langle a_{j,s+1}^{\dagger} a_{j,s} - a_{j+1,s} a_{j,s+1}^{\dagger} \rangle")
    #phys_currents = physical_current(d["densmat"],Lx,Ly; if_plot=true)
end

function prettyplot_obc_all(layers::Int64; kwargs...)
    Lx,Ly = get_lattice_dims_from_layers(layers)
    intstren = get(kwargs, :onsite_strength, 0.0)
    particles = get(kwargs, :particles, Int(sqrt(2^layers)/2))
    hanis = get(kwargs, :hopping_anisotropy, 1.0)
    which_lines = get(kwargs, :which_lines, [1,2,4,8])

    dataloc = get_folder_location("cluster-data/synth-dims/obc")
    params_dict = Dict([("layers",layers),("onsite_strength",intstren),("particles",particles),("hopping_anisotropy",hanis)])

    filename = find_data_file(params_dict,"ttn",dataloc)[1]
    d,m = read_data_jld2(filename,dataloc; output_level=0)

    # occupancy matrix
    occs = get_occupancy(d["densmat"]; if_plot=false)
    fig = figure()
    imshow(transpose(occs),vmin=0.0,vmax=maximum(occs),extent=[0.5,Lx+0.5,0.5,Ly+0.5])
    colorbar()
    title("Onsite Density"; fontsize=16)
    ylabel("Synthetic"; fontsize=16)
    xlabel("Physical"; fontsize=16)
    for i in which_lines
        plot(1:Lx,i .* ones(Lx))
    end
    withlines = length(which_lines) > 0 ? "true" : "false"
    #savefig("~/fzj/main-git/synth-dims/local-figs/prettyoccupancy-layers-$(layers)-withlines-$withlines.png",dpi=300)

    # Correlations
    phys_corrs = physical_correlation(d["densmat"],Lx,Ly; if_plot=true,which_lines=which_lines)
    title("Correlations"; fontsize=16)
    xlabel("Distance (j)"; fontsize=16)
    ylabel(L"\vert \langle a_{1,s}^{\dagger} a_{j,s} \rangle \vert"; fontsize=16)
    #savefig("~/fzj/main-git/synth-dims/local-figs/prettycorrelations-layers-$(layers).png",dpi=300)

    # Currents
    syn_currents = synthetic_current(d["densmat"],Lx,Ly; if_plot=true, which_lines=which_lines)
    title("Currents"; fontsize=16)
    xlabel("Site (j)"; fontsize=16)
    ylabel(L"i\langle a_{j,s+1}^{\dagger} a_{j,s} - a_{j+1,s} a_{j,s+1}^{\dagger} \rangle"; fontsize=16)
    #savefig("~/fzj/main-git/synth-dims/local-figs/prettycurrents-layers-$(layers).png",dpi=300)
end

#prettyplot_obc_all(8; which_lines=[3,4,6,9])























"fin"
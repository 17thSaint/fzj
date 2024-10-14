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

include_other_files(["other-funcs/basic-2d-stuff.jl","review-practice-codes/ttn.jl","review-practice-codes/observables.jl"])

#### Generic Plotting

function plot_occupancy(exp_occ; kwargs...)
	data_dict = get(kwargs, :data_dict, nothing)
	if !isnothing(data_dict)
		exp_occ = data_dict["vals"]
	end
	fig = figure()
	imshow(exp_occ)#,vmin=0.0,vmax=maximum(exp_occ))
	colorbar()
	plot_title = get(kwargs, :plot_title, "")
	title_string = "Occupancy, " * plot_title
	title(title_string)
	xlabel("Synthetic")
	ylabel("Physical")
	
	if get(kwargs, :if_save_fig, false)
		location = get(kwargs, :location, pwd())
		fig_name = get(kwargs, :name, "occs")
		fig_name = check_plot_label(fig_name,"occs")
		save_figure(fig_name; kwargs...)
	end
	return
end

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



# make pretty plots for OBC FQH droplet
function prettyplot_obc_occupancy(layers::Int64; kwargs...)
    intstren = get(kwargs, :onsite_strength, 0.0)
    particles = get(kwargs, :particles, Int(sqrt(2^layers)/2))
    hanis = get(kwargs, :hopping_anisotropy, 1.0)

    dataloc = get_folder_location("cluster-data/synth-dims/obc")
    params_dict = Dict([("layers",layers),("onsite_strength",intstren),("particles",particles),("hopping_anisotropy",hanis)])

    filename = find_data_file(params_dict,"ttn",dataloc)[1]
    d,m = read_data_jld2(filename,dataloc; output_level=0)

    occs = get_occupancy(d["densmat"]; if_plot=true)

    return filename
end

function prettyplot_obc_correlations(layers::Int64; kwargs...)
    intstren = get(kwargs, :onsite_strength, 0.0)
    particles = get(kwargs, :particles, Int(sqrt(2^layers)/2))
    hanis = get(kwargs, :hopping_anisotropy, 1.0)

    dataloc = get_folder_location("cluster-data/synth-dims/obc")
    params_dict = Dict([("layers",layers),("onsite_strength",intstren),("particles",particles),("hopping_anisotropy",hanis)])

    filename = find_data_file(params_dict,"ttn",dataloc)[1]
    d,m = read_data_jld2(filename,dataloc; output_level=0)

end























"fin"
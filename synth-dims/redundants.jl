#####################################################
#=

This file is for old redundant functions for TTNs and 1Deff

Depends on:
    

=#
######################################################

####### 1Deff section #######

function mps_plot_occupancy(occ_mat,L,nflavors; kwargs...)
	title_string = "Occupancy, " * get(kwargs, :plot_title, "")
	fig = figure()
	#plot_surface(1:nflavors,1:L,occ_mat)
	imshow(occ_mat)
	xlabel("Virtual Dim")
	ylabel("Physical Dim")
	colorbar()
	title(title_string)
	
	if_save_fig = get(kwargs, :if_save_fig, false)
	if if_save_fig
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "occs")
		filename = check_plot_label(filename,"occs")
	end
	if_save_fig ? save_figure(filename; location=location) : nothing
	return
end

function mps_plot_occupancy_3d(exp_occ; kwargs...)
	m,n = size(exp_occ)

	x = collect(1:m)  # x-coordinates (1 to m)
	y = collect(1:n)  # y-coordinates (1 to n)

	# Create the corresponding z-values from the matrix
	z = exp_occ

	# Now, we need to turn x, y, z into vectors for scatter3D
	# Use the function `repeat` and `vec` to flatten them
	x_flat = repeat(x, n)     # Repeat x n times
	y_flat = sort(repeat(y, m))    # Repeat each y m times
	z_flat = vec(z)           # Flatten the matrix into a vector

	# Create the 3D scatter plot
	title_string = "Occupancy, " * get(kwargs, :plot_title, "")
	fig = figure()
	scatter3D(x_flat, y_flat, z_flat)
	ylabel("Virtual Dim")
	xlabel("Physical Dim")
	title(title_string)
	return
end

function log_sum(all_values)
	if length(all_values) < 2
		return all_values[1]
	else
		consecutive = [all_values[1],all_values[2]]
		for i in 3:length(all_values) + 1
			added_value = log_add(consecutive[1],consecutive[2])
			consecutive[1] = added_value
			consecutive[2] = i <= length(all_values) ? all_values[i] : 0.0
		end
		return consecutive[1]
	end
end

# unsure about what this is for
function varied_alpha_wavefuncs(seed_wavefunc,seed_params_dict,change=0.001)
	og_alpha = seed_params_dict["phi"] / (2*pi)
	plus_psi = run_mps_new_variable(seed_wavefunc,seed_params_dict,Dict([("phi",2*pi*(og_alpha+change)),("psi_guess",seed_wavefunc)]))
	minus_psi = run_mps_new_variable(seed_wavefunc,seed_params_dict,Dict([("phi",2*pi*(og_alpha-change)),("psi_guess",seed_wavefunc)]))
	return minus_psi,plus_psi
end

#############################






















"fin"
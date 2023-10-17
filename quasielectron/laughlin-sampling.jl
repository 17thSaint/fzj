using LinearAlgebra,Statistics,PyPlot

function start_rand_config(num_parts::Int, p::Int)
    # Calculate the filling fraction
    filling = 1 / (2 * p + 1)

    # Calculate the characteristic length scale
    rm = sqrt(2 * num_parts / filling)

    # Generate random real and imaginary parts in one step
    real_parts = rand(Float64, num_parts) .* rand(-1:2:1, num_parts) .* rm
    imag_parts = rand(Float64, num_parts) .* rand(-1:2:1, num_parts) .* rm

    # Combine real and imaginary parts to create complex numbers
    config = real_parts .- im .* imag_parts

    return config
end

function particles_within_radius(particle_positions, center_point, radius)
    """
    Find the positions of particles within a given radius from a specified center particle (efficient version).

    Parameters:
        particle_positions: Array{Complex{Float64}} - An array of complex positions of particles.
        center_index: Int - The index of the center particle.
        radius: Float64 - The radius within which to find particles.

    Returns:
        Array{Complex{Float64}} - An array of complex positions of particles within the specified radius.
    """
    #center = particle_positions[center_index]
    return [particle for particle in particle_positions if abs(particle - center_point) <= radius]
end

function log_laughlin_wavefunction(z, m)
    """
    Compute the logarithm of the Laughlin wavefunction for a given set of particle positions.
    
    Parameters:
        z: Array{Complex{Float64}} - An array of complex positions of particles.
        m: Int - An integer parameter that determines the Laughlin state (e.g., 1 for Laughlin ν=1/3 state).

    Returns:
        Float64 - The natural logarithm of the Laughlin wavefunction at the given particle positions.
    """
    N = length(z)  # Number of particles
    log_norm_factor = 0.0
    for i in 1:N
        for j in 1:i-1
            log_norm_factor -= 2*m*log(Complex(abs(z[i]) - abs(z[j])))
        end
    end
    
    log_exponent = 0.5 * sum(abs2.(z))
    
    return log_norm_factor + log_exponent
end

function move_particle(num_parts::Int,chosen::Int,step_size::Float64)
	shift_matrix::Vector{ComplexF64} = [0.0+im*0.0 for i = 1:num_parts]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size - im*rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(initial_config::Array,chosen_particle::Int; kwargs...)
	approx_radius = get(kwargs, :approx_radius, 0.1*mean(abs.(initial_config)))
	step_size = get(kwargs, :step_size, 0.001*mean(abs.(initial_config)))
	m = get(kwargs, :m, 3)
	
	shift_matrix = move_particle(length(initial_config),chosen_particle,step_size)
	rand_num = rand(Float64)
	
	approx_start_config = particles_within_radius(initial_config, initial_config[chosen_particle], approx_radius)
	approx_new_config = particles_within_radius(initial_config + shift_matrix, initial_config[chosen_particle], approx_radius)
	
	start_wavefunc = log_laughlin_wavefunction(approx_start_config,m)
	new_wavefunc = log_laughlin_wavefunction(approx_new_config,m)
	
	start_ham = real(start_wavefunc)
	new_ham = real(new_wavefunc)
	check = new_ham - start_ham
	
	if isinf(new_ham) && isinf(start_ham)
		println("Both Inf")
		return true, intial_config, 0, start_wavefunc
	elseif isinf(new_ham)
		if new_ham > 0
			println("New Pos Inf -> Accept")
			return true, intial_config+shift_matrix, 1, new_wavefunc
		else
			#println("New Neg Inf -> Reject")
			return true, initial_config, 0, start_wavefunc
		end
	elseif isinf(start_ham)
		if start_ham > 0
			#println("Start Pos Inf -> Reject")
			return true, initial_config, 0, start_wavefunc
		else
			println("Start Neg Inf -> Accept")
			return true, initial_config+shift_matrix, 1, new_wavefunc
		end
	end
	
	if exp(-1*check) >= rand_num
		#println("Accept: ",new_ham,", ",start_ham,", ",rand_num)
		return true, initial_config+shift_matrix, 1, new_wavefunc
	else
		#println("Reject: ",new_ham,", ",start_ham,", ",rand_num)
		return true, initial_config, 0, start_wavefunc
	end
end

function mc(num_parts::Int,p::Int,steps::Int; kwargs...)
	steps_therm = get(kwargs, :therm_time, 2000)
	m = 2*p + 1
	running_config = start_rand_config(num_parts,p)
	start_count = 0
	index = 1
	starting_check = true
	
	while starting_check
		start_count += 1
		starting_wavefunc = log_laughlin_wavefunction(running_config, m)
		if isinf(real(starting_wavefunc))
			running_config = start_rand_config(num_parts,p)
		else
			#println("Started in $start_count steps")
			starting_check = false
		end
	end
	#
	for i_therm in 1:steps_therm
		for j_therm in 1:num_parts
			movement = acc_rej_move(running_config,j_therm; kwargs..., m = m)
			
			if movement[1]
				running_config = movement[2]
			else
				println("NaN Wavefunc")
				return movement
			end
			
		end
	end
	
	
	#println("Thermalization completed")
	#
	samp_freq = 10
	time_config = fill(0.0+im*0.0,(num_parts,Int(steps/samp_freq)))
	time_wavefunc = fill(0.0+im*0.0,(Int(steps/samp_freq)))
	acc_count = 0
	wavefunc = 0.0*im
	
	time_start = time()
	for i in 1:steps
		for j in 1:num_parts
			movement = acc_rej_move(running_config,j; kwargs..., m = m)
			
			if movement[1]
				acc_count += movement[3]
				running_config = movement[2]
				wavefunc = real(movement[4])
			else
				println("NaN Wavefunc")
				return movement
			end
			
		end
		
		acc_rate = acc_count/(num_parts*i)
		
		if i%samp_freq == 0
			time_config[:,index] = [Complex(running_config[x]) for x in 1:num_parts]
			local_config_time = time_config[:,index]
			time_wavefunc[index] = wavefunc
			index += 1
			#println("Added Data for Sampling Frequency")
		end
		#
		if i%(steps*0.01) == 0
			println("Running 1/$m:"," ",100*i/steps,"%: AccRate = ",acc_rate)
		end
	end
	
	time_end = time()
	total_time = (time_end - time_start)
	acc_rate = acc_count/(num_parts*steps)
	
	return acc_rate,time_config,time_wavefunc,total_time
end

function old_occupancy(time_config)
	num_particles = size(time_config)[1]
	compl = collect(Iterators.flatten([time_config[i,:] for i in 1:num_particles]))
	figure()
	hist2D(real.(compl),-imag.(compl),bins=100)
	colorbar()
end

function get_occupancy(time_config,rm,axis_bins=40)
	occs = zeros(axis_bins,axis_bins)
	max_x,min_x = maximum(real.(time_config))/rm,minimum(real.(time_config))/rm
	max_y,min_y = -1*minimum(imag.(time_config))/rm,-1*maximum(imag.(time_config))/rm
	bins_x = [min_x + (i-1)*(max_x-min_x)/(axis_bins-1) for i in 1:axis_bins]
	bins_y = [min_y + (i-1)*(max_y-min_y)/(axis_bins-1) for i in 1:axis_bins]
	for pos in time_config
		x_site = findfirst(i -> bins_x[i] < real(pos)/rm <= bins_x[i+1],[j for j in 1:axis_bins-1])
		y_site = findfirst(i -> bins_y[i] < -imag(pos)/rm <= bins_y[i+1],[j for j in 1:axis_bins-1])
		if isnothing(x_site) | isnothing(y_site)
			continue
		else
			occs[y_site,x_site] += 1/prod(size(time_config))
		end
	end
	fig = figure()
	imshow(occs)
	colorbar()
	return occs
end

function rad_dist(time_config,rm,axis_bins=40)
	occs = zeros(axis_bins)
	maxr = maximum(abs.(time_config))/rm
	bins = [-maxr + (i-1)*(2*maxr)/(axis_bins-1) for i in 1:axis_bins]
	raddist = zeros(axis_bins)
	for pos in time_config
		r_site = findfirst(i -> bins[i] < abs(pos)*sign(real(pos))/rm <= bins[i+1],[j for j in 1:axis_bins-1])
		if isnothing(r_site)
			continue
		else
			raddist[r_site] += 1
		end
	end
	fig = figure()
	plot(bins,raddist)
	xlabel("Radial Position")
	ylabel("Number Times Found")
	return bins,raddist
end

function radial_distribution(z,rm)
    """
    Calculate the radial distribution of particles from complex positions and create a histogram using PyPlot.

    Parameters:
        z: Array{Complex{Float64}} - An array of complex positions of particles.

    Returns:
        Tuple{Vector{Float64}, Vector{Int}} - A tuple containing two vectors:
            1. Vector of radial distances.
            2. Vector of counts in each radial bin.
    """
    radial_distances = [abs(particle)*sign(real(particle))/rm for particle in z]

    # Set up radial bins (adjust as needed)
    bin_edges = -maximum(radial_distances):0.01*rm:maximum(radial_distances)  # Define bin edges with a bin width of 0.1 units

    # Create a histogram using PyPlot
    plt.figure()
    counts, edges, _ = hist(radial_distances, bins=bin_edges, edgecolor="black")

    xlabel("Radial Distance")
    ylabel("Counts")
    title("Radial Distribution of Particles")

    return edges[2:end], counts
end

#
p = 1
particles = 10
mc_steps = 100000

rm = sqrt(2*particles*(2*p+1))
step_size = 0.1*rm
approx_rad = 2*step_size#0.2*rm
model_paras = (step_size = step_size, approx_radius = approx_rad, m = 2*p+1)
accrate,allconfigs,allpsis,runtime = mc(particles,p,mc_steps; model_paras...)
get_occupancy(allconfigs,rm)
rad_dist(allconfigs,rm)
#radial_distribution(collect(Iterators.flatten([allconfigs[i,:] for i in 1:particles])),rm)
#println("Runtime = ",runtime)





































"fin"

using LinearAlgebra,Statistics,PyPlot,SpecialMatrices

function start_rand_config(num_parts::Int, m::Int)
    # Calculate the filling fraction
    filling = 1 / m

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
            log_norm_factor += m*log(Complex(z[i] - z[j]))
        end
    end
    
    log_exponent = -0.25 * sum(abs2.(z))
    
    return log_norm_factor + log_exponent
end

function get_log_add(a,b)
	if real(a) > real(b)
		ordered::Vector{typeof(a)} = [b,a]
	else
		ordered = [a,b]
	end
	result::ComplexF64 = ordered[2] + log(Complex(1 + exp(ordered[1] - ordered[2])))
	return Complex(result)
end

function deriv_of_slater(config,which_particle=1)
	vdm_deriv = [(i-1)*config[which_particle]^(i-2) for i in 1:length(config)]
	return  -1 * det(Vandermonde(config)) * (( conj(Vandermonde(config)') \ vdm_deriv )[1])
end

function reverse_flux_wavefunction(z,m=3)
	p = Int((m+1)/2)
	num_parts = length(z)
	nt = sum([i for i in 1:num_parts-1])
	wavefunc = 0.0*im
	
	const_part = nt * (log(4*p) + m*log(Complex(det(Vandermonde(z)))))
	
	deriv_part = 0.0*im
	all_derivs = [deriv_of_slater(z,i) for i in 1:num_parts]
	for i in 1:num_parts
		for j in 1:i-1
			new_term = all_derivs[j] - all_derivs[i]
			deriv_part += log(Complex(new_term))
		end
	end
	
	log_exponent = -0.25 * sum(abs2.(z))
	
	return log_exponent + deriv_part + const_part
end

function move_particle(num_parts::Int,chosen::Int,step_size::Float64)
	shift_matrix::Vector{ComplexF64} = [0.0+im*0.0 for i = 1:num_parts]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size - im*rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(initial_config::Array,chosen_particle::Int; kwargs...)
	#approx_radius = get(kwargs, :approx_radius, 0.1*mean(abs.(initial_config)))
	step_size = get(kwargs, :step_size, 0.001*mean(abs.(initial_config)))
	m = get(kwargs, :m, 3)
	start_wavefunc = get(kwargs, :start_wavefunc, nothing)
	wavefunc_type = get(kwargs, :vers, "P")
	wave_function = wavefunc_type == "P" ? log_laughlin_wavefunction : reverse_flux_wavefunction
	
	shift_matrix = move_particle(length(initial_config),chosen_particle,step_size)
	rand_num = log(rand(Float64))
	
	#approx_start_config = particles_within_radius(initial_config, initial_config[chosen_particle], approx_radius)
	#approx_new_config = particles_within_radius(initial_config + shift_matrix, initial_config[chosen_particle], approx_radius)
	
	if isnothing(start_wavefunc)
		start_wavefunc = wave_function(initial_config,m)
		println("Making own wavefunc")
	end
	new_wavefunc = wave_function(initial_config+shift_matrix,m)
	
	start_ham = 2*real(start_wavefunc)
	new_ham = 2*real(new_wavefunc)
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
	
	if check >= rand_num
		#println("Accept: ",new_ham,", ",start_ham,", ",rand_num)
		return true, initial_config+shift_matrix, 1, new_wavefunc
	else
		#println("Reject: ",new_ham,", ",start_ham,", ",rand_num)
		return true, initial_config, 0, start_wavefunc
	end
end

function mc(num_parts::Int,m::Int,steps::Int; kwargs...)
	steps_therm = get(kwargs, :therm_time, 2000)
	opl = get(kwargs, :opl, 1)
	samp_freq = get(kwargs, :samp_freq, 10)
	wavefunc_type = get(kwargs, :vers, "P")
	wave_function = wavefunc_type == "P" ? log_laughlin_wavefunction : reverse_flux_wavefunction
	
	running_config = start_rand_config(num_parts,m)
	start_count = 0
	index = 1
	starting_check = true
	wavefunc = 0.0*im
	
	while starting_check
		start_count += 1
		wavefunc = wave_function(running_config, m)
		if isinf(real(wavefunc))
			running_config = start_rand_config(num_parts,m)
		else
			#println("Started in $start_count steps")
			starting_check = false
		end
	end
	#
	for i_therm in 1:steps_therm
		for j_therm in 1:num_parts
			movement = acc_rej_move(running_config,j_therm; kwargs..., m = m,start_wavefunc=wavefunc)
			
			if movement[1]
				running_config = movement[2]
				wavefunc = movement[4]
			else
				println("NaN Wavefunc")
				return movement
			end
			
		end
	end
	
	
	#println("Thermalization completed")
	#
	time_config = fill(0.0+im*0.0,(num_parts,Int(steps/samp_freq)))
	time_wavefunc = fill(0.0+im*0.0,(Int(steps/samp_freq)))
	acc_count = 0
	
	time_start = time()
	for i in 1:steps
		for j in 1:num_parts
			movement = acc_rej_move(running_config,j; kwargs..., m = m,start_wavefunc=wavefunc)
			
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
		if opl == 1 && i%(steps*0.01) == 0
			println("Running 1/$m:"," ",100*i/steps,"%: AccRate = ",acc_rate)
		end
	end
	
	time_end = time()
	total_time = (time_end - time_start)
	acc_rate = acc_count/(num_parts*steps)
	
	return acc_rate,time_config,time_wavefunc,total_time
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

function rad_dist(time_config,rm; kwargs...)
	axis_bins = get(kwargs, :axis_bins, 40)
	maxr = get(kwargs, :maxr, maximum(abs.(time_config))/rm)
	title_string = get(kwargs, :titlestring, "")
	label_string = get(kwargs, :labelstring, "")
	
	occs = zeros(axis_bins)
	#bins = [-maxr + (i-1)*(2*maxr)/(axis_bins-1) for i in 1:axis_bins]
	bins = [0.0 + (i-1)*(maxr)/(axis_bins-1) for i in 1:axis_bins]
	raddist = zeros(axis_bins)
	for pos in time_config
		#r_site = findfirst(i -> bins[i] < abs(pos)*sign(real(pos))/rm <= bins[i+1],[j for j in 1:axis_bins-1])
		r_site = findfirst(i -> bins[i] < abs(pos)/rm <= bins[i+1],[j for j in 1:axis_bins-1])
		if isnothing(r_site)
			continue
		else
			raddist[r_site] += 1/size(time_config)[1]
		end
	end
	#fig = figure()
	plot(bins,raddist,label=label_string)
	xlabel("Radial Position")
	ylabel("Number Times Found")
	title(title_string)
	legend()
	return bins,raddist
end

function auto_correlation(energies, delta_t)
    average_energy = mean(energies)
    
    points = Int(floor(length(energies)-delta_t))
    
    energy_fluctuations = energies.-average_energy
    
    autocorrelation_top = [0.0 for i in 1:points]
    autocorrelation_bottom = mean(energy_fluctuations.^2)
    
    for i in 1:points
        autocorrelation_top[i] = energy_fluctuations[i]*energy_fluctuations[i+delta_t]
    end
    
    return (mean(autocorrelation_top)/autocorrelation_bottom)
end

function get_autocorr_length(wavefunc_data,samp_freq)
	energy = real.(wavefunc_data)
	full_length = length(wavefunc_data)
	#println("Full Length = $full_length")
	len = full_length
	dts = [1+(i-1)*1 for i in 1:Int(0.1*len)-1]
	autocorr = [0.0 for i in 1:Int(0.1*len)-1]
	for i in 1:Int(0.1*len)-1
		autocorr[i] = auto_correlation(energy,dts[i])
	end
	check_below_tol = autocorr .< [0.01 for i in 1:length(autocorr)]
	#println(autocorr)
	corr_length = samp_freq*dts[findall(check_below_tol)[1]]
	return corr_length,dts,autocorr
end

function get_rightband_count(edge,pos_data,rm)
	local_count = 0
	for pos in pos_data
		if -edge <= real(pos) <= edge && abs(imag(pos)) <= 0.1*rm 
			local_count += 1/prod(size(pos_data))
		end
	end
	return local_count
end

function radial_density_full(pos_data,rm; kwargs...)
	points = get(kwargs, :points, 100)
	rend = get(kwargs, :rend, 1.25*rm)
	typeof(rend) == String ? rend = maximum(abs.(pos_data)) : nothing
	if_plot = get(kwargs, :if_plot, true)
	title_string = get(kwargs, :titlestring, "")
	label_string = get(kwargs, :labelstring, "")
	
	allxs = [(i-1)*rend/(points-1) for i in 1:points]
	raddens = [0.0 for i in 1:points]
	for (i,ed) in enumerate(allxs)
		if i == 1
			shift = 0
		else
			shift = sum(raddens[1:i-1])
		end
		raddens[i] = get_rightband_count(ed,pos_data,rm) - shift
	end
	
	if if_plot
		fig = figure()
		plot(allxs ./ rm,raddens ./ sum(raddens),label=label_string)
		title(title_string)
		xlabel("X Value / rm")
		ylabel("Particle Density")
	end
	
	return allxs,raddens
end


#
m = 3
mc_steps = 1000000
output = 1
sampfreq = 1
particles = 10
rm = sqrt(2*particles*m)
step_size = 0.125*rm
ver = "R"

model_paras = (vers = ver, step_size = step_size, m = m, opl = output, samp_freq = sampfreq)
accrate,allconfigs,allpsis,runtime = mc(particles,m,mc_steps; model_paras...)

raddens = radial_density_full(allconfigs,rm; rend="max")
#rad_dist(allconfigs,rm; axis_bins=100, maxr = 1.25, labelstring = "$(round(approx_rad/rm,digits=2))")
get_occupancy(allconfigs,rm)
#title("$(approx_rad/rm)")
#end
#raddata = rad_dist(allconfigs,rm; axis_bins=100, labelstring = "$particles")
#end
#println("Runtime = ",runtime)
#




































"fin"

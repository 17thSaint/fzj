using LinearAlgebra,Statistics,PyPlot
include("analysis-functions.jl")
include("reverse-flux.jl")
include("laughlin-wavefunc.jl")

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

function get_log_add(a,b)
	if real(a) > real(b)
		ordered::Vector{typeof(a)} = [b,a]
	else
		ordered = [a,b]
	end
	result::ComplexF64 = ordered[2] + log(Complex(1 + exp(ordered[1] - ordered[2])))
	return Complex(result)
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
		start_wavefunc,timesignore = wave_function(initial_config,m)
		println("Making own wavefunc")
	end
	
	new_wavefunc,times = wave_function(initial_config+shift_matrix,m)
	
	check = 0.0
	check_infinity = true
	infinity_steps = 0
	while check_infinity
		
		start_ham = 2*real(start_wavefunc)
		new_ham = 2*real(new_wavefunc)
		check = new_ham - start_ham
		
		#=if isinf(new_ham) && isinf(start_ham)
			println("Both Inf")
			return true, initial_config, 0, start_wavefunc, times
		=#
		if isinf(new_ham)
			#println("New Inf")
			shift_matrix = move_particle(length(initial_config),chosen_particle,step_size)
			new_wavefunc,times = wave_function(initial_config+shift_matrix,m)
			infinity_steps += 1
		#=elseif isinf(start_ham)
			if start_ham > 0
				println("Start Pos Inf -> Reject")
				return true, initial_config, 0, start_wavefunc, times
			else
				println("Start Neg Inf -> Accept")
				return true, initial_config+shift_matrix, 1, new_wavefunc, times
			end
		=#
		elseif !isinf(new_ham) && !isinf(start_ham)
			#=if infinity_steps > 0
				println("Found Finite Config in $infinity_steps Steps")
			end
			=#
			check_infinity = false
		end
	end
	
	if check >= rand_num
		#println("Accept: ",new_ham,", ",start_ham,", ",rand_num)
		return true, initial_config+shift_matrix, 1, new_wavefunc, times
	else
		#println("Reject: ",new_ham,", ",start_ham,", ",rand_num)
		return true, initial_config, 0, start_wavefunc, times
	end
end

function mc(num_parts::Int,m::Int,steps::Int; kwargs...)
	steps_therm = get(kwargs, :therm_time, 1000)
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
		wavefunc = wave_function(running_config, m)[1]
		if isinf(real(wavefunc))
			running_config = start_rand_config(num_parts,m)
		else
			println("Started in $start_count steps")
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
	
	
	println("Thermalization completed")
	#
	time_config = fill(0.0+im*0.0,(num_parts,Int(steps/samp_freq)))
	time_wavefunc = fill(0.0+im*0.0,(Int(steps/samp_freq)))
	acc_count = 0
	
	all_times = [[],[]]
	
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
			
			append!(all_times[1],[movement[5][1]])
			append!(all_times[2],[movement[5][2]])
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
			println("Running 1/$m with $num_parts particles:"," ",100*i/steps,"%: AccRate = ",acc_rate)
		end
	end
	
	time_end = time()
	total_time = (time_end - time_start)
	acc_rate = acc_count/(num_parts*steps)
	
	const_mean  = round(mean(all_times[1]),digits=4)
	const_std  = round(std(all_times[1]),digits=4)
	deriv_mean  = round(mean(all_times[2]),digits=4)
	deriv_std  = round(std(all_times[2]),digits=4)
	println("Times: Const = $const_mean +- $const_std, Deriv = $deriv_mean +- $deriv_std")
	
	return acc_rate,time_config,time_wavefunc,total_time
end




#
m = 3
mc_steps = 100000
output = 1
sampfreq = 1
#everyconfig = []
allraddens = []
for particles in [i for i in 3:9]
#particles = 6
rm = sqrt(2*particles*m)
step_size = 0.125*rm
ver = "R"

model_paras = (vers = ver, step_size = step_size, m = m, opl = output, samp_freq = sampfreq)
#accrate,allconfigs,allpsis,runtime = mc(particles,m,mc_steps; model_paras...)
#append!(everyconfig,[allconfigs])
allconfigs = everyconfig[particles-2]
raddists = rad_dist(allconfigs,rm; axis_bins=100,labelstring="$particles")
#raddens = radial_density_full(allconfigs,rm; rend="max", labelstring = "$particles")
append!(allraddens,[raddists])
end
#rad_dist(allconfigs,rm; axis_bins=100, maxr = 1.25, labelstring = "$(round(approx_rad/rm,digits=2))")
#get_occupancy(allconfigs,rm)
#title("$(approx_rad/rm)")
#end
#raddata = rad_dist(allconfigs,rm; axis_bins=100, labelstring = "$particles")
#end
#println("Runtime = ",runtime)
#




































"fin"

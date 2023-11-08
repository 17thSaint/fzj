using LinearAlgebra,Statistics,PyPlot,LsqFit
include("analysis-functions.jl")
include("laughlin-wavefunc.jl")
include("../other-funcs/data-storage-funcs.jl")

#include("fqh-thesis/cf-wavefunc.jl")

include("reverse-flux.jl")

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
	wave_function = wavefunc_type == "P" ? laughlin_wavefunction_girvinjach : reverse_flux_wavefunction
	
	shift_matrix = move_particle(length(initial_config),chosen_particle,step_size)
	rand_num = log(rand(Float64))
	
	#approx_start_config = particles_within_radius(initial_config, initial_config[chosen_particle], approx_radius)
	#approx_new_config = particles_within_radius(initial_config + shift_matrix, initial_config[chosen_particle], approx_radius)
	
	if isnothing(start_wavefunc)	
		if wavefunc_type == "old-rf"
			times = [0,0]
			acc_sets_matrix = get(kwargs, :accmat, nothing)
			all_pascal = get(kwargs, :pascal, nothing)
			all_deriv_orders = get(kwargs, :derivs, nothing)
			new_wavefunc = get_rf_wavefunc(initial_config,acc_sets_matrix,all_pascal,all_deriv_orders,[0,[0]],true)
		else
			start_wavefunc,timesignore = wave_function(initial_config,m; kwargs...)
		end
		println("Making own wavefunc")
	end
	
	if wavefunc_type == "old-rf"
		times = [0,0]
		acc_sets_matrix = get(kwargs, :accmat, nothing)
		all_pascal = get(kwargs, :pascal, nothing)
		all_deriv_orders = get(kwargs, :derivs, nothing)
		new_wavefunc = get_rf_wavefunc(initial_config+shift_matrix,acc_sets_matrix,all_pascal,all_deriv_orders,[0,[0]],true)
	else
		new_wavefunc,times = wave_function(initial_config+shift_matrix,m; kwargs...)
	end
	
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
			println("New Inf")
			shift_matrix = move_particle(length(initial_config),chosen_particle,step_size)
			new_wavefunc,times = wave_function(initial_config+shift_matrix,m; kwargs...)
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
		#println("Reject: ",check,", ",rand_num)
		return true, initial_config, 0, start_wavefunc, times
	end
end

function mc(num_parts::Int,m::Int,steps::Int; kwargs...)
	steps_therm = get(kwargs, :therm_time, 1000)
	opl = get(kwargs, :opl, 1)
	samp_freq = get(kwargs, :samp_freq, 10)
	if_save_data = get(kwargs, :if_save_data, false)
	wavefunc_type = get(kwargs, :vers, "P")
	wave_function = wavefunc_type == "P" ? laughlin_wavefunction_girvinjach : reverse_flux_wavefunction
	qe_cutoff = get(kwargs, :qe_cutoff, 0)
	qe_loc = get(kwargs, :qe_loc, 0.0+im*0.0)
	
	rm = sqrt(2*num_parts*m)
	
	if qe_cutoff != 0
		println("Using Quasielectron with cutoff at $qe_cutoff / $num_parts")
	end
	
	running_config = start_rand_config(num_parts,m)
	start_count = 0
	index = 1
	starting_check = true
	wavefunc = 0.0*im
	
	#
	while starting_check
		start_count += 1
		if wavefunc_type == "old-rf"
			acc_sets_matrix = get(kwargs, :accmat, nothing)
			all_pascal = get(kwargs, :pascal, nothing)
			all_deriv_orders = get(kwargs, :derivs, nothing)
			wavefunc = get_rf_wavefunc(running_config,acc_sets_matrix,all_pascal,all_deriv_orders,[0,[0]],true)
		else
			wavefunc = wave_function(running_config,m; kwargs...)[1]
		end
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
			rejected = true
			attempts = 0
			while rejected
				movement = acc_rej_move(running_config,j_therm; kwargs..., m = m,start_wavefunc=wavefunc)
				if movement[3] == 1
					running_config = movement[2]
					wavefunc = movement[4]
					rejected = false
				else
					attempts += 1
				end
				if wavefunc_type == "R"
					if attempts > 50
						println("Too many attempts to move during Thermalization")
						return nothing
					end
				end
			end
			#=
			if movement[1]
				running_config = movement[2]
				wavefunc = movement[4]
			else
				println("NaN Wavefunc")
				return movement
			end
			=#
			
		end
	end
	
	
	println("Thermalization completed")
	#
	time_config = fill(0.0+im*0.0,(num_parts,Int(steps/samp_freq)))
	time_wavefunc = fill(0.0+im*0.0,(Int(steps/samp_freq)))
	
	all_times = [[],[]]
	
	i = 0
	steps_done = 0
	time_start = time()
	for i in 1:steps
		rej_count = 0
		for j in 1:num_parts
			rejected = true
			attempts = 0
			while rejected
				movement = acc_rej_move(running_config,j; kwargs..., m = m,start_wavefunc=wavefunc)
				if movement[3] == 1
					running_config = movement[2]
					wavefunc = real(movement[4])
					if wavefunc_type == "R"
						append!(all_times[1],[movement[5][1]])
						append!(all_times[2],[movement[5][2]])
					end
					rejected = false
				else
					attempts += 1
				end
				if wavefunc_type == "R"
					if attempts > 50
						println("Too many attempts to move")
						return nothing
					end
				end
			end
			
			#=
			if movement[1]
				rej_count += attempts
				running_config = movement[2]
				wavefunc = real(movement[4])
			else
				println("NaN Wavefunc")
				return movement
			end
			=#

			rej_count += attempts
						
		end
		
		acc_rate = num_parts/(num_parts+rej_count)
		
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
	#acc_rate = acc_count/(num_parts*steps)
	
	if wavefunc_type == "R"
		const_mean  = round(mean(all_times[1]),digits=4)
		const_std  = round(std(all_times[1]),digits=4)
		deriv_mean  = round(mean(all_times[2]),digits=4)
		deriv_std  = round(std(all_times[2]),digits=4)
		println("Times: Const = $const_mean +- $const_std, Deriv = $deriv_mean +- $deriv_std")
	end
	
	if if_save_data
		location = get(kwargs, :location, "../cluster-data/quasielectron")
		params_dict::Dict{String,Any} = Dict([("mc_steps",steps),("m",m),("particles",num_parts)])
		if qe_cutoff > 0
			params_dict["qe_loc"] = round(real(qe_loc)/rm,digits=4)
			params_dict["qe_cutoff"] = qe_cutoff
		end
		wavefunc_type_name = wavefunc_type == "P" ? "laugh-" : "rfa-"
		filename = get(kwargs, :name, wavefunc_type_name * make_parameters_filename(params_dict))
		metadata_dict = merge(named_tuple_to_dict(kwargs),params_dict)
		metadata = get(kwargs, :metadata, metadata_dict)
		config_data_dict = Dict([("configs",time_config),("wavefuncs",time_wavefunc)])
		println(filename)
		display(metadata)
		write_data_jld2(filename,config_data_dict,location,metadata)
	end
	
	return time_config,time_wavefunc,total_time
end



#
axisbins = 300
m = 1
mc_steps = 1000000
output = 1
sampfreq = 1
num_parts = [i for i in 8:8]
for particles in num_parts
#particles = 4
rm = sqrt(2*particles*m)
step_size = 0.125*rm
ver = "P"

#=
allowed_sets_matrix = get_full_acc_matrix(particles)
full_pasc_tri = [get_pascals_triangle(i)[2] for i in 1:particles]
full_derivs = get_deriv_orders_matrix(particles)
oldRF_paras = (accmat = allowed_sets_matrix, pascal = full_pasc_tri, derivs = full_derivs)
=#

for i in 1:10
model_paras = (therm_time = 1000, vers = ver, step_size = step_size, m = m, opl = output, samp_freq = sampfreq, if_save_data = true, qe_cutoff = particles)
allconfigs,allpsis,runtime = mc(particles,m,mc_steps; model_paras...,qe_loc = rm*i/10)
end

#allconfigs = everyconfig[particles-12]
#raddenss = radial_density_full(allconfigs,rm; points=axisbins, rend="max",if_plot=true,labelstring="with QE")
#corrlength = get_autocorr_length(allpsis,1)[1]
#println("CorrLength = ",corrlength)
#get_occupancy(allconfigs,rm; title_string="N = $particles")
#edge_loc = raddenss[1][findfirst(i -> raddenss[2][i] == maximum(raddenss[2]),1:axisbins)]
#append!(edges,[edge_loc])
#println("Edge = ",edge_loc)
#plot([i for i in 5:particles],edges,"-p",c="b")
#=raddists = rad_dist(allconfigs,rm; axis_bins=axisbins,labelstring="$particles")

if particles > 5
	start_point = findfirst(i->raddists[2][i] > 0.2,Int(axisbins/2):axisbins) + Int(axisbins/2)
	end_point = findfirst(i->raddists[2][i] < 0.2,start_point+10:axisbins) + start_point + 10
	println(raddists[1][start_point],", ",raddists[1][end_point])
	maxval = findfirst(i->raddists[2][i]==maximum(raddists[2][start_point:end_point]),start_point:end_point) + start_point
	p0 = [0.1,raddists[1][maxval],0.5]
	localfit = LsqFit.curve_fit(gauss,raddists[1][start_point:end_point],raddists[2][start_point:end_point],p0)
	append!(allfits,[localfit.param])
	scatter(raddists[1][start_point:end_point],gauss(raddists[1][start_point:end_point],localfit.param))
else
	start_point = findfirst(i->raddists[2][i] >= 0.2,1:axisbins)
	end_point = findfirst(i->raddists[2][i] <= 0.2,start_point+10:axisbins) + start_point+10
	println(raddists[1][start_point],", ",raddists[1][end_point])
	maxval = findfirst(i->raddists[2][i]==maximum(raddists[2][start_point:end_point]),start_point:end_point) + start_point
	p0 = [0.5,raddists[1][maxval],0.5]
	localfit = LsqFit.curve_fit(gauss,raddists[1][start_point:end_point],raddists[2][start_point:end_point],p0)
	append!(allfits,[localfit.param])
	scatter(raddists[1][start_point:end_point],gauss(raddists[1][start_point:end_point],localfit.param))
end

if particles > 5
	start_point = findfirst(i->raddens[2][i] > 0.2,Int(axisbins/2):axisbins) + Int(axisbins/2)
	end_point = findfirst(i->raddens[2][i] < 0.2,start_point+10:axisbins) + start_point + 10
	println(raddens[1][start_point],", ",raddens[1][end_point])
	maxval = findfirst(i->raddens[2][i]==maximum(raddens[2][start_point:end_point]),start_point:end_point) + start_point
	p0 = [0.1,raddens[1][maxval],0.5]
	localfit = LsqFit.curve_fit(gauss,raddens[1][start_point:end_point],raddens[2][start_point:end_point],p0)
	append!(allfits,[localfit.param])
	scatter(raddens[1][start_point:end_point],gauss(raddens[1][start_point:end_point],localfit.param))
else
	start_point = findfirst(i->raddens[2][i] >= 0.2,1:axisbins)
	end_point = findfirst(i->raddens[2][i] <= 0.2,start_point+10:axisbins) + start_point+10
	println(raddens[1][start_point],", ",raddens[1][end_point])
	maxval = findfirst(i->raddens[2][i]==maximum(raddens[2][start_point:end_point]),start_point:end_point) + start_point
	p0 = [0.5,raddens[1][maxval],0.5]
	localfit = LsqFit.curve_fit(gauss,raddens[1][start_point:end_point],raddens[2][start_point:end_point],p0)
	append!(allfits,[localfit.param])
	scatter(raddens[1][start_point:end_point],gauss(raddens[1][start_point:end_point],localfit.param))
end
=#
#append!(allraddens,[raddists])
#
#end
#rad_dist(allconfigs,rm; axis_bins=100, maxr = 1.25, labelstring = "$(round(approx_rad/rm,digits=2))")
#get_occupancy(allconfigs,rm)
#title("$(approx_rad/rm)")
#end
#raddata = rad_dist(allconfigs,rm; axis_bins=100, labelstring = "$particles")
end
#
#println("Runtime = ",runtime)
#
#=fig = figure()
edges = [0.0 for i in 1:length(num_parts)]
for i in 1:length(num_parts)
	parts = num_parts[i]
	rm = sqrt(2*parts*3)
	
	con = everyconfig[i]
	#get_occupancy(con,rm,60; title_string="N = $parts")
	
	xs,rads = radial_density_full(con,rm; rend="max", labelstring="$parts", if_plot=false)
	edge_loc = xs[end]#xs[findfirst(j -> rads[j] == maximum(rads),[k for k in 1:length(rads)])]
	edges[i] = edge_loc
end
#
fig = figure()
plot(num_parts,edges,"-p")
xlabel("Particle Number")
ylabel("Radius of Edge / rm")
title("RF Wavefunction")
=#

#= Location of edges for Laughlin
# for particles = [i for i in 5:20]
partsss1 = [i for i in 5:20]
edges1 = [1.495859155858919, 1.435630612727625, 1.37945930696872, 1.3654778777900691, 1.3864465233636276, 1.3130718120943006, 1.3447262795261632, 1.3038602437553868, 1.3199477842162546, 1.283222341188332, 1.2888736012325286, 1.2742575476804927, 1.2755683525810448, 1.2522972698797814, 1.2562847066293055, 1.2856352734705268]
# for particles = [i for i in 10:30] with 10000 mc_steps
partsss2 = [i for i in 10:30]
edges_2 = [1.268621992836878, 1.323711571746685, 1.2754580172717975, 1.3377087481957572, 1.2428266925375056, 1.2268445567658224, 1.2688307138262522, 1.2453171121437712, 1.2292019961086964, 1.2603686230278608, 1.2667532728302995, 1.2025124522874775, 1.2418971614954482, 1.2010335988960956, 1.1871327744952063, 1.195154274324511, 1.2066232617157169, 1.1958074532675584, 1.2174605053706915, 1.1921377687610828, 1.21118076453031]

plot(partsss2[11:end],edges_2[11:end],"-p",c="b")
plot(partsss1,edges1,"-p",c="b")
xlabel("Particle Number")
ylabel("Radius of Edge / rm")
title("Laughlin Wavefunction")
=#



























"fin"

using NumericalIntegration

function plot_circle(radius,center)
	counts = 100
	xs = [-(radius-10^-4) + (i-1)*2*(radius-10^-4)/(counts-1) for i in 1:counts]
	ys = [sqrt(radius^2 - x^2) for x in xs]
	plot(xs .+ center,ys.+ center,c="r")
	plot(xs .+ center,-1 .* ys .+ center,c="r")
end

function get_occupancy(time_config,rm,axis_bins=40; kwargs...)
	title_string = get(kwargs, :title_string, "")
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
	binwidth = abs(bins_x[2] - bins_x[1])
	fig = figure()
	imshow(occs)
	colorbar()
	plot_circle(1/binwidth,(axis_bins/2)-0.5)
	title(title_string)
	#plot_circle(5*rm,axis_bins/2)
	return occs
end

function rad_dist(time_config,rm; kwargs...)
	axis_bins = get(kwargs, :axis_bins, 40)
	maxr = get(kwargs, :maxr, maximum(abs.(time_config))/rm)
	title_string = get(kwargs, :titlestring, "")
	label_string = get(kwargs, :labelstring, "")
	if_plot = get(kwargs, :if_plot, true)
	
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
	normalization = integrate(bins, raddist)
	if if_plot
		#fig = figure()
		plot(bins,raddist ./ normalization,label=label_string)
		xlabel("Radial Position / rm")
		ylabel("Number Times Found")
		title(title_string)
		legend()
	end
	return bins,raddist ./ normalization
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
		normalization = integrate(allxs ./ rm, raddens)
		#fig = figure()
		plot(allxs ./ rm,raddens ./ normalization,label=label_string)
		title(title_string)
		xlabel("X Value / rm")
		ylabel("Particle Density")
		legend()
	end
	
	return allxs ./ rm,raddens ./ normalization
end




































"fin"

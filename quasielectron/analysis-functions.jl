using NumericalIntegration,LsqFit,PyPlot
include("../other-funcs/data-storage-funcs.jl")

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
	
	normalization = integrate(allxs ./ rm, raddens)
	if if_plot
		fig = figure()
		plot(allxs ./ rm,raddens ./ normalization,label=label_string)
		title(title_string)
		xlabel("X Value / rm")
		ylabel("Particle Density")
		legend()
	end
	
	return allxs ./ rm,raddens ./ normalization
end

function fit_gaussian_raddens(particles,raddens,axisbins; kwargs...)
	if_plot = get(kwargs, :if_plot, true)
	gauss(x,p) = (1/(p[1]*sqrt(2*pi))) .* exp.(-0.5 .* (((x .- p[2]) ./ p[1]).^2)) .+ p[3]
	yaxis_min = 0.1
	if 5 < particles < 17
		fit_parameters = [[],[]]
		
		start_point = findfirst(i->raddens[2][i] > yaxis_min,Int(axisbins/2):axisbins) + Int(axisbins/2)
		end_point = findfirst(i->raddens[2][i] < yaxis_min,start_point+10:axisbins) + start_point + 10
		#println(raddens[1][start_point],", ",raddens[1][end_point])
		maxval = findfirst(i->raddens[2][i]==maximum(raddens[2][start_point:end_point]),start_point:end_point) + start_point
		p0 = [0.1,raddens[1][maxval],0.5]
		localfit = LsqFit.curve_fit(gauss,raddens[1][start_point:end_point],raddens[2][start_point:end_point],p0)
		#append!(allfits,[localfit.param])
		if_plot ? scatter(raddens[1][start_point:end_point],gauss(raddens[1][start_point:end_point],localfit.param)) : nothing
		fit_parameters[2] = localfit.param
		
		start_point_left = findfirst(i->raddens[2][i] > yaxis_min,1:Int(axisbins/2))
		end_point_left = findfirst(i->raddens[2][i] < yaxis_min,start_point_left+10:Int(axisbins/2)) + start_point_left + 10
		#println(raddens[1][start_point],", ",raddens[1][end_point])
		maxval_left = findfirst(i->raddens[2][i]==maximum(raddens[2][start_point_left:end_point_left]),start_point_left:end_point_left) + start_point_left
		p0_left = [0.1,raddens[1][maxval_left],0.5]
		localfit_left = LsqFit.curve_fit(gauss,raddens[1][start_point_left:end_point_left],raddens[2][start_point_left:end_point_left],p0_left)
		#append!(allfits,[localfit.param])
		if_plot ? scatter(raddens[1][start_point_left:end_point_left],gauss(raddens[1][start_point_left:end_point_left],localfit_left.param)) : nothing
		fit_parameters[1] = localfit_left.param
	elseif particles > 16
		fit_parameters = [[],[],[]]
		
		start_point_right = findfirst(i->raddens[2][i] > yaxis_min,Int(2*axisbins/3):axisbins) + Int(2*axisbins/3)
		end_point_right = findfirst(i->raddens[2][i] < yaxis_min,start_point_right+10:axisbins) + start_point_right + 10
		#println(raddens[1][start_point],", ",raddens[1][end_point])
		maxval_right = findfirst(i->raddens[2][i]==maximum(raddens[2][start_point_right:end_point_right]),start_point_right:end_point_right) + start_point_right
		p0_right = [0.1,raddens[1][maxval_right],0.5]
		localfit_right = LsqFit.curve_fit(gauss,raddens[1][start_point_right:end_point_right],raddens[2][start_point_right:end_point_right],p0_right)
		#append!(allfits,[localfit.param])
		if_plot ? scatter(raddens[1][start_point_right:end_point_right],gauss(raddens[1][start_point_right:end_point_right],localfit_right.param)) : nothing
		fit_parameters[3] = localfit_right.param
		
		start_point_center = findfirst(i->raddens[2][i] > yaxis_min,Int(axisbins/3):Int(2*axisbins/3)) + Int(1*axisbins/3)
		end_point_center = findfirst(i->raddens[2][i] < yaxis_min,start_point_center+10:axisbins) + start_point_center + 10
		#println(raddens[1][start_point],", ",raddens[1][end_point])
		maxval_center = findfirst(i->raddens[2][i]==maximum(raddens[2][start_point_center:end_point_center]),start_point_center:end_point_center) + start_point_center
		p0_center = [0.1,raddens[1][maxval_center],0.5]
		localfit_center = LsqFit.curve_fit(gauss,raddens[1][start_point_center:end_point_center],raddens[2][start_point_center:end_point_center],p0_center)
		#append!(allfits,[localfit.param])
		if_plot ? scatter(raddens[1][start_point_center:end_point_center],gauss(raddens[1][start_point_center:end_point_center],localfit_center.param)) : nothing
		fit_parameters[2] = localfit_center.param
		
		start_point_left = findfirst(i->raddens[2][i] > yaxis_min,1:Int(axisbins/3))
		end_point_left = findfirst(i->raddens[2][i] < yaxis_min,start_point_left+10:Int(axisbins/3)) + start_point_left + 10
		#println(raddens[1][start_point],", ",raddens[1][end_point])
		maxval_left = findfirst(i->raddens[2][i]==maximum(raddens[2][start_point_left:end_point_left]),start_point_left:end_point_left) + start_point_left
		p0_left = [0.1,raddens[1][maxval_left],0.5]
		localfit_left = LsqFit.curve_fit(gauss,raddens[1][start_point_left:end_point_left],raddens[2][start_point_left:end_point_left],p0_left)
		#append!(allfits,[localfit.param])
		if_plot ? scatter(raddens[1][start_point_left:end_point_left],gauss(raddens[1][start_point_left:end_point_left],localfit_left.param)) : nothing
		fit_parameters[1] = localfit_left.param
	else
		start_point = findfirst(i->raddens[2][i] >= yaxis_min,1:axisbins)
		end_point = findfirst(i->raddens[2][i] <= yaxis_min,start_point+10:axisbins) + start_point+10
		#println(raddens[1][start_point],", ",raddens[1][end_point])
		maxval = findfirst(i->raddens[2][i]==maximum(raddens[2][start_point:end_point]),start_point:end_point) + start_point
		p0 = [0.5,raddens[1][maxval],0.5]
		localfit = LsqFit.curve_fit(gauss,raddens[1][start_point:end_point],raddens[2][start_point:end_point],p0)
		#append!(allfits,[localfit.param])
		if_plot ? scatter(raddens[1][start_point:end_point],gauss(raddens[1][start_point:end_point],localfit.param)) : nothing
		fit_parameters = localfit.param
	end
	return fit_parameters
end

#=
dataloc = "/home/patrick/fzj/main-git/cluster-data/quasielectron"
paradict = Dict([("m",3),("mc_steps",100000)])
allfiles = find_data_file(paradict,"rfa","jld2",dataloc)
everyconfig = []
particles = []
for file in allfiles
	parts = get_params_dict_from_filename(file)["particles"]
	append!(particles,[parts])
	append!(everyconfig,[read_data_jld2(file,dataloc)[1]["configs"]])
end
#
plotting_particles = []
allfits = []
axisbins = 300
for (i,config) in enumerate(everyconfig)
	parts = particles[i]
	println(parts)
	rm = sqrt(2*parts*3)
	#get_occupancy(config,rm; points=axisbins,title_string="N = $parts")
	raddens = radial_density_full(config,rm; rend="max",points=axisbins,labelstring="$parts",if_plot=false)
	#
	fitparams = fit_gaussian_raddens(parts,raddens,axisbins; if_plot=false)
	if 5 < parts < 17
		append!(allfits,[fitparams[1]])
		append!(allfits,[fitparams[2]])
		append!(plotting_particles,[parts,parts])
	elseif parts > 16
		append!(allfits,[fitparams[1]])
		append!(allfits,[fitparams[2]])
		append!(allfits,[fitparams[3]])
		append!(plotting_particles,[parts,parts,parts])
	else
		append!(allfits,[fitparams])
		append!(plotting_particles,[parts])
	end
	#
end
#=
fig = figure()
scatter(plotting_particles,[allfits[i][1] for i in 1:length(plotting_particles)])
xlabel("Particles")
ylabel("Width of Ring / rm")
=#
#
fig = figure()
scatter(plotting_particles,[allfits[i][2] for i in 1:length(plotting_particles)],label="Reverse")
xlabel("Particles")
ylabel("Location of Ring / rm")
=#

#=

alloccs = []
for (i,parts) in enumerate(particles)
	config = everyconfig[i]
	rm = sqrt(2*parts*3)
	raddens = radial_density_full(config,rm; rend="max",points=axisbins,labelstring="$parts",if_plot=false)
	if 5 < parts < 17
		left_occ = integrate(raddens[1][1:Int(axisbins/2)],raddens[2][1:Int(axisbins/2)])
		right_occ = integrate(raddens[1][Int(axisbins/2):end],raddens[2][Int(axisbins/2):end])
		if parts <= 9
			occs = [left_occ,0.0,right_occ]
		else
			occs = [0.0,left_occ,right_occ]
		end
	elseif parts > 16
		left_occ = integrate(raddens[1][1:Int(axisbins/3)],raddens[2][1:Int(axisbins/3)])
		center_occ = integrate(raddens[1][Int(axisbins/3):Int(2*axisbins/3)],raddens[2][Int(axisbins/3):Int(2*axisbins/3)])
		right_occ = integrate(raddens[1][Int(2*axisbins/3):end],raddens[2][Int(2*axisbins/3):end])
		occs = [left_occ,center_occ,right_occ]
	else
		occs = [0.0,0.0,1.0]
	end
	append!(alloccs,[occs])
end
fig = figure()
labstring = ["origin-disk","ring-2","ring-1"]
for i in 1:3
	plot(particles,[alloccs[j][i] for j in 1:length(particles)],"p",label="$(labstring[i])")
end
legend()
xlabel("Particles")
ylabel("Occupation")
=#
#=
distance_btw_rings = []
local_particles = []
for i in 1:Int(length(plotting_particles)/2)
	if i == 8
		continue
	else
		println(plotting_particles[2*(i-1)+2],", ",plotting_particles[2*(i-1)+1])
		distbtw = abs(allfits[2*(i-1)+2][2] - allfits[2*(i-1)+1][2])
		append!(distance_btw_rings,[distbtw])
		append!(local_particles,[plotting_particles[2*(i-1)+2]])
	end
end
scatter(local_particles,distance_btw_rings)
=#





# Location of edges for Laughlin
#=partsss = [i for i in 5:30]
edges = [0.686007068731582, 0.7526144437351984, 0.7325654184434776, 0.7416893053012075, 0.7509802474404931, 0.7472506221284726, 0.7733763808284947, 0.7831085483799342, 0.8078296756174796, 0.7825334529355564, 0.799055760880487, 0.8025235335176291, 0.8315009730608518, 0.8031209327915373, 0.8308061727475506, 0.8181201376432033, 0.85149171199076, 0.8327787217982823, 0.8378196375522812, 0.8482237544988176, 0.8323290375268354, 0.8477416349471587, 0.8489149479914027, 0.8437840871136896, 0.8522205951470362, 0.8431593803625228]
plot(partsss,edges,"-p",c="r",label="Laughlin")=#

#= for particles = [i for i in 4:20]
partsss1 = [i for i in 4:20]
rmax_edges = [1.55835484,1.495859155858919, 1.435630612727625, 1.37945930696872, 1.3654778777900691, 1.3864465233636276, 1.3130718120943006, 1.3447262795261632, 1.3038602437553868, 1.3199477842162546, 1.283222341188332, 1.2888736012325286, 1.2742575476804927, 1.2755683525810448, 1.2522972698797814, 1.2562847066293055, 1.2856352734705268]
plot(partsss1,rmax_edges,"-p",c="r",label="Laughlin (max radius)")
legend()
=#


























"fin"

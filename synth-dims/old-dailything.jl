if true
include("long-range-ttn.jl")
#include("fqh_effective.jl")
#include("time_evolution.jl")
include("../other-funcs/data-storage-funcs.jl")
#include("reproduce-other.jl")
#using PyPlot,Observers,ITensorTDVP,LsqFit
end
#=
alpha_start = 0.1
alpha_end = 0.15
alpha_count = 5
alpha_vals = [alpha_start + (i-1)*(alpha_end-alpha_start)/(alpha_count-1) for i in 1:alpha_count]
change = 0.0001
density = 5/40
Ls = [8,24,40,56,72,88,96]
=#

function percent_change(values)
    return [(values[i+1] - values[i]) for i in 1:length(values)-1]
end

function find_quadratic_region(values)
	pchange = percent_change(values)
    start_index = findfirst(x -> x > 1e-3, pchange) + 1
    end_index = findfirst(x -> pchange[x+1] < pchange[x],collect(1:length(pchange)-1))
	#println("Start Index = ",start_index,", End Index = ",end_index)
    return start_index, end_index
end

function fit_quadratic(times,values,start_index,end_index)
    x = times[start_index:end_index]
    y = values[start_index:end_index]
    model(x, p) = p[1].*((x .+ p[2]).^2) .+ values[start_index]
    fit = curve_fit(model, x, y, [0.5,0.5])
    return fit.param
end

function find_slope(times,values; kwargs...)
	if_plot = get(kwargs, :if_plot, false)
    start_index, end_index = find_quadratic_region(values)
    a,b = fit_quadratic(times,values,start_index,end_index)

	if if_plot
		fig = figure()
        plot(times, values, "-p")
        plot(times[start_index:end], a.*((times[start_index:end] .+ b).^2) .+ values[start_index])
        title("Quadratic Fit with Slope = $(round(a,digits=4))")
        xlabel("Times")
        ylabel("Value")
		ylim(-0.1*maximum(values),1.5*maximum(values))
    end

    return a
end

function find_sforderparam_transition(all_files,loc; kwargs...)
	if_plot = get(kwargs, :if_plot, false)

	fillings = []
	sf_orderparams = []
	file_params = get_params_dict_from_filename(all_files[1])
	for f in all_files
		data,metadata = read_data_jld2(f,loc)
		psi = data["mps"]
		rho = data["densmat"]
		chi = metadata["chi"]
		if chi == 0.0
			continue
		else
			append!(fillings,[metadata["nbosons"]/(file_params["L"]*file_params["nflavors"]*metadata["chi"])])
			sf_op = abs(momentum_occupation(psi,1,0.0; densmat=rho)[2][1])
			append!(sf_orderparams,[sf_op])
		end
	end

	deriv_sf_orderparams = percent_change(sf_orderparams)
	minval,min_loc = findmin(deriv_sf_orderparams)

	if if_plot
		#fig = figure()
		plot(fillings,sf_orderparams ./ (file_params["L"]*file_params["nflavors"]),"-p",label="L = $(file_params["L"]), nf=$(file_params["nflavors"])")
		#title("SF Order Parameter vs Filling: L = $(file_params["L"]), nflavors = $(file_params["nflavors"])")
		#plot([fillings[min_loc+1],fillings[min_loc+1]],[minimum(sf_orderparams),maximum(sf_orderparams)],label="Transition Point",c="r")
		#legend()
	end

	return fillings[min_loc+1],fillings,sf_orderparams
end

function find_current_transition(all_files,loc; kwargs...) # adding HC part of current so not exactly equal to current
	if_plot = get(kwargs, :if_plot, false)

	fillings = []
	currents = []
	file_params = get_params_dict_from_filename(all_files[1])
	for f in all_files
		data,metadata = read_data_jld2(f,loc)
		psi = data["mps"]
		rho = data["densmat"]
		chi = metadata["chi"]
		if chi == 0.0
			continue
		else
			append!(fillings,[metadata["nbosons"]/(file_params["L"]*file_params["nflavors"]*metadata["chi"])])
			loc_curr = real(get_current(data["mps"];alpha=metadata["chi"],densmat=data["densmat"],if_addhc=true)[1])
			append!(currents,[loc_curr])
		end
	end

	deriv_currents = percent_change(currents)
	minval,min_loc = findmin(deriv_currents)

	if if_plot

		fig = figure()
		plot(fillings,currents,"-p",label="Real")#,label="L = $(file_params["L"]), nf=$(file_params["nflavors"])")
		title("Current vs Filling: L = $(file_params["L"]), nflavors = $(file_params["nflavors"])")
		#plot([fillings[min_loc+1],fillings[min_loc+1]],[minimum(currents),maximum(currents)],label="Transition Point",c="r")
		legend()
	end

	return fillings[min_loc+1],fillings,currents
end

for if_periodic_phys in [true,false]
	for if_periodic_virt in [true,false]
		if if_periodic_phys && !if_periodic_virt
			continue
		end

		if if_periodic_phys && if_periodic_virt
			dataloc = get_folder_location("cluster-data/synth-dims/torus")
			ver = "Torus"
		elseif if_periodic_phys || if_periodic_virt
			dataloc = get_folder_location("cluster-data/synth-dims")
			ver = "Cylinder"
		elseif !if_periodic_phys && !if_periodic_virt
			dataloc = get_folder_location("cluster-data/synth-dims/obc")
			ver = "OBC"
		end

		for (i,stren) in enumerate([0.0,10.0])
			files = find_data_file(Dict([("layers",6),("particles",4),("onsite_strength",stren),("hopping_anisotropy",1.0)]),"ttn",dataloc; output_level=0)
			if ver == "Torus"
				data,metadata = read_data_jld2(files[1],dataloc)
				thisloc = i*5
				println(thisloc)
				occs = get_occupancy(data["ttn"]; densmat=data["densmat"], plot_title="$ver, LR=$stren, $thisloc",if_plot=true)
			else
				for f in files
					data,metadata = read_data_jld2(f,dataloc)
					alpha = metadata["alpha"]
					periodic_shift = alpha != 0.125
					ver_shift = ver == "Cylinder"
					thisloc = 2*ver_shift + periodic_shift + (i-1)*5 + 1
					println(thisloc)
					occs = get_occupancy(data["ttn"]; densmat=data["densmat"], plot_title="$ver, LR=$stren, alpha=$(round(alpha,digits=3)), $thisloc",if_plot=true)
				end
			end

			
			
		end

	end
end




#=
loc = get_folder_location("cluster-data/synth-dims","fzj")
for L in [6,7,8]#10,16
	for nflavors in [4,5,6,7,8]
		for nbosons in [3,4,5,8]
			params_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",nbosons)])
			all_files = find_data_file(params_dict,"mps",loc)
			if length(all_files) < 2
				continue
			end
			chis = []
			sfs = []
			for f in all_files
				data,metadata = read_data_jld2(f,loc)
				rho = data["densmat"]
				append!(chis,[metadata["chi"]])
				append!(sfs,[real(2*sum(rho))/(L*nflavors)])
			end
			plot(chis,sfs,"-p",label="L=$L,nf=$nflavors,nb=$nbosons")
			legend()
		end
	end
	xlabel("Flux Density")
	ylabel("SF Order Parameter")
end
=#

# 6,4,3 = 69
# 6,5,3 = 50
# 8,6,4 = 78


#=exc_loc = get_folder_location("cluster-data/synth-dims/higher-states","fzj")
loc = get_folder_location("cluster-data/synth-dims")

L = 6
nflavors = 5
nbosons = 3

for (L,nflavors,nbosons) in [(8,6,4)]#(6,5,3),(6,4,3),
	params_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",nbosons)])

	if true # this finds all the files and removes them if not converged
		gs_files = find_data_file(params_dict,"mps",loc)
		chis = []
		badfiles = []
		for f in gs_files
			namedata = get_params_dict_from_filename(f)
			filling = round(nbosons/(L*nflavors*namedata["chi"]),digits=4)
			data,metadata = read_data_jld2(f,loc)
			allvars = (metadata["observer"]).nrg_var[2:end]
			if allvars[end] > 1E-6
				#=plot(allvars,label="$filling")
				yscale("log")
				legend()=#
				append!(badfiles,[f])
			elseif namedata["chi"] == 0.0
				append!(badfiles,[f])
			elseif nbosons/(L*nflavors*namedata["chi"]) < 0.2
				append!(badfiles,[f])
			else
				append!(chis,[namedata["chi"]])
			end
		end
		filter!(x -> x ∉ badfiles, gs_files)
	end

	#
	deltas = []
	localchis = []
	fids = []
	#bonddims = []
	bulkdenss = []
	for i in 2:length(gs_files)
		println(i/length(gs_files))
		data1,metadata1 = read_data_jld2(gs_files[i-1],loc)
		data2,metadata2 = read_data_jld2(gs_files[i],loc)
		delta = abs(metadata1["chi"] - metadata2["chi"])/2
		append!(deltas,[delta])
		localchi = (metadata1["chi"] + metadata2["chi"])/2
		append!(localchis,[localchi])
		#fid = (-1/(L*nflavors)) * log(abs2(inner(data1["mps"],data2["mps"])))
		#append!(fids,[fid])
		#append!(bonddims,maxlinkdim(data1["mps"]))
		append!(bulkdenss,[bulk_density(data1["mps"],0,1)])
	end
	#=
	fig = figure()
	plot(nbosons ./ ((L*nflavors) .* localchis),1 .- exp.(-1 .* fids),"-p")
	plot(nbosons ./ ((L*nflavors) .* localchis),deltas,"-p")
	xlabel("Filling")
	ylabel("Infidelity (1-exp(-F))")
	title("L = $L, nflavors = $nflavors, nbosons = $nbosons")
	=#
	fig2 = figure()
	plot(nbosons ./ ((L*nflavors) .* localchis),bulkdenss,"-p")
	xlabel("Filling")
	ylabel("Bond Dim")
	title("L = $L, nflavors = $nflavors, nbosons = $nbosons")
end
=#


#=
loc = get_folder_location("cluster-data/synth-dims","fzj")
exc_loc = get_folder_location("cluster-data/synth-dims/higher-states","fzj")
L = 6
nflavors = 5
nbosons = 3
params_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",nbosons)])
gs_files = find_data_file(params_dict,"mps",loc)
exc_files = find_data_file(params_dict,"mps",exc_loc)

idx = 0
for f_exc in exc_files[21:40]
	data_exc,metadata_exc = read_data_jld2(f_exc,exc_loc)
	psi_exc = data_exc["mps"]
	filling_exc = nbosons/(L*nflavors*metadata_exc["chi"])
	fig = figure()
	max_overlap = 0.0
	for f_gs in gs_files
		data_gs,metadata_gs = read_data_jld2(f_gs,loc)
		psi_gs = data_gs["mps"]
		filling_gs = nbosons/(L*nflavors*metadata_gs["chi"])

		if isinf(filling_gs) || isinf(filling_exc)
			continue
		end

		overlap = abs2(inner(psi_exc,psi_gs))
		overlap > max_overlap ? max_overlap = overlap : nothing
		 
		global idx += 1
		if (filling_gs < 0.3 && filling_exc < 0.3) || (filling_gs > 0.3 && filling_exc > 0.3)
			pc = "b"
		elseif filling_gs == 0.3 || filling_exc == 0.3
			pc = "r"
		else
			pc = "g"
		end
		println(idx)
		scatter([filling_gs],[overlap],c=pc)
	end
	plot([filling_exc,filling_exc],[0,max_overlap],c="k")
end
title("Overlap of excited state at black line with other GSs")
=#

#=
L = 6
nflavors = 5
nbosons = 3

col = ["b","g","r","c","m","y","k","w"]

if false
params_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",nbosons)])
loc = get_folder_location("cluster-data/synth-dims","fzj")
all_files = find_data_file(params_dict,"mps",loc)
display(all_files)
sforderparams = zeros(length(all_files))
fillings = zeros(length(all_files))
corrlengths = [zeros(nflavors) for i in 1:length(all_files)]
#
for (idx,f) in enumerate(all_files)
	name_data = get_params_dict_from_filename(f)
	filling = name_data["nbosons"] / (name_data["L"] * name_data["nflavors"] * name_data["chi"])
	#if filling < 0.4
	#	continue
	#end
	fillings[idx] = filling
	data,metadata = read_data_jld2(f,loc)
	wavefunc = data["mps"]
	if "densmat" in keys(data)
		rho = data["densmat"]
	else
		rho = density_matrix(wavefunc)
	end
	#sforderparams[idx] = abs(momentum_occupation(wavefunc,1,0.0; densmat=rho)[2][1])
	corrlengths[idx][:] = physical_distance_correlation(wavefunc; densmat=rho)[3]
	if idx > 1
		for i in 1:edge_length
			plot(fillings[idx-1:idx],[corrlengths[idx-1][i],corrlengths[idx][i]],"-p",c=col[i])
		end
	end
end
end
fig = figure()
xvals = fillings#num_parts ./ (fillings .* (edge_length^2))
for i in 1:edge_length
	plot(xvals,[corrlengths[j][i] for j in 1:length(all_files)],"-p",label="$i")
end
xlabel("Filling")
ylabel("Physical Correlation Length")
legend()
=#

function integral_part(correlation_length,cutoff_radius,counts=500; kwargs...)
	xs = range(0.1,stop=cutoff_radius,length=counts)
	ys = exp.(-xs ./ correlation_length) ./ (xs .^ 0.5)
	return integrate(xs,ys)
end

if false
	ll = 6
	np = 4
	magfield = 2*np / (2^ll)
	pbc = true
	anises = [1.0,0.8,0.7,0.6,0.5,0.4,0.35]
	trans_strens = [zeros(length(anises)),zeros(length(anises))]
	for (jj,anis) in enumerate(anises)
		pdict = Dict([("layers",ll),("particles",np),("if_periodic_phys",pbc),("lr",7),("alpha",magfield),("hopping_anisotropy",anis)])
		whichfiles = find_data_file(pdict,"ttn",get_folder_location("cluster-data/synth-dims"); output_level=false)
		#
		cdwsfs = [zeros(length(whichfiles)),zeros(length(whichfiles))]
		lrs = zeros(length(whichfiles))
		periods = [1,2]
		theta = 0.0
		for (idx,f) in enumerate(whichfiles)
			data,metadata = read_data_jld2(f,get_folder_location("cluster-data/synth-dims"); output_level=false)
			psi = data["ttn"]
			dens = "densmat" in keys(data) ? data["densmat"] : nothing
			cdwsfs[1][idx] = abs(cdw_structure_factor(dens,(periods[1]*cos(theta),periods[1]*sin(theta)),psi))
			cdwsfs[2][idx] = abs(cdw_structure_factor(dens,(periods[2]*cos(theta),periods[2]*sin(theta)),psi))
			lrs[idx] = 2*metadata["onsite_strength"]
		end
		
		# transition is the first time the derivative is greater than 20% of the maximum derivative
		cdw_derivs = abs.((cdwsfs[2][2:end] .- cdwsfs[2][1:end-1]) ./ cdwsfs[2][1])
		trans_loc = findfirst(x -> cdw_derivs[x] > 0.2*(maximum(cdw_derivs)),1:length(cdw_derivs))
		transition_strength = lrs[trans_loc]
		trans_strens[2][jj] = transition_strength

		# period 1 transition is when the value of cdwsf is within 1% of the maximum and after the period 2 transition
		period1_onset = findfirst(x -> x > trans_loc && isapprox(cdwsfs[1][x],cdwsfs[1][1],atol=0.5*(cdwsfs[1][1] - minimum(cdwsfs[1]))),1:length(cdwsfs[1]))
		trans_strens[1][jj] = lrs[period1_onset]

		#
		fig2,ax1 = subplots()
		ax1.plot(lrs,cdwsfs[1],c="b")
		ax1.set_ylabel("Period $(periods[1])",c="b")
		ax1.scatter(lrs[period1_onset],cdwsfs[1][period1_onset],c="g")
		xlabel("Onsite Strength")
		title("CDW Structure Factor")

		ax2 = ax1.twinx()
		ax2.plot(lrs,cdwsfs[2],c="r")
		ax2.set_ylabel("Period $(periods[2])",c="r")
		ax2.scatter(lrs[trans_loc],cdwsfs[2][trans_loc],c="g")
		#
	end
	fig = figure()
	#plot(trans_strens[1],anises,"-p",label="Period 1 Onset")
	plot(trans_strens[2],anises,"-p",label="Period 2 Onset")
	legend()
	xlabel("Transition Strength")
	ylabel("Hopping Anisotropy")
	yscale("log")
end


if false

ll = 4
np = 2
pbc = true
pdict = Dict([("layers",ll),("particles",np),("if_periodic_phys",pbc),("lr",3)])
whichfiles = find_data_file(pdict,"ttn",get_folder_location("cluster-data/synth-dims"); output_level=false)

#fig = figure("pyplot_subplot_column",figsize=(5,20))

strens = [range(0.01,2.0,length=10); range(0.01,2.0,length=10)[2:end] .+ 2]

datapoints = 20
anises = [0.8,0.7,0.6,0.5]#[1.0,0.8,0.7,0.6,0.5,0.4,0.35,0.3]
howmany = length(anises)
plotting_dict = Dict()
for anis in anises
	plotting_dict[anis] = 0
end

for (idx,f) in enumerate(whichfiles)
	data,metadata = read_data_jld2(f,get_folder_location("cluster-data/synth-dims"); output_level=false)
	psi = data["ttn"]
	dens = "densmat" in keys(data) ? data["densmat"] : nothing
	anis = "hopping_anisotropy" in keys(metadata) ? metadata["hopping_anisotropy"] : 1.0
	uu = metadata["onsite_strength"]
	if uu > 4.0 || anis > 1.0 || anis < 0.1
		continue
	end
	
	#
	
	#=append!(anises,[anis])
	append!(nrgs,[-metadata["energies"][end]])
	scatter([anis],[-metadata["energies"][end]],c="b")
	xlabel("Anisotropy")
	ylabel("Energy")
	xscale("log")=#

	#=cols = ["b","r","g"]

	for ang in [0.0,0.25,0.5]
		if ang == 0.0
			col = "b"
		elseif ang == 0.25
			col = "r"
		else
			col = "g"
		end
		cdwsf = cdw_structure_factor(dens,(cos(ang*pi),sin(ang*pi)),psi)
		if idx == 1
			scatter([anis],[abs(cdwsf)],c=col,label="$(round(ang,digits=2))π")
		else
			scatter([anis],[abs(cdwsf)],c=col)
		end
	end
	xscale("log")
	legend()=#
	#get_occupancy(psi; if_plot=true,plot_title="Hopping Anis = $anis, LR Strength = $uu",densmat=dens)

	if !(anis in anises)
		continue
	end

	if anis == 1.0 && !(uu in strens)
		continue
	end

	uu = round(2*uu,digits=3)

	occs = get_occupancy(psi; plot_title="Hopping Anis = $anis, LR Strength = $uu",densmat=dens,if_plot=false)

	
	plotting_dict[anis] += 1
	shift = findfirst(x -> anises[x] == anis,1:howmany)
	thisloc = Int(howmany*(plotting_dict[anis]-1) + shift)
	#println("Plotting Anis = $anis, LR Stren = $uu t, this loc = ",thisloc)
	subplot(datapoints,howmany,thisloc)
	imshow(occs)
	#=fig = figure()
	xs = collect(Iterators.flatten([i .* ones(8) for i in 1:8]))
	ys = collect(Iterators.flatten([[1,2,3,4,5,6,7,8] for i in 1:8]))
	zs = collect(Iterators.flatten([occs[i,:] for i in 1:8]))
	plot3D(xs,ys,zs,"-p")=#
	#
	if plotting_dict[anis] == 1
		title("H Anis=$anis")
		if anis == anises[end]
			scatter([occs[1,1]],[occs[1,1]],label="LR=$uu t",c="b",s=1.0)
			legend(loc=2,bbox_to_anchor=(1.05,1.0))
		end
	elseif anis == anises[end]
		scatter([occs[1,1]],[occs[1,1]],label="LR=$uu t",c="b",s=1.0)
		legend(loc=2,bbox_to_anchor=(1.05,1.0))
	end
	#
	

	#=if anis > 1.0
		lat = TTNKit.physical_lattice(TTNKit.network(psi))
		new_dens = zeros(2^ll,2^ll) .* im
		for j in 1:8
			for s in 1:8
				site1 = TTNKit.linear_ind(lat,(j,s))
				redone_site1 = TTNKit.linear_ind(lat,(s,j))
				for jj in 1:8
					for ss in 1:8
						site2 = TTNKit.linear_ind(lat,(jj,ss))
						redone_site2 = TTNKit.linear_ind(lat,(ss,jj))
						new_dens[redone_site1,redone_site2] = dens[site1,site2]
					end
				end
			end
		end
		append!(alloccs,[new_dens])
	else
		append!(alloccs,[dens])
	end
	
	fig = figure()
	imshow(abs.(alloccs[end]))
	title("Anisotropy = $anis")
	colorbar()=#

	#=zeromomocc = TTNKit.maxlinkdim(psi)#abs(sum(dens)) / ((2^ll) * np)
	scatter([anis],[zeromomocc],c="b")
	xlabel("Anisotropy")
	ylabel("Zero Momentum Occupancy")
	xscale("log")#
	occs = get_occupancy(psi; if_plot=true,plot_title="Hopping Anis = $anis",densmat=dens)
	fig = figure()
	plot(1:8,occs[4,:],label="Synth")
	plot(1:8,occs[:,4],label="Phys")
	legend()
	title("Filling = 1/2, Anisotropy = $anis")
	ylabel("Occupancy")=#


#=
lat = TTNKit.physical_lattice(TTNKit.network(psi))
phys_length = ll % 2 == 0 ? Int(sqrt(2^ll)) : Int(sqrt(2^(ll+1))/2)
virt_length = ll % 2 != 0 ? Int(sqrt(2^(ll-1))) : Int(sqrt(2^ll))


phys_corrs = zeros(phys_length,virt_length)
for s in 1:virt_length
	for j in 1:phys_length
		site1 = TTNKit.linear_ind(lat,(1,s))
		site2 = TTNKit.linear_ind(lat,(j,s))
		corr_val = dens[site1,site2]
		corr_normalization = sqrt(dens[site1,site1] * dens[site2,site2])
		corr_val /= corr_normalization
		phys_corrs[j,s] = abs(corr_val)
	end
end
fig = figure()
for i in 1:Int(virt_length/2)#2:Int(virt_length/2)
	plot(0:phys_length-1,phys_corrs[:,i],"-p",label="$i")
end
xlabel("Physical Distance")
ylabel("Correlation")
legend()
#title("Greens Func LR = 7, Filling = 1/2")
title("Greens Func LR = 7, Onsite = 10t, Anisotropy = $anis")
yscale("log")
#
virt_corrs = zeros(phys_length,virt_length)
for j in 1:phys_length
	for s in 1:virt_length
		site1 = TTNKit.linear_ind(lat,(j,1))
		site2 = TTNKit.linear_ind(lat,(j,s))
		corr_val = dens[site1,site2]
		corr_normalization = sqrt(dens[site1,site1] * dens[site2,site2])
		corr_val /= corr_normalization
		virt_corrs[j,s] = abs(corr_val)
	end
end
fig2 = figure()
for i in 1:Int(phys_length/2)
	plot(0:virt_length-1,virt_corrs[i,:],"-p",label="$i")
end
xlabel("Virtual Distance")
ylabel("Correlation")
legend()
title("Greens Func along Virtual Dimension for Filling = 1/2")
yscale("log")
#

phys_currents = zeros(phys_length,virt_length)
for j in 1:phys_length
	next_phys = j + 1
	if j == phys_length
		if pbc
			next_phys = 1
		else
			continue
		end
	end
	for s in 1:virt_length
		site2 = TTNKit.linear_ind(lat,(j,s))
		site1 = TTNKit.linear_ind(lat,(next_phys,s))
		current_val = imag(dens[site1,site2] - dens[site2,site1])
		current_normalization = dens[site1,site1] + dens[site2,site2]
		current_val /= current_normalization
		phys_currents[j,s] = real(current_val)
	end
end
fig3 = figure()
for i in 1:phys_length
	plot(1:virt_length,phys_currents[i,:],"-p",label="$i")
end
xlabel("Virtual Dimension")
ylabel("Physical Current")
legend()
#title("Current LR = 7, Filling = 1/2")
title("Current LR = 7, Onsite = 10t, Anisotropy = $anis")
#
virt_currents = zeros(phys_length,virt_length)
for s in 1:virt_length-1
	next_virt = s + 1
	for j in 1:phys_length
		site2 = TTNKit.linear_ind(lat,(j,s))
		site1 = TTNKit.linear_ind(lat,(j,next_virt))
		current_val = imag(dens[site1,site2] - dens[site2,site1])
		current_normalization = dens[site1,site1] + dens[site2,site2]
		current_val /= current_normalization
		virt_currents[j,s] = real(current_val)
	end
end
fig4 = figure()
for i in 1:Int(virt_length/2)
	plot(1:phys_length,virt_currents[:,i],"-p",label="$i")
end
xlabel("Physical Dimension")
ylabel("Virtual Current")
legend()
title("Current along Virtual Dimension for Filling = 1/2")
=#

end
end

if false
	fig = figure()
	col = ["b","g","r","c","m","y","k","w"]
	all_nrgs = [[],[]]
	all_twists = [[],[]]
	for layers in [4,5,6,7,8]
		if layers % 2 == 0
			phys_edge_length = Int(sqrt(2^layers))
			virt_edge_length = phys_edge_length
		else
			phys_edge_length = Int(sqrt(2^(layers+1))/2)
			virt_edge_length = Int(sqrt(2^(layers-1))/2)
		end
		num_parts = layers % 2 == 0 ? Int(sqrt(2^(layers))/2) : Int(sqrt(2^(layers+1))/2)

		if true
			params_dict = Dict([("layers",layers),("particles",num_parts),("if_periodic_phys",true)])
			loc = get_folder_location("cluster-data/synth-dims")
			all_files = find_data_file(params_dict,"ttn",loc; output_level=false)
			#display(all_files)
			#nrgs = []
			#twists = []
			#fillings = zeros(length(all_files))
			intparts = []
			fillings = []
			zeromoms = []
			halffillingvalue = nothing
			sfvalue = nothing
			for (idx,f) in enumerate(all_files)
				name_data = get_params_dict_from_filename(f)
				filling = name_data["particles"] / (2^(name_data["layers"]) * name_data["alpha"])
				onsite_strength = "onsite_strength" in keys(name_data) ? name_data["onsite_strength"] : 0.0
				hopping_anisotropy = "hopping_anisotropy" in keys(name_data) ? name_data["hopping_anisotropy"] : 1.0

				if onsite_strength != 0.0
					continue
				end
				if hopping_anisotropy != 1.0
					continue
				end

				#=
				if isapprox(filling,0.5,atol=0.01)
					col = "b"
				elseif name_data["alpha"] == 0.0
					col = "r"
				elseif isapprox(filling,0.25,atol=0.001)
					col = "g"
				elseif isapprox(filling,0.3,atol=0.01)
					col = "c"
				else
					continue
				end
				=#
				data,metadata = read_data_jld2(f,loc;output_level=false)
				wavefunc = data["ttn"]
				rho = data["densmat"]
				
				#fig = figure()
				#dists,corrs,corrlens = physical_distance_correlation(wavefunc; densmat=rho,if_plot=true)
				#title("Filling = $(round(filling,digits=4))")
				#get_occupancy(wavefunc; plot_title="Filling = $(round(filling,digits=4)), Bond Dim=$(maxlinkdim(wavefunc))",densmat=rho)
				#=if any(corrlens .> 10.0)
					continue
				end=#

				append!(fillings,[filling])
				#allints_parts = [integral_part(corrlens[i],phys_edge_length/2) for i in 1:length(corrlens)]
				zeromomoccs = abs(sum(rho)) / (2^(layers) * num_parts)
				#append!(intparts,[2*sum(allints_parts)])
				append!(zeromoms,[zeromomoccs])
				if isapprox(filling,0.5,atol=0.01)
					halffillingvalue = zeromomoccs
				elseif name_data["alpha"] == 0.0
					sfvalue = zeromomoccs
				end
				#=println("Integral = ",sum(allints_parts)," Zero Mom Occ = ",zeromomoccs)
				if length(zeromoms) > 1
					#plot(fillings[end-1:end],intparts[end-1:end],"-p",c="b")
					plot(fillings[end-1:end],zeromoms[end-1:end],"-p",c="b")
					xlabel("Filling")
					#ylabel("ODLRO")
				end=#
				

				#append!(all_nrgs[layers-3],[metadata["energies"][end]])
				#append!(all_twists[layers-3],[twist])

			end
			scatter(2^layers,sfvalue,c="r")
			scatter(2^layers,halffillingvalue,c="b")
			println("Half filling value = ",halffillingvalue)
			length(zeromoms) > 0 ? scatter(2^layers,minimum(zeromoms),c="c") : nothing
			#all_nrgs[layers-3] .-= all_nrgs[layers-3][1]
		end
		phys_edge_length = layers % 2 == 0.0 ? Int(sqrt(2^layers)) : Int(sqrt(2^(layers+1))/2)
		virt_edge_length = layers % 2 == 0.0 ? Int(sqrt(2^layers)) : Int(sqrt(2^(layers-1))/2)
		#all_nrgs[layers-3] .*= phys_edge_length / virt_edge_length
	end
	#=
	all_hms = 2 .* [all_nrgs[i] ./ (all_twists[i] .^2) for i in 1:2]
	helicity_moduli = [mean(all_hms[i][2:end]) for i in 1:2]
	errors = [std(all_hms[i][2:end]) for i in 1:2]
	plot(all_twists[1],all_nrgs[1],"-p",label="Ns = $(2^4)")
	plot(all_twists[2],all_nrgs[2],"-p",label="Ns = $(2^5)")
	legend()
	xlabel("Twist Angle")
	title("Helicity Modulus = $(round(helicity_moduli[1],digits=4)) ± $(round(errors[1],digits=4))")=#
end


if false
halfvals = [0.17246013515736627,0.10920578094440274,0.0633053552666524,0.02970695637101318]
layers = [4,5,6,7]
expfit(x,p) = p[1] .* exp.(-p[2] .* x) .+ p[3]
fit = curve_fit(expfit,2 .^ layers,halfvals,[0.5,0.5,0.1])
plot(2 .^ layers,halfvals,"p")
xs = 2 .^ range(4,stop=7,length=100)
#plot(xs,expfit(xs,fit.param))
title("Decay Coeff = $(round(fit.param[2],digits=4)), Thermo Limit = $(round(fit.param[3],digits=4))")
yscale("log")
xscale("log")
end

#=
allberries = [zeros(edge_length-2,edge_length-2) for i in 1:length(all_files)]
avgberries = [0.0 for i in 1:length(all_files)]
fillings = [0.0 for i in 1:length(all_files)]
for (idx,f) in enumerate(all_files)
	data,metadata = read_data_jld2(f,loc)
	wavefunc = data["ttn"]
	fillings[idx] = metadata["particles"] / (2^(metadata["layers"]) * metadata["alpha"])
	for i in 3:edge_length-2
		for j in 3:edge_length-2
			println(i,", ",j)
			allberries[idx][i-2,j] = closed_loop(wavefunc,(i,j))[1]
			if allberries[idx][i-2,j] < 0.0
				newval = 2*pi + allberries[idx][i-2,j]
				allberries[idx][i-2,j] = newval
			end
		end
	end
	avgberries[idx] = mean(allberries[idx])
end
#
plot(fillings,avgberries,"-p")
=#

#=
chis = []
all_states = []
for f in all_files
	metadata = read_data_jld2(f,loc)[2]
	append!(chis,[metadata["chi"]])
	append!(all_states,[read_data_jld2(f,loc)[1]["mps"]])
end
end
=#

#=
entrops = []
for (i,psi) in enumerate(all_states)
	u,s,v = svd(psi[Int(L/2)],linkind(psi,Int(L/2)))
	entrop = sum([s[n,n]^2 * log(s[n,n]^2) for n in 1:size(s)[1]])
	append!(entrops,[entrop])
end
plot(5 ./ (chis .* (L*nflavors)) ,log.(entrops),"-p")
=#

#=
currents = zeros(nflavors,length(all_files)) .* im
for (i,psi) in enumerate(all_states)
	ham_params = (if_periodic_phys=true,if_periodic_synth=false,centralflux_strength=centralflux_strength,tilt_strength=0.0)
	currents[:,i] = [calc_deriv(1,psi,s,Int(L/2),nflavors,chi,ham_params) for s in 1:nflavors]
end
#
hall_currents = abs.((currents[1,:] .- currents[nflavors,:])/2)
plot(chis,hall_currents,"-p")
xlabel("Flux Density")
ylabel("Hall Current")
=#



#=
toberemoved = []
for (idx,f) in enumerate(all_files)
	alpha = get_params_dict_from_filename(f)["alpha"]
	if !isapprox(10/(alpha*10*16),1.0,atol=10^-3)
		append!(toberemoved,[idx])
	end
end
deleteat!(all_files,toberemoved)
#display(all_files)

alpha1,alpha2 = get_params_dict_from_filename(all_files[1])["alpha"],get_params_dict_from_filename(all_files[2])["alpha"]
data1,metadata1 = read_data_jld2(all_files[1],loc)
#data2,metadata2 = read_data_jld2(all_files[2],loc)
end

#obs_measures=Dict("denspols" => current_density_polarization)

mdim = 100
if_save_data = false

if_current = false
#time_end = 30.0
time_change = 0.00001
time_count = 10
time_end = time_count * time_change

strens = [0.05,0.025,0.01,0.0075,0.005]
current_strength = 0.0
#for (i,time_change) in enumerate(strens)
#time_end = time_count * time_change
tevo_params = (if_GScheck=true,current_strength=current_strength,if_current=if_current,mdim=mdim,location=loc,if_save_data=if_save_data)
	rez,ham = execute_tevo(all_files[1],time_end,time_change; tevo_params...)
	#slopes[i] = (rez["denspols"].results[end] - rez["denspols"].results[1])/(time_end)
	times = rez["times"].results
	#
	nrg_invs = mean([energy_variance(rez["states"].results[i],ham) for i in 1:length(times)])
	occ_invs = mean([first(occupancy_variance(rez["states"].results[i]; if_plot=false)) for i in 1:length(times)])
	println("Energy Variance = ",nrg_invs,", Occupation Variance = ",occ_invs)
	#
	allspacialpols = spacial_density_polarization(rez["occs"].results)
	spacial_limit = 1.0*get_params_dict_from_filename(all_files[1])["nbosons"]/get_params_dict_from_filename(all_files[1])["L"]
	fig2 = figure()
	plot(times,[spacial_limit for i in 1:length(times)],c="r",label="Spacial Limit")
	plot(times,[-spacial_limit for i in 1:length(times)],c="r")
	plot(times,allspacialpols,"-p")
	xlabel("Time")
	ylabel("Spacial Density Polarization")
	title("Current Strength = $(current_strength)")
	legend()
	#
#end

	#
	alldenspols = density_polarization(rez["occs"].results)
	virtual_limit = 1.0
	fig3 = figure()
	plot(times,[virtual_limit for i in 1:length(times)],c="r",label="Spacial Limit")
	plot(times,[-virtual_limit for i in 1:length(times)],c="r")
	plot(times,alldenspols,"-p",label="$mdim")
	legend()
	xlabel("Time")
	ylabel("Density Polarization")
	=#
	
	#=
	fig4 = figure()
	currents = get_current(rez["states"].results; alpha=alpha1,if_exp_part=true)
	plot(times,currents,"-p",label="$time_change")
	xlabel("Time")
	ylabel("Current")
	legend()
	#
	fig3 = figure()
	plot(times,(1/alpha1) .* (alldenspols ./ currents),"-p",label="$time_change")
	xlabel("Time")
	ylabel("Hall Imbalance")
	=#
#end

#

#=
alphas = []
dbds = []
for i in 1:Int(length(all_files)/2)
	data1,metadata1 = read_data_jld2(all_files[2*i-1],loc)
#	data2,metadata2 = read_data_jld2(all_files[2*i],loc)
	alpha1,alpha2 = metadata1["alpha"],metadata2["alpha"]
	#get_greenfunc(data1["ttn"],"phys"; plot_title="$(round(alpha1,digits=4))")
	#get_greenfunc(data1["ttn"],"virt"; plot_title="$(round(alpha1,digits=4))")
	#get_occupancy(data1["ttn"]; plot_title="$(round(alpha1,digits=4))")
#	append!(alphas,[0.5*(alpha1+alpha2)])
	append!(dbds,[deriv_bulk_dens(data1["ttn"],data2["ttn"],abs(alpha1-alpha2),1,2)])
end
plot(alphas,dbds,"-p")
xlabel("Flux Density")
ylabel("Deriv Bulk Density")
=#
#=
notcon = []
location_dict = Dict()
for (i,f) in enumerate(all_files)
	namedata = get_params_dict_from_filename(f)
	nf,L,nbosons,alpha = namedata["nflavors"],namedata["L"],namedata["nbosons"],namedata["alpha"]
	#filling = nbosons/(L*nf*alpha)
	if "$L,$nf" in collect(keys(location_dict))
		append!(location_dict["$L,$nf"],[i])
	else
		location_dict["$L,$nf"] = [i]
	end
end

for (k,v) in location_dict
	if length(location_dict[k]) < 5
		delete!(location_dict,k)
	end
end
#
toberm = []
for (idx,i) in enumerate(location_dict["16,10"])
	file = all_files[i]
	if occursin("mk",file)
		append!(toberm,[idx])
	end
end
deleteat!(location_dict["16,10"],toberm)
#

allfillings = []
alllocs = []
for (idx,i) in enumerate(location_dict["16,10"])
	file = all_files[i]
	namedata = get_params_dict_from_filename(file)
	alpha = namedata["alpha"]
	#filling = nbosons/(L*nf*alpha)
	append!(allfillings,[alpha])
	append!(alllocs,[idx])
end

noneighbors = []
for (i,v) in enumerate(allfillings)
	alldifs = deleteat!(abs.(v .- allfillings),i)
	if !any(alldifs .< 0.0002)
		append!(noneighbors,[alllocs[i]])
	end
end
deleteat!(location_dict["16,10"],noneighbors)


#
#
alphas = []
dbds = []
for i in 1:Int(length(location_dict["16,10"])/2)
	file1 = all_files[location_dict["16,10"][2*i-1]]
	file2 = all_files[location_dict["16,10"][2*i]]
	ttndict1,metadata1 = read_data_jld2(file1,loc)
	ttndict2,metadata2 = read_data_jld2(file2,loc)
	alpha1, alpha2 = metadata1["phi"]/(2*pi) ,metadata2["phi"]/(2*pi)
	alphadiff = alpha1 - alpha2
	append!(alphas,[0.5*(alpha1+alpha2)])
	dbd = deriv_bulk_dens(ttndict1["mps"],ttndict2["mps"],alphadiff)
	append!(dbds,[dbd])
end
#
ys = collect(minimum(dbds):0.01:maximum(dbds))
plot([1/16 for i in 1:length(ys)],ys)
scatter(alphas,dbds)
=#
#=
allpsis = Dict()
all_dbds = Dict()
all_nrgs = Dict()

for (idx,Lx) in enumerate(Ls)
	
	#wavefuncs = [[] for i in 1:alpha_count]
	#bulkdenses = [0.0 for i in 1:3*alpha_count]
	#alpha_vals = [0.0 for i in 1:length(all_files)]
	#derivbds = [0.0 for i in 1:alpha_count]
	#=
	if false
	for (idx,file) in enumerate(all_files)
		alpha = get_params_dict_from_filename(file)["alpha"]
		alpha_vals[idx] = alpha
		wavefunc = read_data_jld2(file,loc)[1]["mps"]
		bd = bulk_density(wavefunc,10)
		bulkdenses[idx] = bd
	end
	plot(alpha_vals,bulkdenses,"-p")
	end

	if false
	deriv_alphas = []
	for i in 1:Int(length(alpha_vals)/3)
		append!(derivbds,[(bulkdenses[3*(i-1)+3]-bulkdenses[3*(i-1)+1])/(alpha_vals[3*(i-1)+3]-alpha_vals[3*(i-1)+1])])
		append!(deriv_alphas,[alpha_vals[3*(i-1)+2]])
	end
	end
	=#
	nbosons = Int(Lx*5*density)
	center_alpha = Lx != 96 ? round(nbosons/((Lx-1)*5),digits=4) : 0.125
	seed_data = read_data_jld2("mps-if_periodic_phys-false-nflavors-5-alpha-$(center_alpha)-if_periodic_synth-true-nbosons-$(nbosons)-L-$Lx.jld2",loc)
	seed_nrg,minus_nrg,plus_nrg = 0.0,0.0,0.0
	try
		seed_nrg = seed_data[2]["final_energy"]/nbosons
		println("Found Seed NRG = ",seed_nrg)
	catch
		println("No final energy for Lx=$Lx Seed wavefunction")
		new_execution_dict = seed_data[2]
		seed_nrg = execute_mps(new_execution_dict["U1"],new_execution_dict["U2"],new_execution_dict["phi"],new_execution_dict["L"], new_execution_dict["nflavors"],new_execution_dict["nbosons"]; dict_to_symbols(new_execution_dict)...,if_nrg=true,psi_guess=seed_data[1]["mps"])/nbosons
		#seed_nrg = 0.0
	end
	plus_data = read_data_jld2("mps-if_periodic_phys-false-nflavors-5-alpha-$(center_alpha)-if_periodic_synth-true-nbosons-$(nbosons+1)-L-$Lx.jld2",loc)
	try
		plus_nrg = plus_data[2]["final_energy"]/(nbosons+1)
		println("Found Plus NRG = ",plus_nrg)
	catch
		println("No final energy for Lx=$Lx Plus wavefunction")
		#plus_nrg = 0.0
	end
	minus_data = read_data_jld2("mps-if_periodic_phys-false-nflavors-5-alpha-$(center_alpha)-if_periodic_synth-true-nbosons-$(nbosons-1)-L-$Lx.jld2",loc)
	try
		minus_nrg = minus_data[2]["final_energy"]/(nbosons-1)
		println("Found Minus NRG = ",minus_nrg)
	catch
		println("No final energy for Lx=$Lx Minus wavefunction")
		#minus_nrg = 0.0
	end
	
	all_nrgs["$Lx"] = [minus_nrg,seed_nrg,plus_nrg]
	
	#psi_p = run_mps_new_variable(seed_data[1]["mps"],seed_data[2],Dict([("nbosons",seed_data[2]["nbosons"]+1),("location",loc),("mdim",100)]))
	#psi_m = run_mps_new_variable(seed_data[1]["mps"],seed_data[2],Dict([("nbosons",seed_data[2]["nbosons"]-1),("location",loc),("mdim",100)]))

	#=
	for i in 1:alpha_count
		println(i/alpha_count)
		alph = alpha_vals[i]
		println(alph)
		#psi_mid = run_mps_new_variable(seed_data[1]["mps"],seed_data[2],Dict([("phi",2*pi*alph),("location",loc),("mdim",100)]))
		psi_m,psi_p = varied_alpha_wavefuncs(seed_data[1]["mps"],merge(seed_data[2],Dict([("phi",2*pi*alph),("location",loc)])),change)
		wavefuncs[i] = [psi_m,seed_data[1]["mps"],psi_p]
		bd_m = bulk_density(psi_m,Int(ceil(0.1*Lx)))
		bulkdenses[3*(i-1)+1] = bd_m
		bd_mid = bulk_density(seed_data[1]["mps"],Int(ceil(0.1*Lx)))
		bulkdenses[3*(i-1)+2] = bd_mid
		bd_p = bulk_density(psi_p,Int(ceil(0.1*Lx)))
		bulkdenses[3*i] = bd_p
		derivbds[i] = (bd_p - bd_m)/(2*change)
		println("DBD = ",derivbds[i])
	end
	=#
	
	#allpsis["$Lx"] = [psi_m,seed_data[1]["mps"],psi_p]
	#all_dbds["$Lx"] = derivbds
	#fig = figure()
	#plot(alpha_vals,derivbds,"-p",label="$Lx")
	#legend()
end
#

#if true
#fig = figure()
#scatter(new_new_new_alpha_vals,new_new_new_derivbds)
#title("Deriv Bulk Density vs Flux Density")
#end

#fig2 = figure()
#plot(Iterators.flatten([[alpha_vals[i]-change,alpha_vals[i],alpha_vals[i]+change] for i in 1:alpha_count]),bulkdenses,"-p")
#title("Bulk Density vs Flux Density")
=#































"fin"

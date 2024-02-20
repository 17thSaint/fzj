if true
include("long-range-ttn.jl")
#include("fqh_effective.jl")
include("time_evolution.jl")
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

loc = get_folder_location("cluster-data/synth-dims","fzj")
#=exc_loc = get_folder_location("cluster-data/synth-dims/higher-states","fzj")

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

if true
layers = 6
edge_length = Int(sqrt(2^layers))
num_parts = 4

col = ["b","g","r","c","m","y","k","w"]

if true
params_dict = Dict([("layers",layers),("particles",num_parts)])
loc = get_folder_location("cluster-data/synth-dims","fzj")
all_files = find_data_file(params_dict,"ttn",loc)
display(all_files)
sforderparams = zeros(length(all_files))
shortrange = zeros(length(all_files))
longrange = zeros(length(all_files))
fillings = zeros(length(all_files))
corrlengths = [zeros(Int(edge_length)) for i in 1:length(all_files)]
#
for (idx,f) in enumerate(all_files)
	name_data = get_params_dict_from_filename(f)
	filling = name_data["particles"] / (2^(name_data["layers"]) * name_data["alpha"])
	#if filling < 0.4
	#	continue
	#end
	fillings[idx] = filling
	data,metadata = read_data_jld2(f,loc)
	wavefunc = data["ttn"]
	if "densmat" in keys(data)
		println("Inside data")
		rho = data["densmat"]
	else
		println("Need to make")
		continue
		#rho = density_matrix(wavefunc)
	end
	momocc = rho
	middle = Int((layers^2)/2)
	longrange[idx] = abs(2*sum([sum(diag(momocc,i) + diag(momocc,-i)) for i in middle+1:layers^2-1]))
	shortrange[idx] = abs(2*sum([sum(diag(momocc,i) + diag(momocc,-i)) for i in 0:middle]))
	sforderparams[idx] = abs(2*sum(momocc))

	#
	#top = [abs(sum(diag(momocc,i))) for i in 1:layers^2-1]
	#fig = figure()
	#plot(collect(1:length(top)),top,"-p",label="$(round(filling,digits=3))")
	#legend()
	#ylabel("Long Range Correlation")
	#

	#=if idx > 1
		for i in 1:edge_length
			plot(fillings[idx-1:idx],[corrlengths[idx-1][i],corrlengths[idx][i]],"-p",c=col[i])
		end
	end=#
end
end
#
fig = figure()
xvals = fillings#num_parts ./ (fillings .* (edge_length^2))
plot(xvals,shortrange,"-p",label="Short Range")
plot(xvals,longrange,"-p",label="Long Range")
plot(xvals,sforderparams,"-p",label="Total")
xlabel("Filling")
ylabel("Zero Mom Occs")
legend()
#
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

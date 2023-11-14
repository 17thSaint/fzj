if false
include("long-range-ttn.jl")
include("fqh_effective.jl")
using PyPlot
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




#
params_dict = Dict([("if_periodic_phys",false),("if_periodic_synth",true),("nflavors",10),("L",16),("nbosons",10)])
loc = "/home/patrick/fzj/main-git/cluster-data/orsay-sept23"
all_files = find_data_file(params_dict,"mps","jld2",loc)
display(all_files)
#
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
plot([1/8 for i in 1:length(ys)],ys)
scatter(alphas,dbds)
#
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

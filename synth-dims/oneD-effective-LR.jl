include("long-range-ttn.jl")
include("fqh_effective.jl")
#using PyPlot

function fix_filling(L,nflavors,nu)
	prod = L * nu * nflavors
	for nbosons in L/2-1:L-2
		inv_alpha = round(prod/nbosons,digits=5)
		if isinteger(inv_alpha)
			println("Found Alpha = 1/$inv_alpha")
			return Int(nbosons),1/inv_alpha
		end
	end
	println("Not Found")
	return nothing,nothing
end

function make_1deff_filenamedict(model_paras)
	fdict = Dict([("Lphys",model_paras[:L]),("Lsynth",model_paras[:nflavors]),("nbosons",model_paras[:nbosons]),("alpha",model_paras[:alpha]),("hopping_anisotropy",model_paras[:hopping_anisotropy]),("if_periodic_phys",model_paras[:if_periodic_phys]),("if_periodic_synth",model_paras[:if_periodic_synth])])
	model_paras[:if_2ord_pert] ? fdict["U2"] = model_paras[:U2] : nothing
	model_paras[:if_nn_int] ? fdict["U1"] = model_paras[:U1] : nothing
	model_paras[:twist_angle] != [0.0,0.0] ? fdict["tw1"] = model_paras[:twist_angle][1] : nothing
	model_paras[:twist_angle] != [0.0,0.0] ? fdict["tw2"] = model_paras[:twist_angle][2] : nothing
	return fdict
end

function get_1deff_model_params(params_dict::Dict)

	# DMRG parameters
	nrgtol = get(params_dict, "nrgtol", 5E-5)
	minimum_sweeps = get(params_dict, "minimum_sweeps", 10)
	observer = TTNKit.DMRGObserver(;energy_tol=nrgtol,minsweeps=minimum_sweeps)
	cutoff = get(params_dict, "cutoff", 1E-8)
	syms = get(params_dict, "syms", true)
	nswps = get(params_dict, "num_sweeps", 100)
	noise = get(params_dict, "noise", [0.0])
	es_count = get(params_dict, "es_count", 0)


	# Lattice/TTN Parameters
	Lphys = get(params_dict, "Lphys", 4)
	Lsynth = get(params_dict, "Lsynth", 4)
	mdim = get(params_dict, "mdim", 300)
	if_periodic_phys = get(params_dict, "if_periodic_phys", false)
	if_periodic_synth = get(params_dict, "if_periodic_synth", false)
	num_particles = get(params_dict, "particles", 2)


	# Hamiltonian parameters
	if_remapping = get(params_dict, "if_remapping", if_periodic_phys)
	remapping = if_remapping ? remapping_nnn(Lphys) : collect(1:Lphys)
	alpha = get(params_dict, "alpha", nothing)
	flux_direction = get(params_dict,"flux_direction", "phys")
	if if_periodic_synth && !if_periodic_phys
        flux_direction = "synth"
    elseif !if_periodic_synth && if_periodic_phys
        flux_direction = "phys"
    end

	hopping_amplitude = get(params_dict, "ts", 1.0)
	anis = get(params_dict, "hopping_anisotropy", 1.0)

	if_2ord_pert = get(params_dict, "if_2ord_pert", false)
	if_NN = get(params_dict, "if_nn_int", false)
	if_NN ? U1 = get(params_dict, "U1", 1.0) : U1 = 0.0
	if_2ord_pert ? U2 = get(params_dict, "U2", 1.0) : U2 = 0.0

	mag_off = get(params_dict, "mag_off", false)
	centralflux_strength = get(params_dict, "centralflux_strength", 0.0)
	twist_angle = [get(params_dict, "tw1", 0.0),get(params_dict, "tw2", 0.0)]

	if isnothing(alpha)
		filling = get(params_dict, "filling", 0.5)
		phys_shift,synth_shift = !if_periodic_phys,!if_periodic_synth
		alpha = num_particles/(filling*(Lphys - phys_shift)*(Lsynth - synth_shift))
		filling == 0.0 ? alpha = 0.0 : nothing
		mag_off = false
	else
		mag_off = alpha == 0.0
	end
	flux_dir = get(params_dict,"flux_direction","phys")
    if if_periodic_synth && !if_periodic_phys
        flux_dir = "synth"
    elseif !if_periodic_synth && if_periodic_phys
        flux_dir = "phys"
    end
	if_check_fluxes = get(params_dict, "if_check_fluxes", true)
	if_check_fluxes ? flux_dir = check_fluxes(alpha,Lphys,Lsynth,if_periodic_phys,if_periodic_synth,flux_dir) : nothing


	# What to calculate
	if_densmat = get(params_dict, :if_densmat, true)
	save_data = get(params_dict, "if_save_data", true)
	if_cluster = any([occursin("local",pwd()),occursin("Local",pwd()),occursin("geraghty",pwd())])
	if_continuous_saving = false#get(params_dict,"if_continuous_saving",if_cluster || layer_count >= 7)
	save_data ? nothing : if_continuous_saving = false
	es_count = get(params_dict, "es_count", 0)

	if if_periodic_phys && if_periodic_synth
		dataloc = get_folder_location("cluster-data/synth-dims/torus")
	elseif if_periodic_phys || if_periodic_synth
		dataloc = get_folder_location("cluster-data/synth-dims")
	elseif !if_periodic_phys && !if_periodic_synth
		dataloc = get_folder_location("cluster-data/synth-dims/obc")
	end
	if es_count > 0
		dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
	end
	loc = get(params_dict, "dataloc", dataloc)


	# hardware parameters
	if_gpu = get(params_dict, "if_gpu", false)
	
	model_paras_dict = Dict("hopping_anisotropy"=>anis,
						"nbosons"=>num_particles,
						"L"=>Lphys,
						"nflavors"=>Lsynth,
						"remapping"=>remapping,
						"ts"=>hopping_amplitude,
						"syms"=>syms,
						"cutoff"=>cutoff,
						"observer"=>observer,
						"twist_angle"=>twist_angle,
						"if_continuous_saving"=>if_continuous_saving,
						"nrgtol"=>nrgtol,
						"if_densmat"=>if_densmat,
						"es_count"=>es_count,
						"centralflux_strength"=>centralflux_strength,
						"if_periodic_phys"=>if_periodic_phys,
						"if_periodic_synth"=>if_periodic_synth,
						"if_nn_int"=>if_NN,
						"if_2ord_pert"=>if_2ord_pert,
						"U1"=>U1,
						"U2"=>U2,
						"alpha"=>alpha,
						"flux_direction"=>flux_direction,
						"no_magF"=>mag_off,
						"if_gpu"=>if_gpu,
						"if_save_data"=>save_data,
						"mdim"=>mdim,
						"num_sweeps"=>nswps,
						"minimum_sweeps"=>minimum_sweeps,
						"phi"=>alpha,
						"noise"=>noise,
						"output_level"=>1,
						"location"=>loc)
	

	filename = join(["mps",make_parameters_filename(make_1deff_filenamedict(dict_to_symbols(model_paras_dict)))],"-")
	filename = check_plot_label(filename,"mps")
	model_paras_dict["name"] = filename
	if_find_data = get(params_dict, "if_find_data", true)
	if_exists, found_data = if_find_data ? check_data_exists(get_params_dict_from_filename(filename),"mps"; location=model_paras_dict["location"],output_level=false) : (false,nothing)
	return dict_to_symbols(model_paras_dict),found_data
end

function get_all_densities(Lmax; kwargs...)
	smallest_density = get(kwargs, :smallest_density, 0.25)
	smallest_sitecount = get(kwargs, :smallest_sitecount, 20)
	alpha_limit = get(kwargs, :alpha_limit, 0.35)
	nmin = get(kwargs, :nmin, 3)
	number_to_keep = get(kwargs, :number_to_keep, "all")

	edge_length_limit = Int(floor(sqrt(smallest_sitecount)))
	configurations = []
	oneDdensities = []
	for Lx in edge_length_limit:Lmax
		for Ly in edge_length_limit:Lmax
			if Lx*Ly < smallest_sitecount || Lx / Ly > 2 || Ly / Lx > 2
				continue
			end
			for n in nmin:Lx-1
				if n/Lx > smallest_density && n/(0.5 * Lx * Ly) < alpha_limit && !(n/Lx in oneDdensities) && typeof(check_fluxes(n/(0.5 * Lx * Ly),Lx,Ly,true,true,"phys",false)) == String
					push!(configurations,(Lx,Ly,n))
					push!(oneDdensities,n/Lx)
				end
			end
		end
	end

	if number_to_keep != "all"
		if number_to_keep > length(oneDdensities)
			println("Number to keep is greater than the number of configurations")
			return configurations
		elseif number_to_keep == length(oneDdensities)
			return configurations
		end
		keeping_indices = zeros(Int, number_to_keep)
		keeping_indices[end] = length(oneDdensities)
		hopeful_densities = range(minimum(oneDdensities),stop=maximum(oneDdensities),length=number_to_keep)[1:end-1]
		for (idx,hd) in enumerate(hopeful_densities)
			keeping_indices[idx] = abs.(oneDdensities .- hd) |> argmin
		end
		configurations = configurations[keeping_indices]
	end
	configurations = configurations[sortperm([cc[3]/cc[1] for cc in configurations])]
	return configurations
end

nev = 4
cols = ["b","g","r","m","c"]
if nev > length(cols)
	cols = repeat(cols,ceil(Int,nev/length(cols)))
end

open_cores = 5
if typeof(open_cores) != String
	BLAS.set_num_threads(open_cores)
	display(BLAS.get_config())
end

if false
	cd("../cluster-data/synth-dims/excited-states")
	allfiles = readdir()
	cd("../../../synth-dims")
	mpsfiles = allfiles[findall(x -> occursin("mps",x),allfiles)]
	for f in mpsfiles
			d,m = read_data_jld2("../cluster-data/synth-dims/excited-states/"*f; output_level=0)
			if m["L"] * m["nflavors"] < 10
				continue
			end
			if m["maxlinkdim"] == 400
				continue
			else
				println("MaxDim is ",m["maxlinkdim"]," while given is ",m["mdim"])
			end
			nrgs = Vector{Float64}(undef,3)
			nrgs[1] = abs(m["observer"].energies[end] - m["observer"].energies[end-1]) < m["nrgtol"] ? m["observer"].energies[end] : 10000.0
			nrgs[2] = abs(m["observer_1"].energies[end] - m["observer_1"].energies[end-1]) < m["nrgtol"] ? m["observer_1"].energies[end] : 10000.0
			nrgs[3] = "observer_2" in keys(m) ? m["observer_2"].energies[end] : 10000.0
			for i in 1:3
				change = 0.05
            	shift = 0.0#(i - 3/2) * ((0.1*change)/(3/2))
				scatter(m["nbosons"]/m["L"] + shift,(nrgs[i] - nrgs[1])/(1.0),c=cols[i])#m["L"]*m["nflavors"]
			end
	end
	xlabel("1D Density")
	ylabel("NRG - E0")
end


if true
		#tws = range(0.0,stop=1.0,length=2)
		#for tw1 in tws
		#for (idx,alpha) in enumerate(alphs)
		#denssets = [(8,4,4),(4,8,4),(6,8,4),(7,6,6),(8,5,5),(9,5,5)]
		#oneDdensities = [c[end]/c[1] for c in denssets]
		#for (lx,ly,n) in denssets
		#which_configs = get_all_densities(15,smallest_sitecount=100)
		#all_configs = get_all_densities(19,smallest_density=0.5,number_to_keep=30,smallest_sitecount=30)
		#params_dict = make_args_dict(ARGS)
		#which_configs = all_configs[5*params_dict["config_number"]+1:5*params_dict["config_number"]+5]
		#for (lx,ly,n) in which_configs
		#for n in 2:6
			lx,ly,n = 9,6,3#which_configs[i]#6,5,5
			#lx = Int(3*n)
			#ly = Int(2*n)
			params_dict = Dict([("Lphys",lx),("Lsynth",ly),("particles",n),("es_count",nev-1),("nrgtol",1e-6),("mdim",400),("if_periodic_phys",true),("if_periodic_synth",true),("filling",0.5),("if_find_data",false),("if_save_data",false)])
			params_dict["Lphys"] = lx
			params_dict["Lsynth"] = ly
			params_dict["particles"] = n
			model_paras,found_data = get_1deff_model_params(params_dict)
			if isnothing(found_data)
				if model_paras[:es_count] > 0
					psis,rhos,nrgs = excited_states_mps(model_paras[:es_count],model_paras[:phi],model_paras[:L],model_paras[:nflavors],model_paras[:nbosons]; model_paras...,metadata=named_tuple_to_dict(model_paras))
				else
					psi,rho,nrg = execute_mps(model_paras[:phi],model_paras[:L],model_paras[:nflavors],model_paras[:nbosons]; model_paras...,metadata=named_tuple_to_dict(model_paras))
				end
			else-
				#psis = Vector{MPS}(undef,found_data[2]["es_count"]+1)
				#rhos = Vector{Array}(undef,found_data[2]["es_count"]+1)
				nrgs = Vector{Float64}(undef,found_data[2]["es_count"]+1)
				#psis[1] = "mps" in keys(found_data[1]) ? found_data[1]["mps"] : nothing
				#rhos[1] = found_data[1]["densmat"]
				nrgs[1] = found_data[2]["observer"].energies[end]
				for i in 1:found_data[2]["es_count"]
					#psis[i+1] = "mps_$i" in keys(found_data[1]) ? found_data[1]["mps_$(i)"] : nothing
					#rhos[i+1] = found_data[1]["densmat_$(i)"]
					nrgs[i+1] = found_data[2]["observer_$(i)"].energies[end]
				end
			end
			for i in 1:nev
				get_occupancy(psis[i]; remapping=model_paras[:remapping],plot_title="E$i=$(round(nrgs[i],digits=4)) N = $n")
			end
			#s1occ = get_occupancy(psis[1]; if_3d=true, remapping=model_paras[:remapping],plot_title="E1 Phys=$(params_dict["if_periodic_phys"]), Synth=$(params_dict["if_periodic_synth"]) N = $n")
			#s2occ = get_occupancy(psis[2]; if_3d=true, remapping=model_paras[:remapping],plot_title="E2 Phys=$(params_dict["if_periodic_phys"]), Synth=$(params_dict["if_periodic_synth"]) N = $n")

			#plot(model_paras[:observer].energies .- model_paras[:observer].energies[end],label="$(params_dict["if_remapping"])")
			#yscale("log")
			#legend()
			#=xxs = tws
			for i in 1:model_paras[:es_count]+1
				#=change = abs(xxs[1] - xxs[2])
				xval = t
				shift = (i - model_paras[:es_count]/2) * ((0.1*change)/(model_paras[:es_count]/2))
				scatter(xval + shift,nrgs[i],c=cols[i])=#
				scatter(n/lx,nrgs[i] - nrgs[1],c=cols[i])
			end
			xlabel("1D Density")
			ylabel("NRG - E0")=#
		#end
end

if false
	save_nothing = true
	params_dict = Dict()
	#L = 24#get(params_dict, "L", 4)
	#nbosons = 5#get(params_dict, "nbosons", nflavors)
	nflavors = 5#get(params_dict, "nflavors", Int(L/2))
	t1 = get(params_dict, "t1", 1.0)
	t2 = get(params_dict, "t2", 1.0)
	U = get(params_dict, "U", 100)
	U1 = 4*t1^2/U
	U2 = U1/2
	conserve_qns = true
	if_nn_int = false#get(params_dict, "if_nn_int", false)
	if_2ord_pert = false#get(params_dict, "if_2ord_pert", false)
	nsweeps = 100
	mdim = get(params_dict, "mdim", 60)
	noises = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0.0]
	if_save_data = save_nothing ? false : true
	data_loc = "/home/patrick/fzj/main-git/cluster-data/orsay-sept23"
	if_periodic_phys = false
	if_periodic_synth = false
	#nflavors = 9
	#alpha = 23/(24^2)


	dmrg_obs = TTNKit.DMRGObserver(;energy_tol=10^-3,minsweeps=6)

	other_params_dict = Dict([("U",U),("conserve_qns",conserve_qns),("nsweeps",nsweeps),("mdim",mdim),("noise",noises)])
	savefig_data = save_nothing ? false : true
	savefig = save_nothing ? false : true
	if_lines = false

	#nbosons,alpha = fix_filling(L,nflavors,1/2)
	#alpha = 0.0

	density = 5/40
	#Ls = [8]
	L = 1
	count = 5
	alphastart = 0.081
	alphaend = 0.3
	alphas = [0.0]#[alphastart + (i-1)*(alphaend - alphastart)/(count-1) for i in 1:count] .- 0.0001
	#alphas = [alphas; alphas .+ 0.0001]
	display(alphas)
	wavefuncs = []
	rhos = []
	#nbosons = Int(L/2)
	#fillings = ["1/2","2/3","1/3"]
	for (idx,alpha) in enumerate(alphas)
		nbosons = 1#Int(L*nflavors*density)
		#alpha = nbosons/((L-1)*nflavors)
		filename_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",nbosons),("alpha",round(alpha,digits=4)),("if_periodic_synth",if_periodic_synth),("if_periodic_phys",if_periodic_phys)])
		#filename_dict_highdens = Dict([("L",L),("nflavors",nflavors),("nbosons",nbosons_highdens),("alpha",round(alpha,digits=4)),("if_nn_int",if_nn_int),("if_2ord_pert",if_2ord_pert),("if_periodic",if_periodic)])

		datafile_name = make_parameters_filename(filename_dict)
		println(datafile_name)
		#datafile_name_highdens = make_parameters_filename(filename_dict_highdens)

		model_paras = (t1 = t1, t2 = t2, phi = phi, U1 = U1, U2 = U2, L = L, nflavors = nflavors, nbosons = nbosons, if_nn_int = if_nn_int, if_2ord_pert = if_2ord_pert, mdim = mdim, nsweeps = nsweeps, noise = noises, if_save_data = if_save_data, location = data_loc, if_periodic_synth=if_periodic_synth, if_periodic_phys=if_periodic_phys, name=datafile_name, observer = dmrg_obs)

		#model_paras_highdens = (t1 = t1, t2 = t2, phi = phi, U1 = U1, U2 = U2, L = L, nflavors = nflavors, nbosons = nbosons_highdens, if_nn_int = if_nn_int, if_2ord_pert = if_2ord_pert, mdim = mdim, nsweeps = nsweeps, noise = noise, if_save_data = if_save_data, location = data_loc, if_periodic=if_periodic, name=datafile_name_highdens)

		#metadata_dict_highdens = merge(named_tuple_to_dict(model_paras_highdens),other_params_dict)
		metadata_dict = merge(named_tuple_to_dict(model_paras),other_params_dict)

		psi = execute_mps(U1,U2,phi,L,nflavors,nbosons; model_paras...,metadata=metadata_dict)
		append!(wavefuncs,[psi])
		get_occupancy(psi; plot_title="$alpha")
		get_greenfunc(psi)
		get_greenfunc(psi,"phys")
		#densmat = correlation_matrix(psi,"FullDag","FullHat") #./ 2.0
		#append!(rhos,[densmat])
		#=if false
		fig = figure()
		imshow(real.(densmat))
		colorbar()
		title("Virt Dim = $(round(alpha,digits=3))")
		end=#
		#append!(all_wavefuncs,[psi])
		#psi_highdens = execute_mps(U1,U2,phi,L,nflavors,nbosons_highdens; model_paras_highdens...,metadata=metadata_dict_highdens)
		#mrez = momentum_occupation(psi,nbosons,100,10.0; model_paras...,plot_title="Virt Dim = $nflavors",if_log=false,if_plot=true)
		#append!(mom_occs,[mrez[2]])
	end

end




#title("Nflavors = $nflavors")

#=end

fig = figure()
momenta = [0.0 + (i-1)*(10.0 - 0.0)/(200-1) for i in 1:200]
for i in 1:length(mom_occs)
	plot(momenta./pi,mom_occs[i]./nbosons,label="$(Ls[i])")
end
title("Momentum Distribution range PhysDim NF = $nflavors, w/ 2nd Order")
xlabel("p/pi")
ylabel("Occupation / nbosons")
legend()

append!(allall_wavefuncs,[all_wavefuncs])
append!(all_mom_occs,[mom_occs])

end
=#
#=
for i in 1:length(fillings)
	filling = fillings[i]
	nbosons,alpha = fix_filling(L,nflavors,filling)
	isnothing(nbosons) ? continue : nothing
	#alpha = alphas[i]

	phi  = 2π*alpha
	#phi_right  = 2π*(alpha+change)
	#phi_left  = 2π*(alpha-change)
	#filling = nbosons / (alpha * L * nflavors)
	println("Nu = $filling with Alpha = $alpha and NBosons = $nbosons")
	filename_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",nbosons),("alpha",round(alpha,digits=4)),("if_nn_int",if_nn_int),("if_2ord_pert",if_2ord_pert)])
	datafile_name = make_parameters_filename(filename_dict)
	model_paras = (t1 = t1, t2 = t2, phi = phi, U1 = U1, U2 = U2, L = L, nflavors = nflavors, nbosons = nbosons, if_nn_int = if_nn_int, if_2ord_pert = if_2ord_pert, mdim = mdim, nsweeps = nsweeps, noise = noise, if_save_data = if_save_data, location = data_loc, name = datafile_name)
	metadata_dict = merge(named_tuple_to_dict(model_paras),other_params_dict)
	println(datafile_name)
	if true
	psi = execute_mps(U1,U2,phi,L,nflavors,nbosons; model_paras...,metadata=metadata_dict)
	append!(all_wavefuncs,[psi])
	#println("Done Center")
	
	
	#model_paras_left = copy(model_paras)
	#model_paras_left[:phi] = phi_left
	#psi_left = execute_mps(U1,U2,phi_left,L,nflavors; model_paras...,psi_guess = psi, if_save_data=false)
	#psi_right = execute_mps(U1,U2,phi_right,L,nflavors; model_paras...,psi_guess = psi, if_save_data=false)
	plottitle = "Filling = $filling"
	rez51 = get_occupancy(psi;if_plot=true, if_save_data=savefig_data,name=datafile_name,location=data_loc,metadata=metadata_dict,if_save_fig=savefig,plot_title=plottitle)
	rez11 = get_greenfunc(psi,"virt"; if_plot=true,if_lines=if_lines, if_save_data=savefig_data,name=datafile_name,location=data_loc,metadata=metadata_dict,if_save_fig=savefig,plot_title=plottitle)
	rez21 = get_greenfunc(psi,"phys"; if_plot=true,if_lines=if_lines, if_save_data=savefig_data,name=datafile_name,location=data_loc,metadata=metadata_dict,if_save_fig=savefig,plot_title=plottitle)
	end
	
	#=
	bd_l = bulk_density(psi_left)
	append!(bds,[bd_l])
	bd_c = bulk_density(psi)
	append!(bds,[bd_c])
	bd_r = bulk_density(psi_right)
	append!(bds,[bd_r])
	
	dbd = (bd_r - bd_l) / (2*change)
	println(alpha,", ",dbd)
	dbds[i] = dbd
	=#
	
end
#=
rez5 = get_occupancy(all_wavefuncs[1];)
rez = get_greenfunc(all_wavefuncs[1],"virt"; if_plot=true,if_lines=if_lines)#, if_save_data=savefig_data,name=datafile_name,location=data_loc,metadata=metadata_dict,if_save_fig=savefig)
rez2 = get_greenfunc(all_wavefuncs[1],"phys"; if_plot=true,if_lines=if_lines)#, if_save_data=savefig_data,name=datafile_name,location=data_loc,metadata=metadata_dict,if_save_fig=savefig)
#re3 = get_greenfunc(all_wavefuncs[1],"virt"; if_plot=true,rev=true,if_lines=if_lines)
#rez4 = get_greenfunc(all_wavefuncs[1],"phys"; if_plot=true,rev=true,if_lines=if_lines)
=#
#println(derivs,", ",errs)
#for i in 1:2
#errorbar(alphas,derivs[i,:],yerr=[errs[i,:],errs[i,:]],label="$i")
#end
#legend()



=#
























"fin"

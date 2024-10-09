include("long-range-ttn.jl")
include("fqh_effective.jl")
include("plottings.jl")
include("hatsugai-mbcn.jl")
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
	output_level = get(params_dict, "output_level", 1)


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
	if_check_fluxes ? flux_dir = check_fluxes(alpha,Lphys,Lsynth,if_periodic_phys,if_periodic_synth,flux_dir; output_level=output_level,if_ed=false) : nothing


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
	if twist_angle != [0.0,0.0]
		dataloc = get_folder_location("cluster-data/synth-dims/twists")
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
						"output_level"=>output_level,
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

function run_normal_1deffmps(params_dict::Dict; kwargs...)
	model_paras,found_data = get_1deff_model_params(params_dict)
	if isnothing(found_data)
		if model_paras[:es_count] > 0
			psis,rhos,nrgs = excited_states_mps(model_paras[:es_count],model_paras[:phi],model_paras[:L],model_paras[:nflavors],model_paras[:nbosons]; model_paras...,metadata=named_tuple_to_dict(model_paras))
		else
			psis,rhos,nrgs = execute_mps(model_paras[:phi],model_paras[:L],model_paras[:nflavors],model_paras[:nbosons]; model_paras...,metadata=named_tuple_to_dict(model_paras))
		end
	else
		psis = Vector{MPS}(undef,found_data[2]["es_count"]+1)
		rhos = Vector{Array}(undef,found_data[2]["es_count"]+1)
		nrgs = Vector{Float64}(undef,found_data[2]["es_count"]+1)
		psis[1] = "mps" in keys(found_data[1]) ? found_data[1]["mps"] : nothing
		rhos[1] = found_data[1]["densmat"]
		nrgs[1] = found_data[2]["observer"].energies[end]
		for i in 1:found_data[2]["es_count"]
			psis[i+1] = "mps_$i" in keys(found_data[1]) ? found_data[1]["mps_$(i)"] : nothing
			rhos[i+1] = found_data[1]["densmat_$(i)"]
			nrgs[i+1] = found_data[2]["observer_$(i)"].energies[end]
		end
	end

	return psis,rhos,nrgs,model_paras
end


nev = 3
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
	lx,ly,n = 6,6,6
	#ref_multiplets,rm1_name,rm2_name = get_reference_multiplets(lx,ly,n)
	#tws = range(0.0,stop=1.0,length=10)
	#g1s = zeros(ComplexF64,length(tws),length(tws))
	#g2s = zeros(ComplexF64,length(tws),length(tws))
	#oms = zeros(ComplexF64,length(tws),length(tws))
	#for (idx,tw1) in enumerate(tws)
	#for (idx2,tw2) in enumerate(tws)
	tw2 = 0.0
	tw1 = 0.0
		
		params_dict = Dict([("Lphys",lx),("Lsynth",ly),("particles",n),("tw1",tw1),("tw2",tw2),("if_remapping",false),("es_count",nev-1),("nrgtol",1e-6),("mdim",200),("if_periodic_phys",true),("if_periodic_synth",true),("filling",0.5),("if_find_data",false),("if_save_data",false)])
		psis,rhos,nrgs,model_paras = run_normal_1deffmps(params_dict)

		get_occupancy(psis[1]; plot_title=" E1 Lx = $lx, Ly = $ly, n = $n", remapping=model_paras[:remapping])
		get_occupancy(psis[2]; plot_title=" E2 Lx = $lx, Ly = $ly, n = $n", remapping=model_paras[:remapping])

		#plot_spectrum(tws,nrgs,idx,nev,"Theta_x / 2pi",false; plot_title=" Lx = $lx, Ly = $ly, n = $n, tw2 = $tw2")

		#=g1,g2,om = get_hatsugaifull(psis[1],psis[2],ref_multiplets)
		g1s[idx,idx2] = g1
		g2s[idx,idx2] = g2
		oms[idx,idx2] = om=#


	#end
	#end

	#plot_omega(tws,tws,oms; plot_title=" Lx = $lx, Ly = $ly, n = $n",if_mag=true)
	#plot_gamma(tws,tws,g1s,1; plot_title=" Lx = $lx, Ly = $ly, n = $n")
	#plot_gamma(tws,tws,g2s,2; plot_title=" Lx = $lx, Ly = $ly, n = $n")


end


























"fin"

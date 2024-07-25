include("fqh_effective.jl")
include("long-range-ttn.jl")
using PyPlot

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
	
	return dict_to_symbols(model_paras_dict)
end

nev = 3
cols = ["b","g","r","m","c"]
if nev > length(cols)
	cols = repeat(cols,ceil(Int,nev/length(cols)))
end

if true
		#tws = range(0.0,stop=1.0,length=2)
		#for tw1 in tws
		#for (idx,alpha) in enumerate(alphs)
		enssets = [(8,4,4),(4,8,4),(6,8,4),(7,6,6),(8,5,5),(9,5,5)]
		oneDdensities = [c[end]/c[1] for c in denssets]
		for (lx,ly,n) in denssets
				params_dict = Dict([("Lphys",lx),("Lsynth",ly),("particles",n),("es_count",nev),("nrgtol",1e-6),("mdim",200),("if_periodic_phys",true),("if_periodic_synth",true),("filling",0.5),("if_save_data",true)])
				model_paras = get_1deff_model_params(params_dict)
				if model_paras[:es_count] > 0
					psis,rhos,nrgs = excited_states_mps(model_paras[:es_count],model_paras[:phi],model_paras[:L],model_paras[:nflavors],model_paras[:nbosons]; model_paras...,metadata=named_tuple_to_dict(model_paras))
				else
					psi,rho,nrg = execute_mps(model_paras[:phi],model_paras[:L],model_paras[:nflavors],model_paras[:nbosons]; model_paras...,metadata=named_tuple_to_dict(model_paras))
				end
			#xxs = tws
			for i in 1:model_paras[:es_count]+1
				#=change = abs(xxs[1] - xxs[2])
				xval = t
				shift = (i - model_paras[:es_count]/2) * ((0.1*change)/(model_paras[:es_count]/2))
				scatter(xval + shift,nrgs[i],c=cols[i])=#
				scatter(n/lx,nrgs[i] - nrgs[1],c=cols[i])
			end
			xlabel("1D Density")
			ylabel("NRG - E0")
		end
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

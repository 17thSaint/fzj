### A Pluto.jl notebook ###
# v0.19.26

using Markdown
#using InteractiveUtils

# ╔═╡ 07f7c480-198f-11ee-28f5-65891e30bb85
using ITensors

# ╔═╡ 2c9a6330-d820-4637-bbe6-9ce40dbbbb3f
using LinearAlgebra

# ╔═╡ e0c40da6-0ffc-40ae-a31d-f51f6f792dac
#using Plots

# ╔═╡ 641ecd99-2fcc-49a7-8048-6203f8767589
md"
``\ket{0} \Rightarrow \hat{N}\ket{0} = 0``
``\ket{n_S = 1, S = 1,\dots, N_S} \Rightarrow \hat{N} \ket{n_S = 1, S = 1,\dots, N_S} = \ket{n_S = 1, S = 1,\dots, N_S}``
"

# ╔═╡ 2b12ea80-dce7-4b91-b902-9a730516bf49
function ITensors.space(::SiteType"ExtendedHardcore"; nflavors = 1, conserve_qns = false, qnname = "Nbosons")
	if conserve_qns
		return [QN(qnname, 0, -1) => 1, QN(qnname, 1, -1) => nflavors]
	end
	return 1 + nflavors
end

# ╔═╡ 2a6ecb0e-5075-425a-8cc7-fc79423b6abf
begin
	#ITensors.val(::ValName"Vac", s::SiteType"ExtendedHardcore") = 1
	
	function ITensors.val(::ValName{N}, s::SiteType"ExtendedHardcore") where{N}
	    return parse(Int, String(N)) + 1
	end

	
	#ITensors.state(::ValName"Vac", ::SiteType"ExtendedHardcore", s::Index) = vcat(1.0, zeros(dim(s)))
	
	function ITensors.state(::StateName{N}, ::SiteType"ExtendedHardcore", s::Index) where{N}
	    n  = parse(Int, String(N))
	    nd = dim(s)
	    st = zeros(nd)
	
	    st[n + 1] = 1.0
	    return itensor(st, s)
	end
	
end

# ╔═╡ f8382eef-be2a-4846-b378-220face8dc1b
parse_opname(::Type{<:OpName{S}}) where S = string(S)

# ╔═╡ 15f9a000-10a8-448b-b4aa-36d83b208826
function opstring_to_flavor(on::Type{<:OpName})
	on_str = parse_opname((on))

	sup_op_list = ["Cr", "Anh", "Ns", "Mx", "My", "Mz"]
	for sup_op in sup_op_list
		occursin(sup_op, on_str) && return sup_op, parse(Int64, SubString(on_str, (length(sup_op)+1):length(on_str)))
	end

	return on_str, nothing
end

# ╔═╡ bbbb808d-1f35-414f-a429-e6ea5f885512
begin
	function _op(::OpName"Anh", ::SiteType"ExtendedHardcore", d::Int; flavor)
		@assert flavor > 0
	    mat = zeros(d,d)
		mat[1, flavor + 1] = 1 + 0.0*im
		return mat
	end

	function _op(::OpName"Cr", ::SiteType"ExtendedHardcore", d::Int; flavor)
		@assert flavor > 0
	    mat = zeros(d,d)
		mat[flavor + 1,1] = 1 + 0.0*im
		return mat
	end
	
	function _op(::OpName"FullDag", ::SiteType"ExtendedHardcore", d::Int)
	    mat = zeros(d,d)
	    for i in 1:d-1
	    	mat[i + 1,1] = 1 + 0.0*im
	    end
	    return mat
	end
	
	function _op(::OpName"FullHat", ::SiteType"ExtendedHardcore", d::Int)
	    mat = zeros(d,d)
	    for i in 1:d-1
	    	mat[1,i + 1] = 1 + 0.0*im
	    end
	    return mat
	end

	_op(::OpName"I", ::SiteType"ExtendedHardcore", d::Int) = I(d)
	_op(::OpName"Id", ::SiteType"ExtendedHardcore", d::Int) = I(d)
	
	function _op(::OpName"N", ::SiteType"ExtendedHardcore", d::Int)
		mat = zeros(d,d)
		for j in 2:d
			mat[j,j] = 1 + 0.0*im
		end
		return mat
	end
	function _op(::OpName"Ns", ::SiteType"ExtendedHardcore", d::Int; flavor)
		mat = zeros(d,d)
		mat[flavor+1,flavor+1] = 1 + 0.0*im
		
		return mat
	end

	function _op(::OpName"S+", ::SiteType"ExtendedHardcore", d::Int)
		mat = zeros(d,d)
		for j in 2:d-1
			mat[j+1,j] = 1 + 0.0*im
		end
		return mat
	end
	function _op(::OpName"S-", ::SiteType"ExtendedHardcore", d::Int)
		mat = zeros(d,d)
		for j in 2:d-1
			mat[j, j+1] = 1 + 0.0*im
		end
		return mat
	end
	
	sx_12 = 0.5 .* [0 1; 1 0]
	sy_12 = (1/(2*1)) .* [0 1; -1 0]
	sz_12 = (1/(2*1)) .* [1 0; 0 -1]
			
	sx_1 = (1/sqrt(2)) .* [0 1 0; 1 0 1; 0 1 0]
	sy_1 = (1/(sqrt(2)*1)) .* [0 1 0; -1 0 1; 0 -1 0]
	sz_1 = (1/1) .* [1 0 0; 0 0 0; 0 0 -1]
		
	sx_32 = (1/2) .* [0 sqrt(3) 0 0; sqrt(3) 0 2 0; 0 2 0 sqrt(3); 0 0 sqrt(3) 0]
	sy_32 = (1/(2*1)) .* [0 sqrt(3) 0 0; -sqrt(3) 0 2 0; 0 -2 0 sqrt(3); 0 0 -sqrt(3) 0]
	sz_32 =  [3/2 0 0 0; 0 1/2 0 0; 0 0 -1/2 0; 0 0 0 -3/2]
		
	sx_2 = 0.5 .* [0 2 0 0 0; 2 0 sqrt(6) 0 0; 0 sqrt(6) 0 sqrt(6) 0; 0 0 sqrt(6) 0 2; 0 0 0 2 0]
	sy_2 = (0.5*1) .* [0 2 0 0 0; -2 0 sqrt(6) 0 0; 0 -sqrt(6) 0 sqrt(6) 0; 0 0 -sqrt(6) 0 2; 0 0 0 -2 0]
	sz_2 = [2 0 0 0 0; 0 1 0 0 0; 0 0 0 0 0; 0 0 0 -1 0; 0 0 0 0 -2]
	
	sx_52 = 0.5 .* [0 sqrt(5) 0 0 0 0; sqrt(5) 0 sqrt(8) 0 0 0; 0 sqrt(8) 0 sqrt(9) 0 0; 0 0 sqrt(9) 0 sqrt(8) 0 ;0 0 0 sqrt(8) 0 sqrt(5); 0 0 0 0 sqrt(5) 0]
	sy_52 = (0.5) .* [0 sqrt(5) 0 0 0 0; -sqrt(5) 0 sqrt(8) 0 0 0; 0 -sqrt(8) 0 sqrt(9) 0 0; 0 0 -sqrt(9) 0 sqrt(8) 0 ; 0 0 0 -sqrt(8) 0 sqrt(5); 0 0 0 0 -sqrt(5) 0]
	sz_52 = [5/2 0 0 0 0 0; 0 3/2 0 0 0 0; 0 0 1/2 0 0 0; 0 0 0 -1/2 0 0; 0 0 0 0 -3/2 0; 0 0 0 0 0 -5/2]
	
	spin_matrices_dict = Dict([("X,0.5",sx_12),("X,1.0",sx_1),("X,1.5",sx_32),("X,2.0",sx_2),("X,2.5",sx_52),("Y,0.5",sy_12),("Y,1.0",sy_1),("Y,1.5",sy_32),("Y,2.0",sy_2),("Y,2.5",sy_52),("Z,0.5",sz_12),("Z,1.0",sz_1),("Z,1.5",sz_32),("Z,2.0",sz_2),("Z,2.5",sz_52)])
	
	
	function _op(::OpName"Mz", ::SiteType"ExtendedHardcore", d::Int; flavor)
		spin = round((d-1)/2,digits=1)
		flavor1 = Int(round(flavor/10,digits=0))
		flavor2 = Int(flavor - 10*flavor1)
		spin_part = spin_matrices_dict["Z,"*string(spin)]
		anh_mat = zeros(d,d)
		anh_mat[1, flavor2 + 1] = 1
		cr_mat = zeros(d,d)
		cr_mat[flavor1 + 1,1] = 1
		return spin_part .* (cr_mat * anh_mat)
	end
	
	
	function ITensors.op(opn::OpName, st::SiteType"ExtendedHardcore", d::Int)
		op, flavor = opstring_to_flavor(typeof(opn))
		isnothing(flavor) && (return _op(opn, st, d))
		op = OpName{Symbol(op)}()
		_op(op, st, d; flavor = flavor)
		
	end

	function ITensors.op(on::OpName, st::SiteType"ExtendedHardcore", s1::Index, s_tail::Index...; kwargs...)
	  rs = reverse((s1, s_tail...))
	  ds = dim.(rs)
	  opmat = op(on, st, ds...; kwargs...)
	  return itensor(opmat, prime.(rs)..., dag.(rs)...)
	end

	
end

# ╔═╡ 615e311e-b3bd-4cd3-9272-64ba2c806787
md"
``\sum_s b_{s, j}^\dagger b_{s, j+1} = \sum_{s,k} b_{s, j}^\dagger 
	\delta_{s,k}b_{k, j+1}``

``
\sum_{s,k} b_{s, j}^\dagger b_{k, j+1}
``
"

# ╔═╡ fa55709d-a20c-47ee-9802-8ea750588c63
function make_states(L::Int64,nbosons::Int64,nflavors::Int64)
	states = fill("0", L)

	jvisit = Int64[]
	for n in 1:nbosons
		jnext = rand(1:L)
		while jnext ∈ jvisit
			jnext = rand(1:L)
		end
		states[jnext] = string(rand(1:nflavors))
		push!(jvisit, jnext)
	end
	return states
end

function make_states(L::Int64,nflavors::Int64,organization::Vector)
	states = fill("0", L)
	for (i,j) in organization
		states[i] = string(j)
	end
	return states
end

function make_wavefunc(L::Int64,nflavors::Int64,organization::Vector; kwargs...)
	conserve_qns = get(kwargs, :conserve_qns, true)
	sidx = siteinds("ExtendedHardcore", L; conserve_qns = conserve_qns, nflavors = nflavors)
	states = make_states(L,nflavors,organization)
	psi0 = randomMPS(sidx, states)
	return psi0
end

function remapping_nnn(L) # this function remaps the periodic in physical dimension hamiltonian to nearest neighbor hopping only for the physical dimension
    remap = Int.(ones(L))
    shift = L % 2 == 0.0 ? 0 : -1
    for i in 1:Int(ceil(L/2))
         remap[i] = 2*i-1
    end
    remap[Int(ceil(L/2))+1] = L + shift
    for i in Int(ceil(L/2))+2:L
        remap[i] = remap[i-1] - 2
    end

    return remap
end

function make_vacuum(L,nflavors; kwargs...)
	sitetype = get(kwargs, :sitetype, "ExtendedHardcore")
	conserve_qns = get(kwargs, :conserve_qns, true)
	return randomMPS(siteinds("ExtendedHardcore", L; conserve_qns = conserve_qns, nflavors = nflavors), ["0" for i in 1:L])
end

function equal_weight(L,nflavors; kwargs...)
	return ones(nflavors,L) ./ (L * nflavors)
end

function make_particle(L,nflavors,weight_function)

	weights = weight_function(L,nflavors)
	make_parts = OpSum()
	for j in 1:L
		for s in 1:nflavors-1
			weight = weights[s,j]
			make_parts += (weight,"Cr$s",j)
		end
	end
	return make_parts
end

function initialize_mps(psi::MPS,particle_count::Int; weight_function = equal_weight, kwargs...)

	L = length(psi)
	nflavors = size(psi[1])[1]
	sumpart = make_particle(L,nflavors, weight_function)
	creation_mpo = MPO(sumpart,siteinds(psi))
	if particle_count > L
		println("Too many particles: L=$L, nbosons=$particle_count")
		return nothing
	end
	for i in 1:particle_count
		psi = apply(creation_mpo,psi; cutoff=1E-10, kwargs...)
	end
	return psi
end

function hamiltonian(t1, t2, phi, U1, U2, L, nflavors; kwargs...)
	if_nn_int = get(kwargs, :if_nn_int, true)
	if_2ord_pert = get(kwargs, :if_2ord_pert, true)
	if_periodic_phys = get(kwargs, :if_periodic_phys, false)
	if_periodic_synth = get(kwargs, :if_periodic_synth, false)
	if_tilt = get(kwargs, :if_tilt, false)
	tilt_strength = get(kwargs, :tilt_strength, 0.0)
	if_current = get(kwargs, :if_current, false)
	current_strength = get(kwargs, :current_strength, 0.0)
	
	ampo = OpSum()
	for j in 1:L
		for s in 1:nflavors
			# physical dimension hopping
			next_site = j+1
			if j == L
				if if_periodic_phys
					next_site = 1
				else
					continue
				end
			end
			ampo += (-t1 * exp(im*phi*s*1.0), "Cr$s", j, "Anh$s", next_site)
			ampo += (-t1 * exp(-im*phi*s*1.0), "Anh$s", j, "Cr$s", next_site)
		end

		# attractive physical nearest neighbor density interaction
		if if_nn_int
			next_site = j+1
			if j == L
				if if_periodic_phys
					next_site = 1
				else
					continue
				end
			end
			for s in nflavors
				for k in nflavors
					ampo += (-U1/2, "Ns$(s)", j, "Ns$(k)", next_site)
				end
			end
		end
	end
	
	if if_2ord_pert
		for j in 1:L-2
			for s in nflavors
				for k in nflavors
					ampo += (-U2, "Cr$(k)", j, "Ns$(s)", j+1, "Anh$(k)", j+2)
					ampo += (-U2, "Anh$(k)", j, "Ns$(s)", j+1, "Cr$(k)", j+2)
				end
			end
		end
	end
	
	for j in 1:L
		for s in 1:nflavors
			# synthetic dimension hopping
			next_site = s+1
			if s == nflavors
				if if_periodic_synth
					next_site = 1
				else
					continue
				end
			end
			ampo += (-t2 * exp(im*phi*j*0.0), "Cr$(next_site) * Anh$(s)", j)
			ampo += (-t2 * exp(-im*phi*j*0.0), "Cr$(s) * Anh$(next_site)", j)
			#ampo += (-t2 * exp(-im*phi*j), "Anh$(next_site) * Cr$(s)", j)
		end
	end

	if if_tilt
		for j in 1:L
			for s in 1:nflavors
				ampo += (-tilt_strength*j, "Ns$(s)", j)
			end
		end
	end
	
	return ampo
end

function run_again(filename; kwargs...)
	location = get(kwargs, :location, pwd())
	if_densmat = get(kwargs, :if_densmat, true)
	nrgvar_tol = get(kwargs, :nrgvar_tol, 10^-7)

	data,metadata = read_data_jld2(filename,location)
	
	obs = NRGVarObserver(nrgvar_tol,metadata["ham"])
	new_metadata = merge(metadata,Dict([("observer",obs),("nrgvar_tol",nrgvar_tol),("psi_guess",data["mps"])]))

	new_gs,new_densmat = execute_mps(nothing,nothing,metadata["chi"],metadata["L"],metadata["nflavors"],metadata["nbosons"]; dict_to_symbols(new_metadata)...,if_densmat=if_densmat,metadata=new_metadata,mdim=maximum(metadata["mdim"]),running_again=true)
	println("Energy Variance = ",energy_variance(new_gs,metadata["ham"]))
end

function get_occupancy(wavefunc::MPS; kwargs...)
	L,nflavors = get_mps_dims(wavefunc)
	if_squared = get(kwargs, :if_squared, false)
	if_remapping = get(kwargs, :if_remapping, false)

	remap = if_remapping ? remapping_nnn(L) : nothing
	
	if_plot = get(kwargs, :if_plot, true)
	if_save_data = get(kwargs, :if_save_data, false)
	if if_save_data
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "occs")
		filename = check_plot_label(filename,"occs")
		metadata = get(kwargs, :metadata, nothing)
	end
		
	occ_mat = zeros(L,nflavors)
	for s in 1:nflavors
		loc_op = "Ns$(s)"
		if if_squared
			loc_op = "Ns$(s) * Ns$(s)"
		end
		if if_remapping
			linearoccs = expect(wavefunc, loc_op)
			for i in 1:L
				occ_mat[i, s] = linearoccs[remap[i]] 
			end
		else
			occ_mat[:, s] = expect(wavefunc, loc_op)
		end
	end
	if_plot ? mps_plot_occupancy(occ_mat,L,nflavors; kwargs...) : nothing
	data_dict = Dict([("vals",occ_mat)])
	if_save_data ? write_data_jld2(filename,data_dict,location,metadata) : nothing
	
	return occ_mat
end

function execute_mps(U1,U2,phi,L,nflavors,nbosons; kwargs...)
	running_again = get(kwargs, :running_again, false)
	psi_ortho = get(kwargs, :psi_ortho, nothing)
	opl = get(kwargs, :outputlevel, 1)
	conserve_qns = get(kwargs, :conserve_qns, false)
	nsweeps = get(kwargs, :nsweeps, 500)
	psi0 = get(kwargs, :psi_guess, nothing)
	if_parton = get(kwargs, :if_parton, true)
	ham = get(kwargs, :ham, nothing)
	mdim = get(kwargs, :mdim, 100)
	if mdim >= 100 && !running_again
		if mdim >= 400
			mdim = [Int(floor(mdim/4)),Int(floor(mdim/4)),Int(floor(mdim/4)),Int(floor(mdim/2)),Int(floor(mdim/2)),Int(floor(mdim/2)),Int(floor(mdim/2)),Int(floor(mdim/2)), mdim]
		else
			mdim = [Int(floor(mdim/4)),Int(floor(mdim/2)), mdim]
		end
	end
	noise = get(kwargs, :noise, 0.0)
	obs = get(kwargs, :observer, NoObserver())
	if_save_data = get(kwargs, :if_save_data, false)
	if_gpu = get(kwargs, :if_gpu, false)
	t1 = 1.0
	t2 = 1.0
	if_nrg = get(kwargs, :if_nrg, false)
	if_densmat = get(kwargs, :if_densmat, true)

	if if_parton
		psi0 = make_vacuum(L,nflavors; kwargs...)
		psi0 = initialize_mps(psi0, nbosons; maxdim=minimum(mdim))
	end

	metadata = get(kwargs, :metadata, Dict())
	metadata["psi_ortho"] = psi_ortho
	metadata["outputlevel"] = opl
	metadata["psi0"] = psi0
	metadata["if_parton"] = if_parton
	metadata["ham"] = ham
	metadata["mdim"] = mdim
	metadata["noise"] = noise
	metadata["if_gpu"] = if_gpu
	filename = get(kwargs, :name, "mps")
	filename = check_plot_label(filename,"mps")
	
	gs_search_params = copy(metadata)
	delete!(gs_search_params,"ham")
	delete!(gs_search_params,"psi_ortho")
	delete!(gs_search_params,"psi0")
	display(gs_search_params)
	
	if isnothing(psi0) && isnothing(psi_ortho)
		sidx = siteinds("ExtendedHardcore", L; conserve_qns = conserve_qns, nflavors = nflavors)
	elseif !isnothing(psi_ortho)
		sidx = length(psi_ortho) > 1 ? siteinds(psi_ortho[1]) : siteinds(psi_ortho)
	else
		sidx = siteinds(psi0)
	end
	if isnothing(ham)
		H = MPO(hamiltonian(t1,t2,phi,U1,U2,L,nflavors; kwargs...), sidx)
	else
		H = MPO(ham, sidx)
	end
	#println("Built Hams")
	#display(matrix(combiner(dag.(sidx))*prod(H)*combiner(prime.(sidx))))
	if isnothing(psi0)
		states = make_states(L,nbosons,nflavors)
		psi0 = randomMPS(sidx, states)
	end
	if if_gpu
		H = ITensorGPU.cu(H)
		psi0 = ITensorGPU.cu(psi0)
	end
	if !isnothing(psi_ortho)
		if typeof(psi_ortho) != Vector{MPS}
			psi_ortho = [psi_ortho]
		else
			println("Using $(length(psi_ortho)) orthogonal states")
		end
		E, psi = dmrg(H, psi_ortho, psi0; maxdim = mdim, nsweeps = nsweeps, noise = noise, observer = obs, outputlevel=opl, cutoff = 1E-14)
	else
		E, psi = dmrg(H, psi0; maxdim = mdim, nsweeps = nsweeps, noise = noise, observer = obs, outputlevel=opl, cutoff = 1E-14)
	end

	metadata["observer"] = obs
	if_densmat ? densmat = density_matrix(psi) : nothing
	
	if if_save_data
		location = get(kwargs, :location, pwd())
		metadata["final_energy"] = E
		metadata["maxlinkdim"] = maxlinkdim(psi)
		metadata["final_nrg_variance"] = energy_variance(psi,H)
		data_dict::Dict{String,Any} = Dict([("mps",psi)])
		if_densmat ? data_dict["densmat"] = densmat : nothing
		write_data_jld2(filename,data_dict,location,metadata)
	end
	
	if if_nrg
		if if_densmat
			return psi, densmat, E
		else
			return psi, E
		end
	else
		if if_densmat
			return psi, densmat
		else
			return psi
		end
	end
end

function check_convergence(starting_mps,metadata,filename)
	allnrgs = metadata["observer"].energies
	if abs(allnrgs[end] - allnrgs[end-1]) < 10^-2
		println("Converged by DMRGObserver")
		return true
	else
		metadata["if_save_data"] = false
		metadata["nsweeps"] = 10
		u1 = metadata["U1"]
		u2 = metadata["U2"]
		phi = metadata["phi"]
		L = metadata["L"]
		nf = metadata["nflavors"]
		nb = metadata["nbosons"]
		println("Not Converged")
		return false
	end
end

function run_mps_new_variable(seed_wavefunc,seed_params_dict,new_params_dict,location="../cluster-data/orsay-sept23")
	newname = rewrite_filename(seed_params_dict["name"],new_params_dict)
	
	current_loc = pwd()
	cd(location)
		if_dup = check_duplicates("mps-" * newname * ".jld2") != "mps-" * newname * ".jld2"
	cd(current_loc)
	if if_dup
		println("Data file already exists, returning previously found data")
		return read_data_jld2("mps-" * newname * ".jld2",location)[1]["mps"]
	else
		new_params_dict["name"] = newname
		new_execution_dict = merge(seed_params_dict,new_params_dict)
		if "alpha" in keys(new_params_dict)
			new_execution_dict["phi"] = 2*pi*new_params_dict["alpha"]
		end
		new_wavefunc = execute_mps(new_execution_dict["U1"],new_execution_dict["U2"],new_execution_dict["phi"],new_execution_dict["L"], new_execution_dict["nflavors"],new_execution_dict["nbosons"]; dict_to_symbols(new_execution_dict)...,metadata=new_execution_dict)
		return new_wavefunc
	end
end

function varied_alpha_wavefuncs(seed_wavefunc,seed_params_dict,change=0.001)
	og_alpha = seed_params_dict["phi"] / (2*pi)
	plus_psi = run_mps_new_variable(seed_wavefunc,seed_params_dict,Dict([("phi",2*pi*(og_alpha+change)),("psi_guess",seed_wavefunc)]))
	minus_psi = run_mps_new_variable(seed_wavefunc,seed_params_dict,Dict([("phi",2*pi*(og_alpha-change)),("psi_guess",seed_wavefunc)]))
	return minus_psi,plus_psi
end

function get_mps_dims(wavefunc::MPS)
	L = length(wavefunc)
	nflavors = dim(siteind(wavefunc,1)) - 1
	return L,nflavors
end

function occupancy_variance(wavefunc::MPS; kwargs...)
	L,nflavors = get_mps_dims(wavefunc)
	fo_occ = get_occupancy(wavefunc; if_plot=false)
	occ_var = sqrt.(abs.(fo_occ.^2 .- get_occupancy(wavefunc; if_plot=false,if_squared=true)))
	if_plot = get(kwargs, :if_plot, true)
	if_plot ? mps_plot_occupancy(occ_var,L,nflavors; kwargs...) : nothing
	
	return occ_var
end

function mps_plot_occupancy(occ_mat,L,nflavors; kwargs...)
	title_string = "Occupancy, " * get(kwargs, :plot_title, "")
	fig = figure()
	#plot_surface(1:nflavors,1:L,occ_mat)
	imshow(occ_mat)
	xlabel("Virtual Dim")
	ylabel("Physical Dim")
	colorbar()
	title(title_string)
	
	if_save_fig = get(kwargs, :if_save_fig, false)
	if if_save_fig
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "occs")
		filename = check_plot_label(filename,"occs")
	end
	if_save_fig ? save_figure(filename; location=location) : nothing
	return
end

function log_sum(all_values)
	if length(all_values) < 2
		return all_values[1]
	else
		consecutive = [all_values[1],all_values[2]]
		for i in 3:length(all_values) + 1
			added_value = log_add(consecutive[1],consecutive[2])
			consecutive[1] = added_value
			consecutive[2] = i <= length(all_values) ? all_values[i] : 0.0
		end
		return consecutive[1]
	end
end

function entanglement_spectrum(psi::MPS, bond::Int)
    # Split the MPS at the bond using svd
    l, s, r = svd(psi[bond], linkind(psi, bond))

    # Form the reduced density matrix of the left part
    rho = l * s * dag(l)

    # Diagonalize the reduced density matrix to get the eigenvalues
    evals = eigvals(matrix(rho))

    # The entanglement spectrum is given by the negative logarithm of these eigenvalues
    spectrum = evals

    return spectrum
end

function entanglement_entropy(psi::MPS)
	# Get the entanglement spectrum
	spectrum = entanglement_spectrum(psi, Int(floor(length(psi)/2)))

	# Form the entanglement entropy by summing the spectrum
	entropy = -sum(spectrum .* log.(spectrum))

	return abs.(entropy)
end

function density_matrix(wavefunc::MPS; kwargs...)
	L,nflavors = get_mps_dims(wavefunc)
	densmat = zeros(L*nflavors,L*nflavors) .* im
	for s in 1:nflavors
		for sp in 1:nflavors
			#println(s,", ",sp)
			local_mat = correlation_matrix(wavefunc,"Cr$(s)","Anh$(sp)")
			densmat[L*(s-1)+1:L*s,L*(sp-1)+1:L*sp] = local_mat
			#=fig = figure()
			imshow(abs.(densmat))
			colorbar()
			title("s=$s, s'=$sp")
			=#
		end
	end
	return densmat
end

function momentum_occupation(wavefunc::MPS,p_count::Int,p_end::Real,direction="phys"; kwargs...)
	if_neg = get(kwargs, :if_neg, true)
	if_save_data = get(kwargs, :if_save_data, false)
	if_plot = p_count != 1 ? get(kwargs, :if_plot, false) : false
	p_start = get(kwargs, :p_start, 0.0)
	other_p = get(kwargs, :other_p, 0.0)
	densmat = get(kwargs, :densmat, nothing)

	if if_neg
		p_start = -p_end
	end

	phys_length,virt_length = get_mps_dims(wavefunc)

	momenta = range(p_start,stop=p_end,length=p_count)
	mom_occs = zeros(p_count) .* im
	
	for s=1:virt_length, ss=1:virt_length
		if isnothing(densmat)
			corr_val = correlation_matrix(wavefunc,"Cr$(s)","Anh$(ss)")
			corr_val += conj(transpose(corr_val))
		else
			corr_val = densmat[phys_length*(s-1)+1:phys_length*s,phys_length*(ss-1)+1:phys_length*ss]
			corr_val += conj(transpose(corr_val))
		end
		for (i,p) in enumerate(momenta)
			pvec = direction == "virt" ? [other_p,p] : [p,other_p]

			dotprod_mat = zeros(phys_length,phys_length) .* im
			for j=1:phys_length, jj=1:phys_length
				dotprod_mat[j,jj] = exp(im*pi*dot(pvec,[s-ss,j-jj]))
			end

			mom_occs[i] += sum(corr_val .* dotprod_mat)
		end
	end

	if_plot ? plot_momentum_occupation(collect(momenta),abs.(mom_occs); kwargs...) : nothing

	return momenta,mom_occs
end

function momentum_occupation(psi::MPS,p_count::Int,p_end::Real; kwargs...)
	if_neg = get(kwargs, :if_neg, true)
	p_start = get(kwargs, :p_start, 0.0)
	if_plot = get(kwargs, :if_plot, false)

	mom_occs = zeros(p_count,p_count) .* im
	momenta = Matrix{Tuple{Float64,Float64}}(undef,p_count,p_count)

	if if_neg
		p_start = -p_end
	end
	other_momenta = range(p_start*1.0,stop=p_end*1.0,length=p_count)

	for (idx,other_p) in enumerate(other_momenta)
		println(round(100*idx/length(other_momenta),digits=2),"%")
		virt_moms, local_occs = momentum_occupation(psi,p_count,p_end,"phys"; kwargs...,other_p=other_p,if_plot=false)
		momenta[idx,:] = [(virt_moms[i],other_p) for i in 1:length(virt_moms)]
		mom_occs[idx,:] = local_occs
	end

	if_plot ? plot_momentum_occupation(momenta,abs.(mom_occs); kwargs...) : nothing

	return momenta,mom_occs
end

function plot_momentum_occupation(momenta::Vector,mom_occ::Vector; kwargs...)
	title_string = "Momentum Distribution, " * get(kwargs, :plot_title, "")
	plot_label = get(kwargs, :plot_label, "")
	isempty(plot_label) ? fig = figure() : nothing
	plot(momenta,mom_occ,"-p",label=plot_label)
	if_log = get(kwargs, :if_log, false)
	if if_log
		yscale("log")
		xscale("log")
	end
	xlabel("Momenta / pi")
	ylabel("Occupation")
	title(title_string)
end

function plot_momentum_occupation(momenta::Matrix,mom_occ::Matrix; kwargs...)
	title_string = "Momentum Distribution, " * get(kwargs, :plot_title, "")
	plot_label = get(kwargs, :plot_label, "")
	fig = figure()
	imshow(mom_occ, extent=[momenta[1,1][2],momenta[end,1][2],momenta[1,1][1],momenta[1,end][1]])
	title(title_string)
	colorbar()
	xlabel("Virtual Momenta / pi")
	ylabel("Physical Momenta / pi")
end

function find_dist(p1::Tuple{Int,Int}, p2::Tuple{Int,Int}, size::Tuple{Int,Int}, periodic::Tuple{Bool,Bool}=(false, false))
    dx = abs(p1[1] - p2[1])
    dy = abs(p1[2] - p2[2])

    if periodic[1]
        dx = min(dx, size[1] - dx)
    end

    if periodic[2]
        dy = min(dy, size[2] - dy)
    end

    return sqrt(dx^2 + dy^2)
end

function physical_distance_correlation(psi::MPS; kwargs...)
	if_plot = get(kwargs, :if_plot, true)
	if_corr_lengths = get(kwargs, :if_corr_lengths, true)
	if_remapping = get(kwargs, :if_remapping, false)

	remap = if_remapping ? remapping_nnn(L) : nothing

	phys_length,virt_length = get_mps_dims(psi)
	all_corrs = [[] for i in 1:virt_length]
	dists = [i for i in 0:phys_length-1]
	for s in 1:virt_length
		if if_remapping
			linear_corrval = correlation_matrix(psi,"Cr$(s)","Anh$(s)")
			corr_val = zeros(phys_length,phys_length) .* im
			for i in 1:phys_length, j in 1:phys_length
				corr_val[i,j] = linear_corrval[remap[i],remap[j]]
			end
			corr_val ./= expect(psi,"Ns$(s)",sites=1)
		else
			corr_val = correlation_matrix(psi,"Cr$(s)","Anh$(s)")
			corr_val ./= expect(psi,"Ns$(s)",sites=1)
		end
		#corr_val += conj(transpose(corr_val))
		all_corrs[s] = abs.([mean(diag(corr_val,i)) for i in 0:phys_length-1])
		#all_corrs[s] = abs.([diag(corr_val,i)[1] for i in 0:phys_length-1])
	end

	#if_corr_lengths ? corr_lengths = correlation_length(dists,all_corrs; kwargs...) : nothing

	if if_plot
		fig = figure()
		for s in 1:virt_length
			plot(dists,abs.(all_corrs[s]),"-p",label="$s")
		end
		xlabel("Distance")
		ylabel("Correlation")
		title("Physical Distance Correlation")
		legend()
	end

	#if if_corr_lengths
	#	return dists,all_corrs,corr_lengths
	#else
		return dists,all_corrs
	#end
end

function distance_correlation(psi::MPS; kwargs...)
	if_plot = get(kwargs, :if_plot, true)
	if_periodic_phys = get(kwargs, :if_periodic_phys, true)
	if_periodic_virt = get(kwargs, :if_periodic_virt, false)
	densmat = get(kwargs, :densmat, nothing)

	phys_length,virt_length = get_mps_dims(psi)
	all_corrs = []
	dists = []
	for s=1:virt_length, ss=1:virt_length
		println(round(100*s/(virt_length),digits=2),"%")
		if isnothing(densmat)
			corr_val = correlation_matrix(psi,"Cr$(s)","Anh$(ss)")
			corr_val += conj(transpose(corr_val))
		else
			corr_val = densmat[phys_length*(s-1)+1:phys_length*s,phys_length*(ss-1)+1:phys_length*ss]
			corr_val += conj(transpose(corr_val))
		end
		for j=1:phys_length, jj=1:phys_length
			dist_btw = find_dist((s,j),(ss,jj),(virt_length,phys_length),(if_periodic_virt,if_periodic_phys))
			if dist_btw in dists
				append!(all_corrs[findfirst(x -> x == dist_btw,dists)],corr_val[j,jj])
			else
				append!(dists,[dist_btw])
				sort!(dists)
				insert!(all_corrs,findfirst(x -> x == dist_btw,dists),[corr_val[j,jj]])
			end
		end
	end

	corrs = ([mean(all_corrs[i]) for i in 1:length(all_corrs)])
	corr_errors = [std(all_corrs[i]) for i in 1:length(all_corrs)]

	if_plot ? plot_distance_correlation(dists,corrs,corr_errors; kwargs...) : nothing

	return dists,corrs,corr_errors
end

function plot_distance_correlation(dists,corrs,corr_errors; kwargs...)
	title_string = "Distance Correlation, " * get(kwargs, :plot_title, "")
	if_errors = get(kwargs, :if_errors, false)
	if_log = get(kwargs, :if_log, false)
	plot_label = get(kwargs, :plot_label, "")
	isempty(plot_label) ? fig = figure() : nothing
	if_errors ? errorbar(dists,corrs,yerr=[corr_errors,corr_errors]) : plot(dists,abs.(corrs),"-p",label=plot_label)
	if_log ? yscale("log") : nothing
	xlabel("Distance")
	ylabel("Correlation")
	title(title_string)
end

function normalize_densmat(dens_mat::Matrix,part_count::Int; kwargs...)
	if_log = get(kwargs, :if_log, false)
	L = size(dens_mat)[1]
	current_trace = if_log ? log_sum(diag(dens_mat)) : tr(dens_mat)
	if if_log
		shift_mat = Diagonal([log(part_count) - current_trace for i in 1:L])
		norm_densmat = dens_mat + shift_mat
	else
		shift_mat = ones(L,L)
		for i in 1:L
			shift_mat[i,i] = part_count/current_trace
		end
		norm_densmat = dens_mat .* shift_mat
	end
	return norm_densmat
end

function get_greenfunc(wavefunc::MPS,hopping_direction="virt"; kwargs...)
	if_backward = get(kwargs, :rev, false)

	phys_edge_length,virt_edge_length = get_mps_dims(wavefunc)
	start_point = 1
	if if_backward
		start_point = hopping_direction == "virt" ? virt_edge_length : phys_edge_length
	end
	
	all_greens = zeros(virt_edge_length,phys_edge_length).*im
	for s in 1:virt_edge_length
		all_greens[s,:] = hopping_direction == "virt" ? diag(ITensors.correlation_matrix(wavefunc,"Cr$(start_point)","Anh$(s)")) : ITensors.correlation_matrix(wavefunc,"Cr$(s)","Anh$(s)")[start_point,:]
	end
	
	start_norm = hopping_direction == "virt" ? ITensors.expect(wavefunc,"Ns$(start_point)") : [ITensors.expect(wavefunc,"Ns$(i)";sites=start_point) for i in 1:virt_edge_length]
	all_norms = zeros(virt_edge_length,phys_edge_length)
	for s in 1:virt_edge_length
		const_part = hopping_direction == "virt" ? start_norm : [start_norm[s] for i in 1:phys_edge_length]
		all_norms[s,:] = ITensors.expect(wavefunc,"Ns$(s)") .* const_part
	end
	all_greens ./= sqrt.(all_norms)
	all_greens = abs.(all_greens)
	
	if_plot = get(kwargs, :if_plot, true)
	if_plot ? plot_greenfunc(all_greens,hopping_direction; kwargs...) : nothing
	
	if_save_data = get(kwargs, :if_save_data, false)
	if if_save_data
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "$hopping_direction-dir-GF")
		filename = check_plot_label(filename,"$hopping_direction-dir-GF")
		metadata = get(kwargs, :metadata, nothing)
		data_dict = Dict([("vals",all_greens)])
	end
	if_save_data ? write_data_jld2(filename,data_dict,location,metadata) : nothing
	
	
	return all_greens
end

function plot_greenfunc(all_greens,hopping_direction; kwargs...)
	virt_edge_length,phys_edge_length = size(all_greens)
	title_string = "$hopping_direction Greens Function, " * get(kwargs, :plot_title, "")
	if_lines = get(kwargs, :if_lines, false)
	if_backward = get(kwargs, :rev, false)
	fig = figure()
	if if_lines
		if hopping_direction == "virt"
			for i in 1:phys_edge_length
				plot(1:virt_edge_length,all_greens[:,i],"-p",label="$i")
			end
			xlabel("Virtual Dimension")
		else
			for i in 1:virt_edge_length
				plot(1:phys_edge_length,all_greens[i,:],"-p",label="$i")
			end
			xlabel("Physical Dimension")
		end
	else
		imshow(all_greens)
		colorbar()
		xlabel("Physical Dim")
		ylabel("Virtual Dim")
	end
	title(title_string)
	
	if_save_fig = get(kwargs, :if_save_fig, false)
	if if_save_fig
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "$hopping_direction-dir-GF")
		filename = check_plot_label(filename,"$hopping_direction-dir-GF")
	end
	if_save_fig ? save_figure(filename; location=location) : nothing
end

function get_current_phys(wavefunc::MPS; kwargs...)
	if_exp_part = get(kwargs,:if_exp_part,true)
	alpha = get(kwargs,:alpha,0.0)
	densmat = get(kwargs, :densmat, nothing)
	if_addhc = get(kwargs, :if_addhc, false)
	L,nflavors = get_mps_dims(wavefunc)
	m0 = (nflavors)/2

	currents = [0.0*im for i in 1:nflavors]
	for i in 1:nflavors
		if isnothing(densmat)
			fullmat = correlation_matrix(wavefunc,"Cr$(i)","Anh$(i)")
			component,hc_component = diag(fullmat,-1),diag(fullmat,1)
		else
			fullmat = densmat[L*(i-1)+1:L*i,L*(i-1)+1:L*i]
			component,hc_component = diag(fullmat,-1),diag(fullmat,1)
		end
		if if_exp_part
			component .*= exp(im*pi*alpha*(i-m0)/1)
			hc_component .*= exp(-im*pi*alpha*(i-m0)/1)
		end
		currents[i] = if_addhc ? sum(component) + sum(hc_component) : sum(component) - sum(hc_component)
	end
	return sum(currents), currents
end

function get_current_synth(psi::MPS; kwargs...)
	println("Haven't made synthetic dimension current yet")
	return nothing,nothing
end

function get_current(psi::MPS,dir="phys"; kwargs...)
	if dir == "phys"
		return get_current_phys(psi; kwargs...)
	else
		return get_current_synth(psi; kwargs...)
	end
end

function get_current(all_psis::Vector{MPS},dir="phys"; kwargs...)
	if dir == "phys"
		return [get_current_phys(psi; kwargs...)[1] for psi in all_psis]
	else
		return [get_current_synth(psi; kwargs...)[1] for psi in all_psis]
	end
end
#=
function momentum_occupation(wavefunc::MPS,part_count::Int,p_count::Int64,p_end=8.0,p_start=0.0; kwargs...)
	num_sites = length(siteinds(wavefunc))
	dimension = dim(siteinds(wavefunc)[1])
	momenta = [p_start + (i-1)*(p_end - p_start)/(p_count-1) for i in 1:p_count]#[pi*i/(num_sites+1) for i in 1:p_count]
	mom_occ = zeros(p_count)*im
	for i in 1:p_count
		println(round(100*i/p_count,digits=1),"%")
		momentum = momenta[i]
		exp_vect = zeros(num_sites,num_sites) .* im
		pos_occ = zeros(num_sites,num_sites) .* im
		for j in 1:num_sites
			#exp_vect[:,j] = [-1 * cos(momentum*(l+j)) + cos(momentum*(j-l)) for l in 1:num_sites]
			exp_vect[:,j] = [exp(im*momentum*(j-l)) for l in 1:num_sites]
		end
		pos_occ = correlation_matrix(wavefunc,"FullDag","FullHat")
		mom_occ[i] = sum(exp_vect .* pos_occ) / (num_sites*(dimension-1))
	end
	
	if_plot = get(kwargs, :if_plot, true)
	if_plot ? plot_momentum_occupation(momenta,real.(mom_occ),part_count; kwargs...) : nothing
	
	if_save_data = get(kwargs, :if_save_data, false)
	if if_save_data
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "momdist")
		filename = check_plot_label(filename,"momdist")
		metadata = get(kwargs, :metadata, nothing)
		data_dict = Dict([("vals",mom_occ),("moms",momenta)])
	end
	if_save_data ? write_data_jld2(filename,data_dict,location,metadata) : nothing
	
	return momenta,mom_occ
end
=#

function bulk_density(wavefunc::MPS,bulk_width_phys=1,bulk_width_virt=1; kwargs...)
	if isnothing(wavefunc)
		occ_mat = get(kwargs, :occ_mat, nothing)
	else
		occ_mat = get_occupancy(wavefunc;if_plot=false)
	end
	size(occ_mat)[1] == size(occ_mat)[2] ? bulk_width_virt = bulk_width_phys : nothing
	num_particles = sum(occ_mat)
	bulk_occ_mat = occ_mat[1+bulk_width_phys:end-bulk_width_phys,1+bulk_width_virt:end-bulk_width_virt]
	#imshow(bulk_occ_mat)
	bulk_density = sum(bulk_occ_mat)/prod(size(bulk_occ_mat))
	return bulk_density
end

function get_densdens_corrs(wavefunc::MPS,distances=nothing; kwargs...)
	phys_edge_length,virt_edge_length = get_mps_dims(wavefunc)
	chosen_dim,other_dim = phys_edge_length,virt_edge_length
	if isnothing(distances)
		distances = [i for i in 1:chosen_dim-1]
	end
	densdens_corr = zeros(length(distances),chosen_dim)
	corr_errors = zeros(length(distances),chosen_dim)
	for i in 1:other_dim
		corr_mat = real.(ITensors.correlation_matrix(wavefunc,"Ns$(i)","Ns$(i)"))
		densdens_corr[:,i] = [mean(diag(corr_mat,j)) for j in 1:length(distances)]
		corr_errors[:,i] = [std(diag(corr_mat,j)) for j in 1:length(distances)]
	end
	
	if_plot = get(kwargs, :if_plot, true)
	get_avgs = get(kwargs, :avgs, true)
	
	if_save_data = get(kwargs, :if_save_data, false)
	if if_save_data
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "densdens")
		filename = check_plot_label(filename,"densdens")
		metadata = get(kwargs, :metadata, nothing)
	end
	
	
	if get_avgs
		avg_denscorr = [mean(densdens_corr[i,:]) for i in 1:length(distances)]
		avg_errs = [mean(corr_errors[i,:]) for i in 1:length(distances)]
		if_plot ? plot_denscorr(avg_denscorr,avg_errs,distances; kwargs...) : nothing
		data_dict = Dict([("vals",avg_denscorr),("errs",avg_errs),("dists",distances)])
		if_save_data ? write_data_jld2(filename,data_dict,location,metadata) : nothing
		return avg_denscorr,avg_errs,distances
	else
		if_plot ? plot_denscorr(densdens_corr,corr_errors,distances; kwargs...) : nothing
		data_dict = Dict([("vals",densdens_corr),("errs",corr_errors),("dists",distances)])
		if_save_data ? write_data_jld2(filename,data_dict,location,metadata) : nothing
		return densdens_corr,corr_errors,distances
	end
end

function plot_denscorr(denscorrs,corr_errors,distances; kwargs...)
	if length(size(denscorrs)) > 1
		fig = figure()
		for i in 1:size(denscorrs)[2]
			errorbar(distances,denscorrs[:,i],yerr=[corr_errors[:,i],corr_errors[:,i]],label="$i")
		end
	else
		fig = figure()
		errorbar(distances,denscorrs,yerr=[corr_errors,corr_errors],label="AVG")
	end
	
	legend()
	xlabel("Distances")
	yscale("log")
	ylabel("Corr")
	title(title_string)
	
	if_save_fig = get(kwargs, :if_save_fig, false)
	if if_save_fig
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "densdens")
		filename = check_plot_label(filename,"densdens")
	end
	if_save_fig ? save_figure(filename; location=location) : nothing
	
	return denscorrs,corr_errors
end

mutable struct NRGVarObserver <: AbstractObserver
    var_tol::Float64
    local_ham
    nrg_var::Vector{Float64}
 
    NRGVarObserver(var_tol=0.0,local_ham=10.0) = new(var_tol,local_ham,[1000.0])
 end

function ITensors.checkdone!(o::NRGVarObserver;kwargs...)
  sw = kwargs[:sweep]
  psi = kwargs[:psi]
  ham = o.local_ham
  if o.nrg_var[end] < o.var_tol && length(o.nrg_var) > 5
    #println("Stopping DMRG after sweep $sw")
    return true
  elseif length(o.nrg_var) > 10 && o.nrg_var[end] < o.var_tol*1E1 && std(o.nrg_var[end-10:end])/o.nrg_var[end] < 0.03
  	# if variance has not changed by more than 1% in the last 10 sweeps and is close to tolerance, stop
  	return true
  end
  # Otherwise, update last_energy and keep going
  append!(o.nrg_var,[energy_variance(psi,ham)])
  return false
end

function ITensors.measure!(o::NRGVarObserver; kwargs...)
    #display(kwargs)
    half_sweep = kwargs[:half_sweep]
    bond = kwargs[:bond]
    outputlevel = kwargs[:outputlevel]
  
    if bond == 1 && half_sweep == 2 && outputlevel > 0
		percent_change = length(o.nrg_var) > 10 ? round(100*std(o.nrg_var[end-10:end]/o.nrg_var[end]),digits=2) : round(100*std(o.nrg_var[2:end])/o.nrg_var[end],digits=2)
		println("The energy variance is $(round(o.nrg_var[end],digits=9)) for tolerance $(o.var_tol), AVG = $percent_change %")
    end
end

mutable struct NRGErrorObserver <: AbstractObserver
	var_tol::Float64
	local_ham
	nrgs::Vector{Float64}
    nrg_var::Vector{Float64}
 
    NRGErrorObserver(var_tol=0.0,local_ham=10.0) = new(var_tol,local_ham,[10000.0,1000.0],[0.0])
 end

function ITensors.checkdone!(o::NRGErrorObserver;kwargs...)
	sw = kwargs[:sweep]
	#psi = kwargs[:psi]
	#ham = o.local_ham
	if o.nrg_var[end] < o.var_tol && abs(o.nrgs[end] - o.nrgs[end-1]) < o.nrg_var[end]
	  return true
	end
	#o.nrgs[1] = o.nrgs[2]
	return false
end
  
function ITensors.measure!(o::NRGErrorObserver; kwargs...)
	  half_sweep = kwargs[:half_sweep]
	  bond = kwargs[:bond]
	  outputlevel = kwargs[:outputlevel]
	
	  if bond == 1 && half_sweep == 2 && outputlevel > 0
		psi = kwargs[:psi]
	    ham = o.local_ham
	    append!(o.nrg_var,[energy_variance(psi,ham)])
	    append!(o.nrgs,[real(calculate_energy(psi,ham))])
		if outputlevel > 0
			println("The energy variance is $(round(o.nrg_var[end],digits=10)) with energy change $(round(abs(o.nrgs[end] - o.nrgs[end-1]),digits=10))")
		end
	end
end

#virt_link = ind(vac[1],2)
#cr_site1 = creation_physical(siteind(vac,1),nflavors)
#M = [I matrix(cr_site1)]
#site1 = ITensor(M,inds(cr_site1)) # is (1,2) therefore need to double output dimension
#


# ╔═╡ 12309987-d529-4820-bf06-5c3407a977b3
begin
if false
	L = 16
	nflavors = 8
	nbosons = 8
	
	t1 = 1
	t2 = 0.5
	α  = 1/8
	Φ  = 2π*α
	U = 10
	U1 = 4*t1^2/U
	U2 = U1/2

	conserve_qns = true

	model_paras = (t1 = t1, t2 = t2, Φ = Φ, U1 = U1, U2 = U2, L = L,
				   nflavors = nflavors)

	sidx = siteinds("ExtendedHardcore", L; conserve_qns = conserve_qns, nflavors = nflavors)
	
	
	H = MPO(hamiltonian(; model_paras...), sidx)

	states = fill("0", L)

	jvisit = Int64[]
	for n in 1:nbosons
		jnext = rand(1:L)
		while jnext ∈ jvisit
			jnext = rand(1:L)
		end
		states[jnext] = string(rand(1:nflavors))
		push!(jvisit, jnext)
	end

	ψ₀ = randomMPS(sidx, states)

	E, ψ = dmrg(H, ψ₀; maxdim = [25,25,50,50, 100], nsweeps = 10, noise = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0])
end
end

# ╔═╡ 07470842-1097-44af-8a07-4dcff45115bc
#expect(ψ, "N") |> plot

# ╔═╡ e501cd1e-82f3-4d65-a5ce-a87743f1d6fd
let
if false
	occ_mat = zeros(L, nflavors)

	for s in 1:nflavors
		occ_mat[:, s] = expect(ψ, "Ns$(s)")
	end

	plot_surface(1:nflavors,1:L,occ_mat)
	
	#plot(occ_mat[6,:], legend = false)
end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
ITensors = "9136182c-28ba-11e9-034c-db9fb085ebd5"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
ITensors = "~0.3.35"
Plots = "~1.38.16"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.3"
manifest_format = "2.0"
project_hash = "37ffba895fd2bfb832929008b1b42cfd8e0fa5b9"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Requires", "Test"]
git-tree-sha1 = "954634616d5846d8e216df1298be2298d55280b2"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.32"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "76289dc51920fdc6e0013c872ba9551d54961c24"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.2"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables"]
git-tree-sha1 = "e28912ce94077686443433c2800104b061a827ed"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.39"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.BitIntegers]]
deps = ["Random"]
git-tree-sha1 = "abb894fb55122b4604af0d460d3018e687a60963"
uuid = "c3b6d118-76ef-56ca-8cc7-ebb389d030a1"
version = "0.3.0"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "e30f2f4e20f7f186dc36529910beaedc60cfa644"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.16.0"

[[deps.ChangesOfVariables]]
deps = ["LinearAlgebra", "Test"]
git-tree-sha1 = "f84967c4497e0e1955f9a582c232b02847c5f589"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.7"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "be6ab11021cd29f0344d5c4357b163af05a48cba"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.21.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "4e88377ae7ebeaf29a047aa1ee40826e0b708a5d"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.7.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "96d823b94ba8d187a6d8f0826e731195a74b90e9"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.2.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "738fec4d684a9a6ee9598a8bfee305b26831f28c"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.2"

[[deps.ContextVariablesX]]
deps = ["Compat", "Logging", "UUIDs"]
git-tree-sha1 = "25cc3803f1030ab855e383129dcd3dc294e322cc"
uuid = "6add18c4-b38d-439d-96f6-d6bc489c04c5"
version = "0.1.3"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "cf25ccb972fec4e4817764d01c82386ae94f77b4"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.14"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "e82c3c97b5b4ec111f3c1b55228cebc7510525a2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.25"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "e90caa41f5a86296e014e148ee061bd6c3edec96"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.9"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.ExprTools]]
git-tree-sha1 = "c1d06d129da9f55715c6c212866f5b1bddc5fa00"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.9"

[[deps.ExternalDocstrings]]
git-tree-sha1 = "1224740fc4d07c989949e1c1b508ebd49a65a5f6"
uuid = "e189563c-0753-4f5e-ad5c-be4293c83fb4"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FLoops]]
deps = ["BangBang", "Compat", "FLoopsBase", "InitialValues", "JuliaVariables", "MLStyle", "Serialization", "Setfield", "Transducers"]
git-tree-sha1 = "ffb97765602e3cbe59a0589d237bf07f245a8576"
uuid = "cc61a311-1640-44b5-9fba-1b764f453329"
version = "0.2.1"

[[deps.FLoopsBase]]
deps = ["ContextVariablesX"]
git-tree-sha1 = "656f7a6859be8673bf1f35da5670246b923964f7"
uuid = "b9860ae5-e623-471e-878b-f6a53c775ea6"
version = "0.1.1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Folds]]
deps = ["Accessors", "BangBang", "Baselet", "DefineSingletons", "Distributed", "ExternalDocstrings", "InitialValues", "MicroCollections", "Referenceables", "Requires", "Test", "ThreadedScans", "Transducers"]
git-tree-sha1 = "638109532de382a1f99b1aae1ca8b5d08515d85a"
uuid = "41a02a25-b8f0-4f67-bc48-60067656b558"
version = "0.2.8"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Functors]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "478f8c3145bb91d82c2cf20433e8c1b30df454cc"
uuid = "d9f16b24-f501-4c13-a1f2-28368ffc5196"
version = "0.4.4"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "2d6ca471a6c7b536127afccfa7564b5b39227fe0"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.5"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "8b8a2fd4536ece6e554168c21860b6820a8a83db"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "19fad9cd9ae44847fe842558a744748084a722d1"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HDF5]]
deps = ["Compat", "HDF5_jll", "Libdl", "Mmap", "Random", "Requires", "UUIDs"]
git-tree-sha1 = "c73fdc3d9da7700691848b78c61841274076932a"
uuid = "f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f"
version = "0.16.15"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "4cc2bb72df6ff40b055295fdef6d92955f9dede8"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.2+2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "7f5ef966a02a8fdf3df2ca03108a88447cb3c6f0"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.9.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.ITensors]]
deps = ["Adapt", "BitIntegers", "ChainRulesCore", "Compat", "Dictionaries", "Functors", "HDF5", "IsApprox", "KrylovKit", "LinearAlgebra", "LinearMaps", "NDTensors", "PackageCompiler", "Pkg", "Printf", "Random", "Requires", "SerializedElementArrays", "SimpleTraits", "StaticArrays", "Strided", "TimerOutputs", "TupleTools", "Zeros", "ZygoteRules"]
git-tree-sha1 = "6e234938fff8e8e3e5d2a1b87976748ebb32a82e"
uuid = "9136182c-28ba-11e9-034c-db9fb085ebd5"
version = "0.3.35"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "edd1c1ac227767c75e8518defdf6e48dbfa7c3b0"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.10"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.IsApprox]]
deps = ["Dictionaries", "LinearAlgebra"]
git-tree-sha1 = "db9c41f1ea43dce8d57b56cf36758cd135a6c2db"
uuid = "28f27b66-4bd8-47e7-9110-e2746eb8bed7"
version = "0.1.7"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "f377670cda23b6b7c1c0b3893e37451c5c1a2185"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.5"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.JuliaVariables]]
deps = ["MLStyle", "NameResolution"]
git-tree-sha1 = "49fb3cb53362ddadb4415e9b73926d6b40709e70"
uuid = "b14d175d-62b4-44ba-8fb7-3064adc8c3ec"
version = "0.2.4"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "1a5e1d9941c783b0119897d29f2eb665d876ecf3"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f689897ccbe049adb19a065c495e75f372ecd42b"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.4+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearMaps]]
deps = ["ChainRulesCore", "LinearAlgebra", "SparseArrays", "Statistics"]
git-tree-sha1 = "62f9b2762cc107667b137af621e951f52e020a0f"
uuid = "7a12625a-238d-50fd-b39a-03d52299707e"
version = "3.10.2"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "c3ce8e7420b3a6e071e0fe4745f5d4300e37b13f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.24"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.MLStyle]]
git-tree-sha1 = "bc38dff0548128765760c79eb7388a4b37fae2c8"
uuid = "d8e11817-5142-5d16-987a-aa16d5891078"
version = "0.4.17"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "629afd7d10dbc6935ec59b32daeb33bc4460a42e"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.4"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NDTensors]]
deps = ["Adapt", "Compat", "Dictionaries", "FLoops", "Folds", "Functors", "HDF5", "LinearAlgebra", "Random", "Requires", "SimpleTraits", "SplitApplyCombine", "StaticArrays", "Strided", "TimerOutputs", "TupleTools"]
git-tree-sha1 = "5ba23649531a5638b3f260a78187710f852deac1"
uuid = "23ae76d9-e61a-49c4-8f12-3f1a16adf9cf"
version = "0.1.51"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NameResolution]]
deps = ["PrettyPrint"]
git-tree-sha1 = "1a0fa0e9613f46c9b8c11eee38ebb4f590013c5e"
uuid = "71a1bf82-56d0-4bbc-8a3c-48b961074391"
version = "0.1.5"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "51901a49222b09e3743c65b8847687ae5fc78eb2"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1aa4b74f80b01c6bc2b89992b861b5f210e665b5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.21+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "d321bf2de576bf25ec4d3e4360faca399afca282"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PackageCompiler]]
deps = ["Artifacts", "LazyArtifacts", "Libdl", "Pkg", "Printf", "RelocatableFolders", "TOML", "UUIDs"]
git-tree-sha1 = "1a6a868eb755e8ea9ecd000aa6ad175def0cc85b"
uuid = "9b87118b-4619-50d2-8e1e-99f35a4d4d9d"
version = "2.1.7"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "4b2e829ee66d4218e0cef22c0a64ee37cf258c29"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.7.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "f92e1315dadf8c46561fb9396e525f7200cdc227"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.5"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "75ca67b2c6512ad2d0c767a7cfc55e75075f8bbc"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.16"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "9673d39decc5feece56ef3940e5dafba15ba0f81"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.1.2"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "7eb1686b4f04b82f96ed7a4ea5890a4f0c7a09f1"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.0"

[[deps.PrettyPrint]]
git-tree-sha1 = "632eb4abab3449ab30c5e1afaa874f0b98b586e4"
uuid = "8162dcfd-2161-5ef2-ae6c-7681170c5f98"
version = "0.2.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "0c03844e2231e12fda4d0086fd7cbe4098ee8dc5"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SerializedElementArrays]]
deps = ["Serialization"]
git-tree-sha1 = "8e73e49eaebf73486446a3c1eede403bff259826"
uuid = "d3ce8812-9567-47e9-a7b5-65a6d70a3065"
version = "0.1.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "7beb031cf8145577fbccacd94b8a8f4ce78428d3"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.0"

[[deps.SplitApplyCombine]]
deps = ["Dictionaries", "Indexing"]
git-tree-sha1 = "48f393b0231516850e39f6c756970e7ca8b77045"
uuid = "03a91e81-4c3e-53e1-a0a4-9c0c8f19dd66"
version = "1.2.2"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "832afbae2a45b4ae7e831f86965469a24d1d8a83"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.26"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "75ebe04c5bed70b91614d684259b661c9e6274a4"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.0"

[[deps.Strided]]
deps = ["LinearAlgebra", "StridedViews", "TupleTools"]
git-tree-sha1 = "b32eadf6ac726a790567fdc872b63117712e16a8"
uuid = "5e0ebb24-38b0-5f93-81fe-25c709ecae67"
version = "2.0.1"

[[deps.StridedViews]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "59cc024139c20d1ed8400c419c6fe608637d583d"
uuid = "4db3bf67-4bd7-4b4e-b153-31dc3fb37143"
version = "0.1.2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadedScans]]
deps = ["ArgCheck"]
git-tree-sha1 = "ca1ba3000289eacba571aaa4efcefb642e7a1de6"
uuid = "24d252fe-5d94-4a69-83ea-56a14333d47a"
version = "0.1.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "f548a9e9c490030e545f72074a41edfd0e5bcdd7"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.23"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "9a6ae7ed916312b41236fcef7e0af564ef934769"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.13"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "a66fb81baec325cf6ccafa243af573b031e87b00"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.77"

[[deps.TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "ba4aa36b2d5c98d6ed1f149da916b3ba46527b2b"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.14.0"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "ed8d92d9774b077c53e1da50fd81a36af3744c1c"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zeros]]
deps = ["Test"]
git-tree-sha1 = "7eb4fd47c304c078425bf57da99a56606150d7d4"
uuid = "bd1ec220-6eb4-527a-9b49-e79c3db6233b"
version = "0.3.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "49ce682769cd5de6c72dcf1b94ed7790cd08974c"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.5+0"

[[deps.ZygoteRules]]
deps = ["ChainRulesCore", "MacroTools"]
git-tree-sha1 = "977aed5d006b840e2e40c0b48984f7463109046d"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.3"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "868e669ccb12ba16eaf50cb2957ee2ff61261c56"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.29.0+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╠═07f7c480-198f-11ee-28f5-65891e30bb85
# ╠═2c9a6330-d820-4637-bbe6-9ce40dbbbb3f
# ╠═e0c40da6-0ffc-40ae-a31d-f51f6f792dac
# ╠═641ecd99-2fcc-49a7-8048-6203f8767589
# ╠═2b12ea80-dce7-4b91-b902-9a730516bf49
# ╠═2a6ecb0e-5075-425a-8cc7-fc79423b6abf
# ╠═f8382eef-be2a-4846-b378-220face8dc1b
# ╠═15f9a000-10a8-448b-b4aa-36d83b208826
# ╠═bbbb808d-1f35-414f-a429-e6ea5f885512
# ╟─615e311e-b3bd-4cd3-9272-64ba2c806787
# ╠═fa55709d-a20c-47ee-9802-8ea750588c63
# ╠═12309987-d529-4820-bf06-5c3407a977b3
# ╠═07470842-1097-44af-8a07-4dcff45115bc
# ╠═e501cd1e-82f3-4d65-a5ce-a87743f1d6fd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

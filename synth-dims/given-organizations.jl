include("long-range-ttn.jl")

function get_stripe_state(period,lat)
	phys_edge_length,virt_edge_length = sort(size(lat))
	states = fill("0",virt_edge_length*phys_edge_length)
	filled_virt_sites = [1+period*(k-1) for k in 1:Int(ceil(virt_edge_length/period))]
	filled_virt_sites[end] > virt_edge_length ? pop!(filled_virt_sites) : nothing
	for i in 1:phys_edge_length
		for j in filled_virt_sites
			states[TTNKit.linear_ind(lat,(j,i))] = "1"
		end
	end
	return states
end

function get_org_wavefunc(layers; kwargs...)
	states = get(kwargs, :states, nothing)
	period = get(kwargs, :period, 2)

	net = build_HH_net(layers)
	lat = TTNKit.physical_lattice(net)
	
	states == nothing ? states = get_stripe_state(period,lat) : nothing
	
	ttn = TTNKit.ProductTreeTensorNetwork(net,states)
	return ttn,net,lat
end

function get_func()
	eigsolve_tol = TTNKit.DEFAULT_TOL_DMRG
	eigsolve_krylovdim = TTNKit.DEFAULT_KRYLOVDIM_DMRG
	eigsolve_maxiter = TTNKit.DEFAULT_MAXITER_DMRG
	ishermitian = TTNKit.DEFAULT_ISHERMITIAN_DMRG
	eigsolve_which_eigenvalue = TTNKit.DEFAULT_WHICH_EIGENVALUE_DMRG
	func = (action,T) -> TTNKit.eigsolve(action, T, 1,
				    eigsolve_which_eigenvalue;
				    ishermitian=ishermitian,
				    tol=eigsolve_tol,
				    krylovdim=eigsolve_krylovdim,
				    maxiter=eigsolve_maxiter)
	return func
end

function get_energy(layers; kwargs...)
	pos = (1,1)
	
	alpha = get(kwargs, :alpha, 0.0)
	wavefunc = get(kwargs, :wavefunc, nothing)
	pTPO = get(kwargs, :ptpo, nothing)
	
	
	if wavefunc == nothing
		wavefunc,net,lat = get_org_wavefunc(layers; kwargs...)
	else
		net = TTNKit.network(wavefunc)
		lat = TTNKit.physical_lattice(net)
	end
	println("Got Wavefunc")
	
	if pTPO == nothing
		ham_operator = long_range_HH_ham(net,0.01,alpha; kwargs...)
		ham = TTNKit.TPO(ham_operator,lat)
		pTPO = TTNKit.ProjectedTensorProductOperator(wavefunc,ham)
	end
	println("Made pTPO")

	action = TTNKit.∂A(pTPO, pos)
	val,tn = get_func()(action,wavefunc[pos])
	nrg = val[1]
	return nrg
end


#=
layers = 5
lr = 1
period = 2
mag_off = true
states = nothing
alpha = 0.0


org_energy = get_energy(layers; scaling="flat",limit=1.0,cliff=true,if_periodic=false,no_magF=mag_off,scaling_dist=lr,states=states,period=period,alpha=alpha)
println(org_energy)
=#
















"fin"

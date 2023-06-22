include("long-range-ttn.jl")
using PyPlot

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

function get_stripe_state_mixedlaststripe(lat)
	period = 3
	states = get_stripe_state(period,lat)
	
	states[TTNKit.linear_ind(lat,(7,3))] = "0"
	states[TTNKit.linear_ind(lat,(7,4))] = "0"
	states[TTNKit.linear_ind(lat,(6,3))] = "1"
	states[TTNKit.linear_ind(lat,(6,4))] = "1"
	
	return states
end

function get_stripe_state_doublewidth(period,lat)
	phys_edge_length,virt_edge_length = sort(size(lat))
	states = fill("0",virt_edge_length*phys_edge_length)
	filled_virt_sites = [1+(period+1)*(k-1) for k in 1:Int(ceil(virt_edge_length/period))]
	filled_virt_sites[end] > virt_edge_length ? pop!(filled_virt_sites) : nothing
	for i in 1:phys_edge_length
		for j in filled_virt_sites
			states[TTNKit.linear_ind(lat,(j,i))] = "1"
			states[TTNKit.linear_ind(lat,(j+1,i))] = "1"
		end
	end
	return states
end

function make_checkers(lat)
	phys_edge_length,virt_edge_length = sort(size(lat))
	states = fill("0",virt_edge_length*phys_edge_length)
	states[1] = "1"
	for i in 1:phys_edge_length
		for j in 1:virt_edge_length
			if isodd(i) && isodd(j)
				states[TTNKit.linear_ind(lat,(j,i))] = "1"
			elseif iseven(i) && iseven(j)
				states[TTNKit.linear_ind(lat,(j,i))] = "1"
			end
		end
	end
	return states
end

function get_org_wavefunc(layers; kwargs...)
	states = get(kwargs, :states, nothing)
	period = get(kwargs, :period, 2)
	net = get(kwargs, :net, nothing)
	
	if isnothing(net)
		net = build_HH_net(layers)
	end
	lat = TTNKit.physical_lattice(net)
	
	states == nothing ? states = get_stripe_state(period,lat) : nothing
	
	ttn = TTNKit.ProductTreeTensorNetwork(net,states)
	return ttn#,net,lat
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

function make_ham(net,alpha=0.0; kwargs...)
	ham_operator = long_range_HH_ham(net,0.01,alpha; kwargs...)
	lat = TTNKit.physical_lattice(net)
	ham = TTNKit.TPO(ham_operator,lat)
	return ham
end

function get_energy(layers; kwargs...)
	pos = (1,1)
	
	alpha = get(kwargs, :alpha, 0.0)
	wavefunc = get(kwargs, :wavefunc, nothing)
	ham = get(kwargs, :ham, nothing)
	net = get(kwargs, :net, nothing)
	
	if isnothing(wavefunc)
		if isnothing(net)
			wavefunc,net,lat = get_org_wavefunc(layers; kwargs...)
		else
			wavefunc,net,lat = get_org_wavefunc(layers; net=net, kwargs...)
		end
	else
		if isnothing(net)
			net = TTNKit.network(wavefunc)
		end
		lat = TTNKit.physical_lattice(net)
	end
	println("Got Wavefunc")
	
	if isnothing(ham)
		ham = make_ham(net,alpha; kwargs...)
	end
	pTPO = TTNKit.ProjectedTensorProductOperator(wavefunc,ham)
	println("Made pTPO")

	action = TTNKit.∂A(pTPO, pos)
	val,tn = get_func()(action,wavefunc[pos])
	nrg = val[1]
	return nrg
end


if false
layers = 5
lr = 2
mag_off = true
alpha = 0.0
net = build_HH_net(layers)
lat = TTNKit.physical_lattice(net)
ham = make_ham(net,alpha; scaling_dist=lr, scaling="flat",limit=1.0,cliff=true,if_periodic=false,no_magF=mag_off)
end
#
period = 3
for states in [get_stripe_state(period,lat),get_stripe_state_mixedlaststripe(lat)]
	num_particles = sum(parse.(Int,states))
	specific = get_org_wavefunc(layers; states=states,period=period,net=net)
	org_energy = get_energy(layers; ham=ham,net=net,wavefunc=specific,alpha=alpha)
	occs = get_occupancy(specific; plot_title="Energy/part = $(org_energy/num_particles)")
end
#















"fin"

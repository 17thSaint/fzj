using TTNKit,PyPlot

#=
Need to figure out how sweeps works
=#

function get_flattened_index(b_list)
	return sum(b_list .* [2^(length(b_list) - i) for i in 1:length(b_list)]) + 1
end

function get_xy(site_number,side_length)
	y = Int(floor(site_number/side_length) + 1) 
	x = (site_number % side_length)
	if x == 0
		return x+1,y-1
	end
	return x,y
end

function get_site_number(x, y, side_length)
    site_number = (y - 1) * side_length + x
    return Int(site_number)
end

function get_ydir_greenfunc(edge_length,ttn; kwargs...)
	all_yvals = zeros(edge_length,edge_length)
	all_greens = zeros(edge_length,edge_length)
	for starting_site in 1:edge_length
		adag = "S+"
		reg_sites = (starting_site)
		prime_sites = [get_site_number(starting_site,i,edge_length) for i in 1:edge_length]
		ahat = "S-"

		all_yvals[starting_site,:] = [get_xy(prime_sites[i],edge_length)[2] for i in 1:length(prime_sites)] 
		for i in 1:length(prime_sites)
			all_greens[starting_site,i] = TTNKit.correlation(ttn,adag,ahat,reg_sites,(prime_sites[i]))
			norm_reg = TTNKit.correlation(ttn,adag,ahat,reg_sites,reg_sites)
			norm_prime = TTNKit.correlation(ttn,adag,ahat,(prime_sites[i]),(prime_sites[i]))
			all_greens[starting_site,i] /= sqrt(norm_reg * norm_prime)
		end
	end
	
	if get(kwargs, :if_plot, true)
		subplot_num = get(kwargs, :subplot_number, 111)
		if subplot_num > 150
			subplot(subplot_num)
		else
			fig = figure()
		end
		
		for i in 1:edge_length
			plot(all_yvals[i,:],all_greens[i,:],"-p",label="x=$i")
		end
	yscale("log")
	title(get(kwargs, :plot_title, "Spatial Green's Function, Edge Count = $edge_length"))
	xlabel("Y")
	ylabel("Correlation")
	legend()
	end
	
	return all_yvals,all_greens
end

function get_current_yfunc(edge_length,ttn; kwargs...)
	adag = "S+"
	ahat = "S-"
	if_periodic = get(kwargs, :if_periodic, true)
	if if_periodic
		all_yvals = zeros(edge_length,edge_length)
		all_currents = zeros(edge_length,edge_length)
	else
		all_yvals = zeros(edge_length-1,edge_length)
		all_currents = zeros(edge_length-1,edge_length)
	end
	for x in 1:size(all_currents)[1]
		x_upper = x + 1
		if x == edge_length
			if if_periodic
				x_upper = 1
			else
				continue
			end
		end
		for y in 1:size(all_currents)[2]
			all_yvals[x,y] = y
			all_currents[x,y] = TTNKit.correlation(ttn,adag,ahat,(get_site_number(x_upper,y,edge_length)),(get_site_number(x,y,edge_length))) - TTNKit.correlation(ttn,ahat,adag,(get_site_number(x_upper,y,edge_length)),(get_site_number(x,y,edge_length)))
			norm_left = TTNKit.expect(ttn,"S+ * S-",get_site_number(x,y,edge_length)) + TTNKit.expect(ttn,"S+ * S-",get_site_number(x_upper,y,edge_length))
			all_currents[x,y] /= norm_left
		end
	end
	all_currents = round.(all_currents,digits=10)
	if get(kwargs, :if_plot, true)
		subplot_num = get(kwargs, :subplot_number, 111)
		if subplot_num > 150
			subplot(subplot_num)
		else
			fig = figure()
		end
		
		for i in 1:size(all_currents)[1]
			plot(all_yvals[i,:],2 .* all_currents[i,:],"-p",label="x=$i")
		end
	title(get(kwargs, :plot_title, "Current, Edge Count = $edge_length"))
	xlabel("Y")
	ylabel("Current")
	legend()
	end
	
	return all_yvals,all_currents
end


edge_sites = 4
maxdim = 2
num_sweeps = 1
noise = 0.0

if false

square = TTNKit.BinaryNetwork((edge_sites,edge_sites), TTNKit.ITensorNode, "SpinHalf")
lat = TTNKit.physical_lattice(square)
num_sites = length(lat)
println("Finished Building Network")
#sh = TTNKit.SimpleSweepHandler(chain,num_sweeps)

states = fill("Dn",num_sites)
#states[Int(ceil(num_sites/2))] = "Dn"
ttn = TTNKit.ProductTreeTensorNetwork(square,states;orthogonalize=true)
ttn = TTNKit.increase_dim_tree_tensor_network_zeros(ttn, maxdim = maxdim)
println("Added States")

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


js, gs = -1.0, -2.0
ising = TTNKit.TransverseFieldIsing(J = js, g = gs)
tpo = TTNKit.Hamiltonian(ising,lat)
proj_tpo = TTNKit.ProjectedTensorProductOperator(ttn,tpo)
println("Finished Making Hamiltonian")

sp = TTNKit.SimpleSweepHandler(ttn,proj_tpo,func,num_sweeps,[maxdim],[noise],TTNKit.NoExpander())
TTNKit.sweep(ttn,sp);

end

rez = get_current_yfunc(edge_sites,ttn)









"fin"

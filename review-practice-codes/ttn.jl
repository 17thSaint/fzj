using TTNKit,PyPlot

#=
Need to figure out how sweeps works
=#

thetax_1,thetay_1 = 0.2, 0.12
thetax_2,thetay_2 = 0.64,0.56

function get_flattened_index(b_list)
	return sum(b_list .* [2^(length(b_list) - i) for i in 1:length(b_list)]) + 1
end

function get_xy(site_number,side_length)
	y = Int(floor(site_number/side_length) + 1) 
	x = (site_number % side_length)
	if x == 0
		return x+1,y-1
	elseif x > site_length | y > site_length
		println("ERROR: Outside Square")
		return x,y
	end
	return x,y
end

function get_site_number(x, y, side_length)
    site_number = (y - 1) * side_length + x
    if site_number > side_length^2
    	println("ERROR: Outside Square")
    	return Int(site_number)
    end
    return Int(site_number)
end

function get_ydir_greenfunc(edge_length,ttn; kwargs...)
	adag = "Adag"#"S+"
	ahat = "A"#"S-"
	all_yvals = zeros(edge_length,edge_length)
	all_greens = im.*zeros(edge_length,edge_length)
	for x in 1:edge_length
		for y in 1:edge_length
			site_left = get_site_number(x,1,edge_length)
			site_right = get_site_number(x,y,edge_length)
			
			all_greens[x,y] = TTNKit.correlation(ttn,adag,ahat,(site_left),(site_right))
			all_yvals[x,y] = y
			
			if true
			norm_reg = TTNKit.correlation(ttn,adag,ahat,(site_left),(site_left))
			norm_prime = TTNKit.correlation(ttn,adag,ahat,(site_right),(site_right))
			all_greens[x,y] /= sqrt(norm_reg * norm_prime)
			end
		end
	end
	all_greens = abs.(all_greens)
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
		title_string = "Y Spatial Green's Function, " * get(kwargs, :plot_title, "Edge Count = $edge_length")
		title(title_string)
		xlabel("Y")
		ylabel("Correlation")
		legend()
		#
		fig = figure()
		imshow(all_greens)
		xlabel("Y")
		ylabel("X")
		colorbar()
		title(title_string)
		#
	end
	
	return all_yvals,all_greens
end

function get_xdir_greenfunc(edge_length,ttn; kwargs...)
	adag = "Adag"#"S+"
	ahat = "A"#"S-"
	all_xvals = zeros(edge_length,edge_length)
	all_greens = im.*zeros(edge_length,edge_length)
	for x in 1:edge_length
		for y in 1:edge_length
			site_left = get_site_number(1,y,edge_length)
			site_right = get_site_number(x,y,edge_length)
			
			all_greens[x,y] = TTNKit.correlation(ttn,adag,ahat,(site_left),(site_right))
			all_xvals[x,y] = x
			
			if true
			norm_reg = TTNKit.correlation(ttn,adag,ahat,(site_left),(site_left))
			norm_prime = TTNKit.correlation(ttn,adag,ahat,(site_right),(site_right))
			all_greens[x,y] /= sqrt(norm_reg * norm_prime)
			end
		end
	end
	all_greens = abs.(all_greens)
	if get(kwargs, :if_plot, true)
		subplot_num = get(kwargs, :subplot_number, 111)
		if subplot_num > 150
			subplot(subplot_num)
		else
			fig = figure()
		end
		
		for i in 1:edge_length
			plot(all_xvals[:,i],all_greens[:,i],"-p",label="y=$i")
		end
		yscale("log")
		title_string = "X Spatial Green's Function, " * get(kwargs, :plot_title, "Edge Count = $edge_length")
		title(title_string)
		xlabel("Y")
		ylabel("Correlation")
		legend()
		#
		fig = figure()
		imshow(all_greens)
		xlabel("Y")
		ylabel("X")
		colorbar()
		title(title_string)
		#
	end
	
	return all_xvals,all_greens
end

function get_current_yfunc(edge_length,ttn; kwargs...)
	adag = "Adag"#"S+"
	ahat = "A"#"S-"
	norm_string = "Adag * A"#"S+ * S-"
	if_periodic = get(kwargs, :if_periodic, true)
	if if_periodic
		all_yvals = zeros(edge_length,edge_length)
		all_currents = im.*zeros(edge_length,edge_length)
	else
		all_yvals = zeros(edge_length-1,edge_length)
		all_currents = im.*zeros(edge_length-1,edge_length)
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
			norm_left = TTNKit.expect(ttn,norm_string,get_site_number(x,y,edge_length)) + TTNKit.expect(ttn,norm_string,get_site_number(x_upper,y,edge_length))
			all_currents[x,y] *= 1/norm_left
		end
	end
	all_currents = real.(round.( 2 .* all_currents,digits=10))
	if get(kwargs, :if_plot, true)
		subplot_num = get(kwargs, :subplot_number, 111)
		if subplot_num > 150 && subplot_num % 10 != 1
			subplot(subplot_num)
		else
			fig = figure()
		end
		
		for i in 1:size(all_currents)[1]
			plot(all_yvals[i,:],all_currents[i,:],"-p",label="x=$i")
		end
		title_string = "Current, " * get(kwargs, :plot_title, "Edge Count = $edge_length")
		title(title_string)
		xlabel("Y")
		ylabel("Current")
		legend()
		
	end
	return all_yvals,all_currents
end

function get_current_xfunc(edge_length,ttn; kwargs...)
	adag = "Adag"#"S+"
	ahat = "A"#"S-"
	norm_string = "Adag * A"#"S+ * S-"
	if_periodic = get(kwargs, :if_periodic, true)
	if if_periodic
		all_xvals = zeros(edge_length,edge_length)
		all_currents = im.*zeros(edge_length,edge_length)
	else
		all_xvals = zeros(edge_length-1,edge_length)
		all_currents = im.*zeros(edge_length-1,edge_length)
	end
	for y in 1:size(all_currents)[2]
		y_upper = y + 1
		if y == edge_length
			if if_periodic
				y_upper = 1
			else
				continue
			end
		end
		for x in 1:size(all_currents)[1]
			all_xvals[x,y] = x
			all_currents[x,y] = TTNKit.correlation(ttn,adag,ahat,(get_site_number(x,y_upper,edge_length)),(get_site_number(x,y,edge_length))) - TTNKit.correlation(ttn,ahat,adag,(get_site_number(x,y_upper,edge_length)),(get_site_number(x,y,edge_length)))
			norm_left = TTNKit.expect(ttn,norm_string,get_site_number(x,y,edge_length)) + TTNKit.expect(ttn,norm_string,get_site_number(x,y_upper,edge_length))
			all_currents[x,y] *= 1/norm_left
		end
	end
	all_currents = real.(round.( 2 .* all_currents,digits=10))
	if get(kwargs, :if_plot, true)
		subplot_num = get(kwargs, :subplot_number, 111)
		if subplot_num > 150 && subplot_num % 10 != 1
			subplot(subplot_num)
		else
			fig = figure()
		end
		
		for i in 1:size(all_currents)[2]
			plot(all_xvals[:,i],all_currents[:,i],"-p",label="y=$i")
		end
		title_string = "Current, " * get(kwargs, :plot_title, "Edge Count = $edge_length")
		title(title_string)
		xlabel("X")
		ylabel("Current")
		legend()
		
	end
	return all_xvals,all_currents
end


function get_chain_hofstadter(edge_length,u_strength,t_strength,phi; kwargs...)
	onsite = TTNKit.OpSum()
	interaction = TTNKit.OpSum()
	if_periodic = get(kwargs, :if_periodic, true)
	for i in 1:edge_length
		next_site = i+1
		if i == edge_length && if_periodic
			next_site = 1
		end
		
		onsite += (u_strength,"N",i,"N - Id",i)
		
		if i == edge_length && !if_periodic
			continue
		else#*exp(im*2*pi*phi*i)
			interaction += (-t_strength,"Adag",next_site,"A",i)
			interaction += (-t_strength,"Adag",i,"A",next_site)
		end
	end
	return onsite + interaction
end

function get_xy_coeffs(x,y,edge_length,t_strength,phi,thetax=thetax_1,thetay=thetay_1)
	x_coeff = -t_strength * exp(==(x,edge_length) * -im * 2 * pi * thetax)
	y_coeff = -t_strength * exp(im * 2 * pi * (phi * x - ==(y,edge_length) * thetay))
	return x_coeff, y_coeff
end

function get_inter_coeff(s1,s2,edge_length,t_strength,phi; kwargs...)
	if s1[1] == s2[1]
		thetay = get(kwargs, :thetay, thetay_2)
		return -t_strength * exp(im*2*pi*(phi*s1[1] - ==(edge_length,s1[2])*thetay))
	elseif s1[2] == s2[2]
		thetax = get(kwargs, :thetax, thetax_2)
		return -t_strength * exp(-im*2*pi* ==(edge_length,s1[1]) *thetax)
	end
end

function get_hofstadter_interacting_hamilt(edge_length,u_strength,t_strength,phi; kwargs...)
	onsite = TTNKit.OpSum()
	resulting_ham = [onsite]
	if_periodic = get(kwargs, :if_periodic, true)
	if_hopping = get(kwargs, :if_hopping, true)
	
	if if_hopping
		interaction = TTNKit.OpSum()
		lat = TTNKit.Square(edge_length,edge_length)
		for (s1,s2) in TTNKit.nearest_neighbours(lat,collect(1:TTNKit.number_of_sites(lat)); periodic=if_periodic)
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)
			
			#s1_coord = get_site_number(s1_coord_2d[1],s1_coord_2d[2],edge_length)
			#s2_coord = get_site_number(s2_coord_2d[1],s2_coord_2d[2],edge_length)
			
			coeff = get_inter_coeff(s1_coord,s2_coord,edge_length,t_strength,phi)
			interaction += (coeff,"Adag",s1_coord,"A",s2_coord)
			interaction += (conj(coeff),"Adag",s2_coord,"A",s1_coord)
		end
		append!(resulting_ham,[interaction])
	end
	for x in 1:edge_length
		for y in 1:edge_length
			#site_num = get_site_number(x,y,edge_length)
			#onsite += (u_strength,"N",site_num,"N - Id",site_num)
			onsite += (u_strength,"N",(x,y),"N - Id",(x,y))
		end
	end
	
	#=
	for x in 1:edge_length
		x_upper = x+1
		if if_periodic && x == edge_length
			x_upper = 1
		end
		for y in 1:edge_length
			y_upper = y+1
			if if_periodic && y == edge_length
				y_upper = 1
			end
			site_number = (get_site_number(x,y,edge_length))
			onsite += (u_strength,"N",site_number,"N - Id",site_number)
			
			if x == edge_length && !if_periodic
				continue
			elseif y == edge_length && !if_periodic
				continue
			else
				x_coeff,y_coeff = get_xy_coeffs(x,y,edge_length,t_strength,phi,thetax,thetay)
				
				next_x = (get_site_number(x_upper,y,edge_length))
				next_y = (get_site_number(x,y_upper,edge_length))
				interaction += (x_coeff,"Adag",next_x,"A",site_number)
				interaction += (y_coeff,"Adag",site_number,"A",next_y)
				
				interaction += (conj(x_coeff),"Adag",site_number,"A",next_x)
				interaction += (conj(y_coeff),"Adag",next_y,"A",site_number)
			end
		end
	end
	=#
	if length(resulting_ham) > 1
		return sum(resulting_ham)
	else
		return resulting_ham[1]
	end
end

function count_filled_states(states_vector)
	count = 0
	for i in 1:length(states_vector)
		if states_vector[i] == "1"
			count += 1
		end
	end
	return count
end

function fill_states(particle_count,site_count)
	states = fill("0",site_count)
	for i in 1:particle_count
		if i > site_count
			println("Too Many Particles, stopping at $site_count")
			return states
		end
		while count_filled_states(states) < i
			site = rand((1:site_count))
			if states[site] != "1"
				states[site] = "1"
			end
		end
	end
	return states
end

function build_full_harperhofstadter(edge_length,particle_count,u_strength,t_strength,filling; kwargs...)
	max_dim = get(kwargs, :max_dim, particle_count+1)

	square = TTNKit.BinaryNetwork((edge_sites,edge_sites), TTNKit.ITensorNode, "Boson")
	lat = TTNKit.physical_lattice(square)
	num_sites = length(lat)
	println("Finished Building Network")

	states = fill_states(particle_count,num_sites)
	println("Built States Vector")
	ttn = TTNKit.ProductTreeTensorNetwork(square,states;orthogonalize=true)
	ttn = TTNKit.increase_dim_tree_tensor_network_zeros(ttn, maxdim = max_dim)
	println("Added States")

	
	phi = num_particles/(filling * (edge_sites^2))
	ham_operator = get_hofstadter_interacting_hamilt(edge_length,u_strength,t_strength,phi; if_periodic=get(kwargs, :if_periodic, true),if_hopping=get(kwargs, :if_hopping, true))
	ham = TTNKit.TPO(ham_operator,lat)
	#ham = TTNKit.Hamiltonian(ham_operator,lat; mapping=TTNKit.hilbert_curve(lat))
	println("Made TPO")
	proj_tpo = TTNKit.ProjectedTensorProductOperator(ttn,ham)
	println("Finished Making Hamiltonian")
	#
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
	num_sweeps = get(kwargs, :num_sweeps, 1)
	noise = get(kwargs, :noise, 0.0)
	sp = TTNKit.SimpleSweepHandler(ttn,proj_tpo,func,num_sweeps,[max_dim],[noise],TTNKit.NoExpander())
	TTNKit.sweep(ttn,sp);
	#=
	fin_time = get(kwargs, :fin_time, 0.1)
	timestep = get(kwargs, :timestep, fin_time/25)
	sp = TTNKit.tdvp(ttn,ham; timestep=timestep, finaltime=fin_time)
	=#

	return ttn,ham,sp
end

function get_occupancy(ttn,edge_length; kwargs...)
	exp_occ = zeros(edge_length,edge_length)
	for x in 1:edge_length
		for y in 1:edge_length
			exp_occ[x,y] = TTNKit.expect(ttn,"N",(x,y))
		end
	end
	if get(kwargs, :if_plot, true)
		fig = figure()
		imshow(exp_occ)
		colorbar()
		title_string = "Occupancy, " * get(kwargs, :plot_title, "")
		title(title_string)
	end
	return exp_occ
end

function get_avg_occupancy(avg_count,edge_length,particle_count,u_strength,t_strength,filling; kwargs...)
	if_periodic = get(kwargs, :if_periodic, true)
	avg_occ_mat = zeros(edge_length,edge_length)
	for i in 1:avg_count
		ttn = build_full_harperhofstadter(edge_length,particle_count,u_strength,t_strength,filling; num_sweeps=3, if_periodic=if_periodic,max_dim=20)[1]
		avg_occ_mat += get_occupancy(ttn,edge_length; if_plot=false) ./ avg_count
	end
	avg_occ_mat = round.(avg_occ_mat,digits=3)
	fig = figure()
	imshow(avg_occ_mat)
	colorbar()
	title_string = "Avg Occupancy, " * get(kwargs, :plot_title, "")
	title(title_string)
	return avg_occ_mat
end

function get_particles_needed(edge_length;kwargs...)
	phi = get(kwargs, :phi, 0.0)
	nu = get(kwargs, :nu, 0.0)
	parts_needed = Int(ceil(phi * nu * (edge_length^2)))
	#min_dim = parts_needed + 1
	return parts_needed
end

function get_periodic_title_string(if_periodic)
	if if_periodic
		return "PBC"
	else
		return "OBC"
	end
end


#final_time = 0.1
if_periodic = true
bc_string = get_periodic_title_string(if_periodic)
edge_sites = 8
tot_sites = edge_sites^2
us = 100.0
ts = 1.0
#phi_val = 1/16
nu = 1/2
num_particles = 20#get_particles_needed(edge_sites; phi=phi_val, nu=nu)
#max_dim = num_particles + 1
println("Using $num_particles particles on $tot_sites sites")

gs_ttn, harphof_ham, hh_sp = build_full_harperhofstadter(edge_sites,num_particles,us,ts,nu; num_sweeps=3, if_periodic=if_periodic,max_dim=20)
rez = get_ydir_greenfunc(edge_sites,gs_ttn; plot_title="$bc_string")
rez2 = get_xdir_greenfunc(edge_sites,gs_ttn; plot_title="$bc_string")
rez3 = get_current_yfunc(edge_sites,gs_ttn; plot_title="$bc_string")
rez4 = get_current_xfunc(edge_sites,gs_ttn; plot_title="$bc_string")


"fin"

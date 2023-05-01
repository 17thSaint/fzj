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
	elseif x > side_length | y > side_length
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
	direct = get(kwargs, :direction, "norm")
	all_yvals = zeros(edge_length,edge_length)
	all_greens = im.*zeros(edge_length,edge_length)
	for x in 1:edge_length
		for y in 1:edge_length
			if direct == "norm"
				site_left = get_site_number(x,1,edge_length)
				site_right = get_site_number(x,y,edge_length)
			else
				site_left = get_site_number(x,y,edge_length)
				site_right = get_site_number(x,edge_length,edge_length)
			end
			
			all_greens[x,y] = TTNKit.correlation(ttn,adag,ahat,(site_left),(site_right))
			all_yvals[x,y] = y
			
			if false
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
	direct = get(kwargs, :direction, "norm")
	all_xvals = zeros(edge_length,edge_length)
	all_greens = im.*zeros(edge_length,edge_length)
	for x in 1:edge_length
		for y in 1:edge_length
			if direct == "norm"
				site_left = get_site_number(1,y,edge_length)
				site_right = get_site_number(x,y,edge_length)
			else
				site_left = get_site_number(x,y,edge_length)
				site_right = get_site_number(edge_length,y,edge_length)
			end
			
			all_greens[x,y] = TTNKit.correlation(ttn,adag,ahat,(site_left),(site_right))
			all_xvals[x,y] = x
			
			if false
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
	if get(kwargs, :no_magF, false)
		phi = 0.0
	end
	if s1[1] == s2[1]
		thetay = get(kwargs, :thetay, thetay_2)
		#=if ==(edge_length,s1[2])
			println("Using ThetaY")
		end
		=#
		return -t_strength * exp(im*2*pi*(phi*s1[1] - ==(edge_length,s1[2])*thetay))
	elseif s1[2] == s2[2]
		thetax = get(kwargs, :thetax, thetax_2)
		#=if ==(edge_length,s1[1])
			println("Using ThetaX")
		end
		=#
		return -t_strength * exp(-im*2*pi* ==(edge_length,s1[1]) *thetax)
	end
end

function get_hofstadter_interacting_hamilt(edge_length,u_strength,t_strength,phi; kwargs...)
	resulting_ham = []
	if_periodic = get(kwargs, :if_periodic, true)
	if_hopping = get(kwargs, :if_hopping, true)
	if_chem = get(kwargs, :if_chem, false)
	no_magF = get(kwargs, :no_magF, false)
	chem_strength = get(kwargs, :chem_strength, 0.0)
	#
	if if_hopping
		hopping = TTNKit.OpSum()
		lat = TTNKit.Square(edge_length,edge_length)
		for (s1,s2) in TTNKit.nearest_neighbours(lat,collect(1:TTNKit.number_of_sites(lat)); periodic=if_periodic)
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)
			
			#s1_coord = get_site_number(s1_coord_2d[1],s1_coord_2d[2],edge_length)
			#s2_coord = get_site_number(s2_coord_2d[1],s2_coord_2d[2],edge_length)
			
			coeff = get_inter_coeff(s1_coord,s2_coord,edge_length,t_strength,phi; kwargs...)
			hopping += (coeff,"Adag",s1_coord,"A",s2_coord)
			hopping += (conj(coeff),"Adag",s2_coord,"A",s1_coord)
		end
		append!(resulting_ham,[hopping])
	end
	#
	onsite = TTNKit.OpSum()
	for x in 1:edge_length
		for y in 1:edge_length
			#site_num = get_site_number(x,y,edge_length)
			#onsite += (u_strength,"N",site_num,"N - Id",site_num)
			onsite += (u_strength/2,"Adag * Adag * A * A",(x,y))
		end
	end
	append!(resulting_ham,[onsite])
	
	if if_chem
		chem = TTNKit.OpSum()
		for x in 1:edge_length
			for y in 1:edge_length
				chem -= (chem_strength,"N",(x,y))
			end
		end
		append!(resulting_ham,[chem])
	end
	
	#=
	if if_hand
	interaction = TTNKit.OpSum()
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
			
			if x == edge_length && !if_periodic
				continue
			elseif y == edge_length && !if_periodic
				continue
			else
				x_coeff = get_inter_coeff((x_upper,y),(x,y),edge_length,t_strength,phi)
				y_coeff = get_inter_coeff((x,y),(x,y_upper),edge_length,t_strength,phi)
				
				interaction += (x_coeff,"Adag",(x_upper,y),"A",(x,y))
				interaction += (y_coeff,"Adag",(x,y),"A",(x,y_upper))
				
				interaction += (conj(x_coeff),"Adag",(x,y),"A",(x_upper,y))
				interaction += (conj(y_coeff),"Adag",(x,y_upper),"A",(x,y))
			end
		end
	end
	append!(resulting_ham,[interaction])
	end
	=#
	if length(resulting_ham) > 1
		#display(sum(resulting_ham))
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
			site = rand((1:site_count))#rand([(i-1)*Int(sqrt(site_count))+Int(sqrt(site_count)) for i in 1:Int(sqrt(site_count))])
			if states[site] != "1"
				states[site] = "1"
			end
		end
	end
	return states
end

function build_full_harperhofstadter(edge_length,particle_count,u_strength,t_strength,filling; kwargs...)
	max_dim = get(kwargs, :max_dim, particle_count+1)
	num_sweeps = get(kwargs, :num_sweeps, 1)
	if_sweep = get(kwargs, :if_sweep, true)

	square = TTNKit.BinaryNetwork((edge_length,edge_length), TTNKit.ITensorNode, "Boson",conserve_number=false)
	lat = TTNKit.physical_lattice(square)
	num_sites = length(lat)
	println("Finished Building Network")

	states = fill_states(particle_count,num_sites)
	println("Built States Vector")
	ttn = TTNKit.ProductTreeTensorNetwork(square,states;orthogonalize=true)
	ttn = TTNKit.increase_dim_tree_tensor_network_zeros(ttn, maxdim = max_dim)
	println("Added States")

	#get_occupancy(ttn,edge_sites; plot_title="Starting")
	
	phi = particle_count/(filling * (edge_length^2))
	ham_operator = get_hofstadter_interacting_hamilt(edge_length,u_strength,t_strength,phi; kwargs...)
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
	
	noise = get(kwargs, :noise, 0.0)
	sp = TTNKit.SimpleSweepHandler(ttn,proj_tpo,func,num_sweeps,[max_dim],[noise],TTNKit.NoExpander())
	if if_sweep
		TTNKit.sweep(ttn,sp;outputlevel=1);
	end
	#=
	fin_time = get(kwargs, :fin_time, 0.1)
	timestep = get(kwargs, :timestep, fin_time/25)
	sp = TTNKit.tdvp(ttn,ham; timestep=timestep, finaltime=fin_time)
	=#

	return ttn,ham,sp
end

function plot_grid(edge_length)
	for i in 1:edge_length
		constant = [i for j in 1:edge_length]
		change = [k for k in 1:edge_length]
		plot(change,constant,"-pr")
		plot(constant,change,"-pr")
	end
	return
end

function plot_path(path,edge_length; kwargs...)
	path_length = length(path)
	xs = zeros(path_length)
	ys = zeros(path_length)
	for i in 1:path_length
		xs[i] = path[i][1]
		ys[i] = path[i][2]
	end
	fig = figure()
	plot_grid(edge_length)
	plot(xs,ys,"-pk";markersize=15.0,linewidth=7.0)
	plot([xs[1]],[ys[1]],"-pg";markersize=15.0)
	for i in 1:Int(ceil(length(xs)/2))
		arrow((xs[i+1]+xs[i])/2,(ys[i+1]+ys[i])/2,(xs[i+1]-xs[i])/4,(ys[i+1]-ys[i])/4,width=0.15)
	end
	xlabel("X")
	ylabel("Y")
	xlim((0.5,edge_length+0.5))
	ylim((0.5,edge_length+0.5))
	title_string = "Likely Path, " * get(kwargs, :plot_title, "")
	title(title_string)
	return xs,ys
end

function get_stop_path_redo(path,allowed=5)
	if length(path) < allowed + 1
		return path
	else
		return path[end-allowed:end]
	end
end

function find_likely_path(ttn,starting_site; kwargs...)
	net = TTNKit.network(ttn)
	lat = TTNKit.physical_lattice(net)
	edge_length = Int(sqrt(TTNKit.number_of_sites(lat)))
	path_length = get(kwargs, :path_length, Int(4*edge_length))
	path = [starting_site]
	if_periodic = get(kwargs, :periodic, false)
	all_neighbors = TTNKit.nearest_neighbours(lat,collect(1:edge_length^2); periodic=if_periodic)
	for i in 1:path_length
		if i == 1
			current_site_num = get_site_number(starting_site[1],starting_site[2],edge_length)
		end
		#current_site_num = get_site_number(current_site[1],current_site[2],edge_length)
		all_probs = []
		next_sites = []
		stopredo = get_stop_path_redo(path)
		for (s1,s2) in all_neighbors
			if s1 == current_site_num
				if !(TTNKit.coordinate(lat,s2) in stopredo)
					transition_prob = abs2.(TTNKit.correlation(ttn,"Adag","A",(s1),(s2)))
					append!(all_probs,[transition_prob])
					append!(next_sites,[s2])
				end
			elseif s2 == current_site_num
				if !(TTNKit.coordinate(lat,s1) in stopredo)
					transition_prob = abs2.(TTNKit.correlation(ttn,"Adag","A",(s2),(s1)))
					append!(all_probs,[transition_prob])
					append!(next_sites,[s1])
				end
			end
		end
		if length(all_probs) < 1
			println("Got Stuck at $i Steps")
			if get(kwargs, :if_plot, true)
				plot_path(path,edge_length; kwargs...)
			end
			return path
		end
		next_site = next_sites[findfirst(x->x==maximum(all_probs),all_probs)]
		next_coord = TTNKit.coordinate(lat,next_site)
		#println(next_coord)
		append!(path,[next_coord])
		global current_site_num = next_site
	end
	if get(kwargs, :if_plot, true)
		plot_path(path,edge_length; kwargs...)
	end
	return path
end

function get_path_direction(path,edge_length)
	corner_index_list = findall(x->x==(edge_length-1,edge_length-1),path)
	if length(corner_index_list) < 1
		corner_index_list = findall(x->x==(2,2),path)
		if length(corner_index_list) < 1
			println("Path never reached the corner")
			plot_path(path,edge_length; plot_title="Not at Corner")
			return path
		end
	end
	corner_index = corner_index_list[1]
	x_changes = [path[i+1][1] - path[i][1] for i in corner_index:corner_index+edge_length-4]
	y_changes = [path[i+1][2] - path[i][2] for i in corner_index:corner_index+edge_length-4]
	if all(x_changes.==0.0)
		return "CW"
	elseif all(y_changes.==0.0)
		return "CCW"
	else
		if corner_index <= edge_length - 3
			corner_index = corner_index_list[2]
		end
		x_changes = []
		x_changes = [path[i+1][1] - path[i][1] for i in corner_index-edge_length+3:corner_index-1]
		y_changes = []
		y_changes = [path[i+1][2] - path[i][2] for i in corner_index-edge_length+3:corner_index-1]
		if all(y_changes.==0.0)
			return "CW"
		elseif all(x_changes.==0.0)
			return "CCW"
		else
			println("Not a straight path along the edge")
			plot_path(path,edge_length; plot_title="No straight edge path")
			return x_changes,y_changes
		end
	end
end

function plot_paths_directions(paths,edge_length; kwargs...)
	cws_xs = []
	cws_ys = []
	ccws_xs = []
	ccws_ys = []
	for i in 1:length(paths)
		dir = get_path_direction(paths[i],edge_length)
		if dir == "CW"
			append!(cws_xs,[paths[i][1][1]])
			append!(cws_ys,[paths[i][1][2]])
		elseif dir == "CCW"
			append!(ccws_xs,[paths[i][1][1]])
			append!(ccws_ys,[paths[i][1][2]])
		end
	end
	fig = figure()
	plot_grid(edge_length)
	plot(cws_xs,cws_ys,"pb",markersize=15.0,label="CW")
	plot(ccws_xs,ccws_ys,"pg",markersize=15.0,label="CCW")
	xlabel("X")
	ylabel("Y")
	legend()
	xlim((0.5,edge_length+0.5))
	ylim((0.5,edge_length+0.5))
	title_string = "Path Direction" * get(kwargs, :plot_title, "")
	title(title_string)
	return cws_xs,cws_ys,ccws_xs,ccws_ys
end

function get_all_sites_paths_and_plot(ttn,edge_length; kwargs...)
	if_periodic = get(kwargs, :if_periodic, false)
	paths = []
	for i in 1:edge_length
		for j in 1:edge_length
			start = (i,j)
			println(start)
			rez = find_likely_path(ttn,start; periodic=if_periodic, if_plot=false)
			append!(paths,[rez])
		end
	end
	direction_results = plot_paths_directions(paths,edge_sites)
	return direction_results,paths
end

function get_occupancy(ttn,edge_length; kwargs...)
	#=
	exp_occ = zeros(edge_length,edge_length)
	for x in 1:edge_length
		for y in 1:edge_length
			exp_occ[x,y] = TTNKit.expect(ttn,"N",(x,y))
		end
	end
	=#
	exp_occ = abs.(TTNKit.expect(ttn,"N"))
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
		ttn = build_full_harperhofstadter(edge_length,particle_count,u_strength,t_strength,filling; kwargs...)[1]
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

function rewrite_inds(tensor,ref_tensor)
	num_layers = TTNKit.number_of_layers(TTNKit.network(tensor))
	for i in 1:num_layers
		num_tensors = length(tensor.data[i])
		for j in 1:num_tensors
			old_inds = TTNKit.inds(tensor.data[i][j])
			new_inds = TTNKit.inds(ref_tensor.data[i][j])
			if all([TTNKit.tags(old_inds[i1])==TTNKit.tags(new_inds[i1]) for i1 in 1:length(old_inds)])
				tensor.data[i][j] = TTNKit.replaceinds(tensor.data[i][j],old_inds,new_inds)
			else
				println("Index Tags Don't Match")
			end
		end
	end
	return tensor
end


function get_particles_needed(edge_length;kwargs...)
	phi = get(kwargs, :phi, 1/edge_length)
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

function localinner(ttn1::TTNKit.TreeTensorNetwork{N, T}, old_ttn2::TTNKit.TreeTensorNetwork{N, T},nexttt=false) where{N<:TTNKit.BinaryNetwork,T}

    net = TTNKit.network(ttn1)
    
    ttn2 = rewrite_inds(old_ttn2,ttn1)

    elT = promote_type(eltype(ttn1), eltype(ttn2))
    # check in case if symmetric the Top node for qn correspondence
    if !(TTNKit.sectortype(net) == Int64)
        fl1 = flux(ttn1[TTNKit.number_of_layers(net), 1])
        fl2 = flux(ttn2[TTNKit.number_of_layers(net), 2])
        if fl1 != fl2
        	println("Flux issue")
        	return zero(elT)
        end
    end

    # contruct the network starting from the first layer upwards
    #ns = number_of_sites(net)
    
    phys_lat = TTNKit.physical_lattice(net)
    if nexttt
	    println(size(res))
    end
    res = map(phys_lat) do nd
    	TTNKit.delta(TTNKit.hilbertspace(nd), TTNKit.prime(TTNKit.hilbertspace(nd)))
    end


    for ll in TTNKit.eachlayer(net)
        nt = TTNKit.number_of_tensors(net,ll)
        res_new = Vector{T}(undef, nt)
        for pp in TTNKit.eachindex(net, ll)
            childs_idx = TTNKit.getindex.(TTNKit.child_nodes(net, (ll,pp)),2)
            tn1 = ttn1[ll,pp]
            tn2 = ttn2[ll,pp]
            rpre1 = res[childs_idx[1]]
            rpre2 = res[childs_idx[2]]
            #=
            if prod(size(rpre1)) > 2^5
            	println("Stop here: ",typeof(res))
            	return res
            end
            =#
            res_new[pp] = TTNKit._dot_inner(tn1, tn2, rpre1, rpre2)
        end
        res = res_new
    end
    # better exception
    length(res) == 1 || error("Tree Tensor Contraction don't leed to a single resulting tensor.")
    res = res[1]
    
    sres = TTNKit.ITensors.scalar(res)

    return abs2(sres)
end

#final_time = 0.1
if_per = false
mag_off = false
evolve = true
chemical = true
mu = 0.5
bc_string = get_periodic_title_string(if_per)
edge_sites = 8
tot_sites = edge_sites^2
us = 1.0
ts = 0.2
#phi_val = 1/16
nu = 1/2
#for num_particles in [1,5,10]
num_particles = get_particles_needed(edge_sites; nu=nu)
mdim = 50
nswps = 2
println("Using $num_particles particles on $tot_sites sites")

#for iter in 1:1
gs_ttn, harphof_ham, hh_sp = build_full_harperhofstadter(edge_sites,num_particles,us,ts,nu; num_sweeps=nswps, if_periodic=if_per,max_dim=mdim, if_chem=chemical, chem_strength=mu, no_magF=mag_off, if_sweep=evolve)

get_all_sites_paths_and_plot(gs_ttn; if_periodic=if_per)

#get_occupancy(gs_ttn,edge_sites)
#end
#=
rez = get_ydir_greenfunc(edge_sites,gs_ttn; plot_title="N=$num_particles, $bc_string")
rez2 = get_xdir_greenfunc(edge_sites,gs_ttn; plot_title="N=$num_particles, $bc_string")
rez3 = get_ydir_greenfunc(edge_sites,gs_ttn; plot_title="N=$num_particles, $bc_string", direction="rev")
rez4 = get_xdir_greenfunc(edge_sites,gs_ttn; plot_title="N=$num_particles, $bc_string", direction="rev")
#
#
rez3 = get_current_yfunc(edge_sites,gs_ttn; plot_title="N=$num_particles, $bc_string")
rez4 = get_current_xfunc(edge_sites,gs_ttn; plot_title="N=$num_particles, $bc_string")
=#
#end
"fin"

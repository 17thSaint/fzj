using TTNKit,Statistics,NBInclude
#cd("/home/patrick/Downloads")
@nbinclude("parton-model-syms.ipynb")
include("../other-funcs/data-storage-funcs.jl")
include("../other-funcs/cluster-execution-funcs.jl")
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

#=
function get_site_number(x, y, side_length)
    site_number = (y - 1) * side_length + x
    if site_number > side_length^2
    	println("ERROR: Outside Square")
    	return Int(site_number)
    end
    return Int(site_number)
end
=#


function get_site_count(ttn)
	layers = Int(TTNKit.number_of_layers(ttn))
	return 2^layers
end

function get_lattice_dims(input_data)
	num_layers = TTNKit.number_of_layers(input_data)
	if num_layers % 2 != 0
		physical_edge_length = Int(sqrt(2^(num_layers-1)))
		virt_edge_length = physical_edge_length * 2
		return physical_edge_length,virt_edge_length
		
	else
		edge_length = Int(sqrt(2^num_layers))
		return edge_length,edge_length
	end
end

function get_ydir_greenfunc(ttn; kwargs...)
	phys_edge_length,virt_edge_length = get_lattice_dims(ttn)
	#edge_length = Int(sqrt(get_site_count(ttn)))
	adag = "Adag"#"S+"
	ahat = "A"#"S-"
	direct = get(kwargs, :direction, "norm")
	all_yvals = zeros(phys_edge_length,virt_edge_length)
	all_greens = im.*zeros(phys_edge_length,virt_edge_length)
	lat = TTNKit.physical_lattice(TTNKit.network(ttn))
	for x in 1:phys_edge_length
		for y in 1:virt_edge_length
			if direct == "norm"
				site_left = get_site_number(1,x,virt_edge_length,phys_edge_length)
				site_right = get_site_number(y,x,virt_edge_length,phys_edge_length)
			else
				site_left = get_site_number(y,x,virt_edge_length,phys_edge_length)
				site_right = get_site_number(phys_edge_length,x,virt_edge_length,phys_edge_length)
			end
			
			all_greens[x,y] = TTNKit.correlation(ttn,adag,ahat,(site_left),(site_right))
			all_yvals[x,y] = y
			
			if true
			norm_reg = TTNKit.correlation(ttn,adag,ahat,(site_left),(site_left))
			norm_prime = TTNKit.correlation(ttn,adag,ahat,(site_right),(site_right))
			if norm_reg == 0.0
				println("Normalization is zero at $x,$y")
			end
			all_greens[x,y] /= sqrt(norm_reg * norm_prime)
			end
		end
	end
	all_greens = abs.(all_greens)
	
	if_plot = get(kwargs, :if_plot, true)
	if_save_data = get(kwargs, :if_save_data, false)
	if_save_fig = get(kwargs, :if_save_fig, false)
	
	if_plot || if_save_fig ? plot_greenfunc(all_yvals,all_greens,virt_edge_length,"Y"; kwargs...) : nothing
	if_save_data ? save_greenfunc(all_yvals,all_greens,"Y"; kwargs...) : nothing
	
	return all_yvals,all_greens
end

function plot_greenfunc(all_vals,all_greens,virt_edge_length,direction; kwargs...)
	data_dict = get(kwargs, :data_dict, nothing)
	if !isnothing(data_dict)
		all_vals,all_greens = data_dict["coords"],data_dict["greens"]
	end
	title_string = "$direction Spatial Green's Function, " * get(kwargs, :plot_title, "Virt Edge Count = $virt_edge_length")
	other_direction = "X"
	if direction == "X"
		all_vals = transpose(all_vals)
		all_greens = transpose(all_greens)
		other_direction = "Y"
	end
	fig1 = figure()
	for i in 1:length(all_vals[:,1])
		plot(all_vals[i,:],all_greens[i,:],"-p",label="x=$i")
	end
	yscale("log")
	title(title_string)
	xlabel("Y")
	ylabel("Correlation")
	legend()
	#
	fig2 = figure()
	imshow(all_greens)
	xlabel("$direction")
	ylabel("$other_direction")
	colorbar()
	title(title_string)
	
	if get(kwargs, :if_save_fig, false)
		filename = get(kwargs, :name, "$direction-dir-GF")
		save_figure(filename; kwargs...)
	end
	return
end

function save_greenfunc(all_vals,all_greens,direction; kwargs...)
	filename = get(kwargs, :name, "$direction-dir-GF")
	location = get(kwargs, :location, pwd())
	greenfunc_data = Dict([("coords",all_vals),("greens",all_greens)])
	write_data_hdf5(filename,greenfunc_data,location=location; kwargs...)
	return
end

function get_xdir_greenfunc(ttn; kwargs...)
	phys_edge_length,virt_edge_length = get_lattice_dims(ttn)
	#edge_length = Int(sqrt(get_site_count(ttn)))
	adag = "Adag"#"S+"
	ahat = "A"#"S-"
	direct = get(kwargs, :direction, "norm")
	all_xvals = zeros(phys_edge_length,virt_edge_length)
	all_greens = im.*zeros(phys_edge_length,virt_edge_length)
	for x in 1:phys_edge_length
		for y in 1:virt_edge_length
			if direct == "norm"
				site_left = get_site_number(y,1,virt_edge_length,phys_edge_length)
				site_right = get_site_number(y,x,virt_edge_length,phys_edge_length)
			else
				site_left = get_site_number(y,x,virt_edge_length,phys_edge_length)
				site_right = get_site_number(y,virt_edge_length,virt_edge_length,phys_edge_length)
			end
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
	
	if_plot = get(kwargs, :if_plot, true)
	if_save_data = get(kwargs, :if_save_data, false)
	if_save_fig = get(kwargs, :if_save_fig, false)
	
	if_plot || if_save_fig ? plot_greenfunc(all_xvals,all_greens,virt_edge_length,"X"; kwargs...) : nothing
	if_save_data ? save_greenfunc(all_xvals,all_greens,"X"; kwargs...) : nothing
	
	return all_xvals,all_greens
end

function get_current_yfunc(ttn; kwargs...)
	phys_edge_length,virt_edge_length = get_lattice_dims(ttn)
	edge_length = virt_edge_length
	#edge_length = Int(sqrt(get_site_count(ttn)))
	adag = "Adag"#"S+"
	ahat = "A"#"S-"
	norm_string = "Adag * A"#"S+ * S-"
	if_periodic = get(kwargs, :if_periodic, false)
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
			all_currents[x,y] = TTNKit.correlation(ttn,adag,ahat,(get_site_number(y,x_upper,virt_edge_length,phys_edge_length)),(get_site_number(y,x,virt_edge_length,phys_edge_length))) - TTNKit.correlation(ttn,ahat,adag,(get_site_number(y,x_upper,virt_edge_length,phys_edge_length)),(get_site_number(y,x,virt_edge_length,phys_edge_length)))
			norm_left = TTNKit.expect(ttn,norm_string,get_site_number(y,x,virt_edge_length,phys_edge_length)) + TTNKit.expect(ttn,norm_string,get_site_number(y,x_upper,virt_edge_length,phys_edge_length))
			all_currents[x,y] *= 1/norm_left
		end
	end
	all_currents = real.(round.( 2 .* all_currents,digits=10))
	
	if_plot = get(kwargs, :if_plot, true)
	if_save_data = get(kwargs, :if_save_data, false)
	if_save_fig = get(kwargs, :if_save_fig, false)
	
	if_plot || if_save_fig ? plot_current(all_yvals,all_currents,edge_length,"Y"; kwargs...) : nothing
	if_save_data ? save_current(all_yvals,all_currents,"Y"; kwargs...) : nothing
	
	return all_yvals,all_currents
end

function plot_current(all_vals,all_currents,edge_length,direction; kwargs...)
	data_dict = get(kwargs, :data_dict, nothing)
	if !isnothing(data_dict)
		all_vals,all_currents = data_dict["coords"],data_dict["currents"]
	end
	if direction == "X"
		all_vals = transpose(all_vals)
		all_currents = transpose(all_currents)
	end
	
	fig = figure()
	for i in 1:size(all_currents)[1]
		plot(all_vals[i,:],all_currents[i,:],"-p",label="x=$i")
	end
	title_string = "Current, " * get(kwargs, :plot_title, "Edge Count = $edge_length")
	title(title_string)
	xlabel("$direction")
	ylabel("Current")
	legend()
	
	if get(kwargs, :if_save_fig, false)
		filename = get(kwargs, :name, "$direction-dir-current")
		save_figure(filename; kwargs...)
	end
	
	return
end

function save_current(all_vals,all_currents,direction; kwargs...)
	filename = get(kwargs, :name, "$direction-dir-current")
	location = get(kwargs, :location, pwd())
	current_data = Dict([("coords",all_vals),("currents",all_currents)])
	write_data_hdf5(filename,current_data,location=location; kwargs...)
	return
end

function get_current_xfunc(ttn; kwargs...)
	phys_edge_length,virt_edge_length = get_lattice_dims(ttn)
	edge_length = phys_edge_length
	#edge_length = Int(sqrt(get_site_count(ttn)))
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
			all_currents[x,y] = TTNKit.correlation(ttn,adag,ahat,(get_site_number(y_upper,x,virt_edge_length,phys_edge_length)),(get_site_number(y,x,virt_edge_length,phys_edge_length))) - TTNKit.correlation(ttn,ahat,adag,(get_site_number(y_upper,x,virt_edge_length,phys_edge_length)),(get_site_number(y,x,virt_edge_length,phys_edge_length)))
			norm_left = TTNKit.expect(ttn,norm_string,get_site_number(y,x,virt_edge_length,phys_edge_length)) + TTNKit.expect(ttn,norm_string,get_site_number(y_upper,x,virt_edge_length,phys_edge_length))
			all_currents[x,y] *= 1/norm_left
		end
	end
	all_currents = real.(round.( 2 .* all_currents,digits=10))
	
	if_plot = get(kwargs, :if_plot, true)
	if_save_data = get(kwargs, :if_save_data, false)
	if_save_fig = get(kwargs, :if_save_fig, false)
	
	if_plot || if_save_fig ? plot_current(all_xvals,all_currents,edge_length,"X"; kwargs...) : nothing
	if_save_data ? save_current(all_xvals,all_currents,"X"; kwargs...) : nothing
	
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

function get_inter_coeff(s1,s2,t_strength,phi,edge_length_x,edge_length_y; kwargs...) 
	#
	if get(kwargs, :no_magF, false)
		phi = 0.0
	end
	if s1[1] == s2[1]
		thetay = get(kwargs, :thetay, thetay_2)
		#=if ==(edge_length,s1[2])
			println("Using ThetaY")
		end
		=#
		return -t_strength * 1 * exp(im*2*pi*(phi*s1[1])) #- ==(edge_length_y,s1[2])*thetay))
	elseif s1[2] == s2[2]
		thetax = get(kwargs, :thetax, thetax_2)
		#=if ==(edge_length,s1[1])
			println("Using ThetaX")
		end
		=#
		return -t_strength * 1 #* exp(-im*2*pi* ==(edge_length_x,s1[1]) *thetax)
	else
		return 0.0
	end
	#=
	no_magF = get(kwargs, :no_magF, false)
	if no_magF
		phi = 0.0
	end
	if s1[1] == s2[1]
		return -t_strength
	else
		return -t_strength * exp(im * 2 * pi * phi * s1[1])
	end
	=#
end

function v_central(sz, Vj)
    V = zeros(sz)
    V[sz[1]÷2, sz[2]÷2] = Vj
    V[sz[1]÷2, sz[2]÷2 + 1] = Vj
    V[sz[1]÷2 + 1, sz[2]÷2] = Vj
    V[sz[1]÷2 + 1, sz[2]÷2 + 1] = Vj
    return V
end
#no_v(sz) = zeros(sz)

function get_hofstadter_interacting_hamilt(net,t_strength,phi; kwargs...)
	resulting_ham = []
	if_periodic = get(kwargs, :if_periodic, true)
	if_hopping = get(kwargs, :if_hopping, true)
	if_chem = get(kwargs, :if_chem, false)
	if_onsite = get(kwargs, :if_onsite, true)
	if_pinning_pot = get(kwargs, :if_pinning_pot, false)
	vpinning = get(kwargs, :vpinning, 2.5)
	no_magF = get(kwargs, :no_magF, false)
	chem_strength = get(kwargs, :chem_strength, 0.0)
	u_strength = get(kwargs, :u_strength, 1.0)
	
	lat = TTNKit.physical_lattice(net)
	edge_length_x,edge_length_y = size(lat)
	#
	if if_hopping
		hopping = TTNKit.OpSum()
		for (s1,s2) in TTNKit.nearest_neighbours(lat,collect(1:TTNKit.number_of_sites(lat)); periodic=if_periodic)
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)
			
			#s1_coord = get_site_number(s1_coord_2d[1],s1_coord_2d[2],edge_length)
			#s2_coord = get_site_number(s2_coord_2d[1],s2_coord_2d[2],edge_length)
			
			coeff = get_inter_coeff(s1_coord,s2_coord,t_strength,phi,edge_length_x,edge_length_y; kwargs...)
			hopping += (coeff,"Adag",s1_coord,"A",s2_coord)
			hopping += (conj(coeff),"Adag",s2_coord,"A",s1_coord)
		end
		append!(resulting_ham,[hopping])
	end
	if if_onsite
		onsite = TTNKit.OpSum()
		for i in TTNKit.eachindex(lat)
			onsite += (u_strength/2,"Adag * Adag * A * A",TTNKit.coordinate(lat,i))
		end
		append!(resulting_ham,[onsite])
	end
	#=
	for x in 1:edge_length
		for y in 1:edge_length
			#site_num = get_site_number(x,y,edge_length)
			#onsite += (u_strength,"N",site_num,"N - Id",site_num)
			onsite += (u_strength/2,"Adag * Adag * A * A",(x,y))
		end
	end
	=#
	
	if if_chem && chem_strength != 0.0
		chem = TTNKit.OpSum()
		for i in TTNKit.eachindex(lat)
			chem -= (chem_strength,"N",TTNKit.coordinate(lat,i))
		end
		#=
		for x in 1:edge_length
			for y in 1:edge_length
				chem -= (chem_strength,"N",(x,y))
			end
		end
		=#
		append!(resulting_ham,[chem])
	end
	
	if if_pinning_pot
		Vj = v_central(size(lat), vpinning)
		pinning_pot = TTNKit.OpSum()
		for p in TTNKit.coordinates(lat)
	    		pinning_pot += (Vj[p[1],p[2]], "N", p)
		end
		append!(resulting_ham,[pinning_pot])
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

function get_position_dims(ttn)
	for pos in TTNKit.NodeIterator(TTNKit.network(ttn))
		println("Position $pos: ",TTNKit.ITensors.dims(ttn[pos]))
	end
	return
end

function count_filled_states(states_vector)
	count = 0
	for i in 1:length(states_vector)
		if states_vector[i] != "0"
			count += parse(Int64,states_vector[i])
		end
	end
	return count
end

function localsweep(psi0::TTNKit.TreeTensorNetwork, sp::TTNKit.AbstractSweepHandler; kwargs...)
    
    obs = get(kwargs, :observer, TTNKit.NoObserver())

    outputlevel = get(kwargs, :outputlevel, 1)

    # now start with the sweeping protocol
    TTNKit.initialize!(sp)
    #sp = SimpleSweepProtocol(net, n_sweeps)
    for sw in TTNKit.sweeps(sp)
        if outputlevel ≥ 2 
            println("Start sweep number $(sw)")
            flush(stdout)
        end
        t_p = time()
        for pos in sp
            TTNKit.update!(sp, pos)
            println("Position $pos")
            get_position_dims(sp.ttn)
            TTNKit.measure!(
                obs;
                sweep_handler=sp,
                pos=pos,
                outputlevel=outputlevel
            )
        end
        t_f = time()
        if outputlevel ≥ 1
            print("Finsihed sweep $sw. ")
            #@printf("Needed Time %.3fs\n", t_f - t_p)
            # additional info string provided by the sweephandler
            TTNKit.info_string(sp, outputlevel)
            #@printf("\n")
            flush(stdout)
        end
    end
    return sp
end

function fill_states(particle_count,site_count,max_occupation)
	states = fill("0",site_count)
	for i in 1:particle_count
		#
		if i > site_count * max_occupation
			println("Too Many Particles, stopping at $site_count each with $max_occupation particles")
			return states
		end
		#
		while count_filled_states(states) < i
			site = rand((1:site_count))#rand([(i-1)*Int(sqrt(site_count))+Int(sqrt(site_count)) for i in 1:Int(sqrt(site_count))])
			#if states[site] != "1"
			#	states[site] = "1"
			#end
			#
			current_occupation = parse(Int64,states[site])
			if current_occupation < max_occupation
				states[site] = string(current_occupation + 1)
			end
		end
	end
	return states
end

function initialize_ttn(ttn,maxdim,particle_count)
	phys_edge_length,virt_edge_length = get_lattice_dims(ttn)
	site_count = TTNKit.number_of_sites(TTNKit.network(ttn))
	wf_coefs = create_wavefunction(Float64,size(TTNKit.physical_lattice(TTNKit.network(ttn))))
	for i in 1:particle_count
		ttn = patron_application!(ttn,wf_coefs,"Adag";maxdim=maxdim)
	end
	return ttn
end

function check_if_frozen(ttn)
	occs = get_occupancy(ttn; if_plot=false)
	if any(isapprox.(occs,0.0,atol=10^-10))
		return true,"frozen"
	elseif any(round.(get_ydir_greenfunc(ttn;if_plot=false)[2],digits=3).==0.0)#sum(occs.==1.0) > length(occs)/3
		return true,"variables"
	else
		return false,"none"
	end
end

function do_sweep(ttn,ham,sweep_type,particle_count; kwargs...)

	opl = get(kwargs, :output_level, 0)
	max_dim = get(kwargs, :max_dim, particle_count+1)
	num_sweeps = get(kwargs, :num_sweeps, 1)
	noise = get(kwargs, :noise, 0.0)
	expander = get(kwargs, :expander, TTNKit.NoExpander())
	
	#println("PreSweep Link Dim = ",TTNKit.maxlinkdim(ttn))
	#get_position_dims(ttn)
	if sweep_type == "dmrg"
		sp = TTNKit.dmrg(ttn,ham; expander=expander, number_of_sweeps=num_sweeps, maxdims=max_dim, noise=noise, output_level=opl)
	elseif sweep_type == "simple"
		proj_tpo = TTNKit.ProjectedTensorProductOperator(ttn,ham)
		#println("Finished Making Hamiltonian")
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
		
		#ttnc = TTNKit.copy(ttn)
		#ttnc = TTNKit.move_ortho!(ttnc,(TTNKit.number_of_layers(TTNKit.network(ttnc)),1))
		sp = TTNKit.SimpleSweepHandler(ttnc,proj_tpo,func,num_sweeps,[max_dim],[noise],expander)
		#println("Sweep Built Link Dim = ",TTNKit.maxlinkdim(sp.ttn))
		TTNKit.sweep(ttnc,sp;outputlevel=opl);
		return ttnc,ham,sp
		#println("PostSweep TTN Link Dim = ",TTNKit.maxlinkdim(ttn))
		#println("PostSweep SP-TTN Link Dim = ",TTNKit.maxlinkdim(sp.ttn))
	end
	
	return ttn,ham,sp
end

function warming(ttn,ham,sp,particle_count,warming_limit; kwargs...)
	
	max_dim = get(kwargs, :max_dim, particle_count+1)
	num_sweeps = 3#get(kwargs, :num_sweeps, 1)
	noise = get(kwargs, :noise, 0.0)
	expander = get(kwargs, :expander, TTNKit.NoExpander())
	sweep_type = get(kwargs, :sweep_type, "dmrg")

	warming_count = 1
	frozen = true
	global old_data = [ttn,ham,sp]
	while frozen && warming_count < warming_limit
		new_maxdim = Int((warming_count+10)*max_dim/10)
		reexpanded_ttn = TTNKit.adjust_tree_tensor_dimensions(old_data[1],new_maxdim)
		new_ttn, new_ham, new_sp = do_sweep(reexpanded_ttn,ham,sweep_type,particle_count; kwargs...)
		println("Max Dim = ",TTNKit.maxlinkdim(new_sp.ttn),", Expected = $new_maxdim")
		if_frozen,why = check_if_frozen(new_sp.ttn)
		if if_frozen
			get_occupancy(new_sp.ttn; plot_title="Attempt $warming_count")
			warming_count += 1
			global old_data = [new_sp.ttn,new_ham,new_sp]
		else
			println("Stable Result Found in $warming_count Attempts")
			return new_sp.ttn,new_ham,new_sp
		end
	end
	println("Hit warming limit, still frozen")
	return old_data
end

function get_excess_particles(part_count,site_count)
	if part_count > 0.5*site_count
		return Int(abs(part_count - site_count))
	else
		return part_count
	end
end

function throwout_therm_time(times)
	if times[1] > 5 * mean(times[2:end])
		return times[2:end]
	else
		return times
	end
end

function build_full_harperhofstadter(num_layers,particle_count,t_strength,filling; kwargs...)
	num_sites = 2^num_layers
	max_dim = get(kwargs, :max_dim, particle_count+1)
	num_sweeps = get(kwargs, :num_sweeps, 3)
	sweep_iter = get(kwargs, :sweep_iter, 1)
	if_sweep = get(kwargs, :if_sweep, true)
	sweep_type = get(kwargs, :sweep_type, "simple")
	noise = get(kwargs, :noise, 0.0)
	expander = get(kwargs, :expander, TTNKit.NoExpander())
	max_occ = get(kwargs, :max_occ, Int(round(particle_count/(num_sites))+1) )
	u_strength = get(kwargs, :u_strength, 100.0)
	warming_limit = get(kwargs, :warming_limit, 100)
	conserve_qns = get(kwargs, :syms, true)
	#excess_particles = get_excess_particles(particle_count,num_sites)
	#phi = get(kwargs, :phi, excess_particles/(filling * (num_sites)))
	phi = get(kwargs, :phi, particle_count/(filling * (num_sites)))
	ham_operator = get(kwargs, :ham_op, nothing)
	net = get(kwargs, :ttn_net, nothing)

	if isnothing(net)
		net = TTNKit.BinaryRectangularNetwork(num_layers, TTNKit.ITensorNode, "Boson";conserve_qns=conserve_qns,dim=max_occ+1)
	end
	lat = TTNKit.physical_lattice(net)

	println("Finished Building Network")

	states = fill("0", num_sites)#fill_states(particle_count-1,num_sites,max_occ)
	old_ttn = TTNKit.ProductTreeTensorNetwork(net,states)
	#ttn = TTNKit.increase_dim_tree_tensor_network_zeros(ttn, maxdim = max_dim)
	#ttn = TTNKit.adjust_tree_tensor_dimensions(old_ttn,max_dim)
	#println("Starting Link Dim = ",TTNKit.maxlinkdim(old_ttn))
	ttn = initialize_ttn(old_ttn,max_dim,particle_count)
	#println("Adjusted Link Dim = ",TTNKit.maxlinkdim(ttn))
	println("Added States")
	
	#get_occupancy(ttn,edge_sites; plot_title="Starting")
	
	if isnothing(ham_operator)
		ham_operator = get_hofstadter_interacting_hamilt(net,t_strength,phi; kwargs...)
	end
	
	ham = TTNKit.TPO(ham_operator,lat)
	println("Built Hamiltonian")
	sp = 0
	times = []
	if if_sweep
		for i in 1:sweep_iter
			time_start = time()
			new_ttn, new_ham, new_sp = do_sweep(ttn,ham,sweep_type,particle_count; kwargs...)
			time_end = time()
			append!(times,[time_end - time_start])
			#return sp.ttn, ham, sp
			if_frozen,why = check_if_frozen(new_sp.ttn)
			if !if_frozen
				#get_position_dims(sp.ttn)
				#return new_sp.ttn, new_ham, new_sp
				ttn,ham,sp = new_ttn,new_ham,new_sp
			else
				if why == "frozen"
					println("Frozen on First Attempt, Starting Warming")
				elseif why == "variables"
					println("Bad Variables on First Attempt, Starting Reset")	
				end
				warmed_results = warming(new_ttn,new_ham,new_sp,particle_count,warming_limit;kwargs...)
				#return warmed_results
				ttn,ham,sp = warmed_results
			end
		end
		return ttn, ham, sp#, throwout_therm_time(times)
	end

	return ttn,ham,"no sweep"
end

function plot_grid(virt_edge_length,phys_edge_length)
	for i in 1:virt_edge_length
		constant = [i for j in 1:phys_edge_length]
		change = [k for k in 1:virt_edge_length]
		plot(change,constant,"-pr")
		plot(constant,change,"-pr")
	end
	return
end

function plot_path(path,phys_edge_length,virt_edge_length; kwargs...)
	path_length = length(path)
	xs = zeros(path_length)
	ys = zeros(path_length)
	for i in 1:path_length
		xs[i] = path[i][1]
		ys[i] = path[i][2]
	end
	fig = figure()
	plot_grid(virt_edge_length,phys_edge_length)
	plot(xs,ys,"-pk";markersize=15.0,linewidth=7.0)
	plot([xs[1]],[ys[1]],"-pg";markersize=15.0)
	for i in 1:Int(ceil(length(xs)/2))
		arrow((xs[i+1]+xs[i])/2,(ys[i+1]+ys[i])/2,(xs[i+1]-xs[i])/4,(ys[i+1]-ys[i])/4,width=0.15)
	end
	xlabel("X")
	ylabel("Y")
	xlim((0.5,phys_edge_length+0.5))
	ylim((0.5,virt_edge_length+0.5))
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

function normalize_transition_probs(given_probs)
	og_sum = sum(given_probs)
	return (1/og_sum) .* given_probs
end

function get_rand_next_site(probabilities)
	normed_probs = normalize_transition_probs(probabilities)
	x = rand()
	s = 0.0
	for i in 1:length(normed_probs)
        	s += normed_probs[i]
        	if s >= x
            		return i
        	end
    	end
    	println("Error: Probabilities didn't work")
    	return
end

function find_path(ttn,starting_site; kwargs...)
	net = TTNKit.network(ttn)
	lat = TTNKit.physical_lattice(net)
	phys_edge_length,virt_edge_length = get_lattice_dims(ttn)
	#edge_length = Int(sqrt(TTNKit.number_of_sites(lat)))
	path = [starting_site]
	
	if_periodic = get(kwargs, :periodic, false)
	likely_path = get(kwargs, :likely_path, true)
	rand_path = get(kwargs, :rand_path, !likely_path)
	path_length = get(kwargs, :path_length, Int(4*virt_edge_length))	
	
	
	all_neighbors = TTNKit.nearest_neighbours(lat,collect(1:phys_edge_length*virt_edge_length); periodic=if_periodic)
	for i in 1:path_length
		if i == 1
			current_site_num = get_site_number(starting_site[2],starting_site[1],virt_edge_length,phys_edge_length)
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
		if likely_path
			next_site = next_sites[findfirst(x->x==maximum(all_probs),all_probs)]
		elseif rand_path
			next_site = next_sites[get_rand_next_site(all_probs)]
		end
		next_coord = TTNKit.coordinate(lat,next_site)
		#println(next_coord)
		append!(path,[next_coord])
		global current_site_num = next_site
	end
	if get(kwargs, :if_plot, true)
		plot_path(path,phys_edge_length,virt_edge_length; kwargs...)
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
			println(corner_index_list,", ",path)
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

function plot_paths_directions(cws_xs,cws_ys,ccws_xs,ccws_ys,edge_length; kwargs...)
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
	return
end

function make_paths_directions(paths,edge_length; kwargs...)
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
	if get(kwargs, :if_plot, true)
		plot_paths_directions(cws_xs,cws_ys,ccws_xs,ccws_ys,edge_length; kwargs...)
	end
	return cws_xs,cws_ys,ccws_xs,ccws_ys
end

function get_all_sites_paths_and_plot(ttn; kwargs...)
	edge_length = Int(sqrt(get_site_count(ttn)))
	if_periodic = get(kwargs, :if_periodic, false)
	likely_path = get(kwargs, :likely_path, false)
	paths = []
	for i in 1:edge_length
		for j in 1:edge_length
			start = (i,j)
			println(start)
			rez = find_path(ttn,start; likely_path=likely_path, periodic=if_periodic, if_plot=false)
			append!(paths,[rez])
		end
	end
	return paths
	direction_results = make_paths_directions(paths,edge_length; kwargs...)
	return direction_results,paths
end

function get_2part_corr(ttn,particle_count; kwargs...)
	corr = sum(abs.(TTNKit.expect(ttn,"N * N - N")))/particle_count
	return corr
end

function get_occupancy(ttn; kwargs...)
	#=
	exp_occ = zeros(edge_length,edge_length)
	for x in 1:edge_length
		for y in 1:edge_length
			exp_occ[x,y] = TTNKit.expect(ttn,"N",(x,y))
		end
	end
	=#
	exp_occ = abs.(TTNKit.expect(ttn,"N"))
	
	if_save_data = get(kwargs, :if_save_data, false)
	if_save_fig = get(kwargs, :if_save_fig, false)
	if_plot = get(kwargs, :if_plot, true)
	
	
	if_plot	|| if_save_fig ? plot_occupancy(exp_occ; kwargs...) : nothing
	if_save_data ? save_occupancy(exp_occ; kwargs...) : nothing
		
	return exp_occ
end

function save_occupancy(exp_occ; kwargs...)
	location = get(kwargs, :location, pwd())
	filename = get(kwargs, :name, "occs")
	metadata = get(kwargs, :metadata, nothing)
	occs_data_dict = Dict([("vals",exp_occ)])
	write_data_hdf5(filename,occs_data_dict,location,metadata)
	return
end

function plot_occupancy(exp_occ; kwargs...)
	data_dict = get(kwargs, :data_dict, nothing)
	if !isnothing(data_dict)
		exp_occ = data_dict["vals"]
	end
	fig = figure()
	imshow(exp_occ)
	colorbar()
	plot_title = get(kwargs, :plot_title, "")
	title_string = "Occupancy, " * plot_title
	title(title_string)
	
	if get(kwargs, :if_save_fig, false)
		location = get(kwargs, :location, pwd())
		fig_name = get(kwargs, :name, "occs")
		save_figure(fig_name; kwargs...)
	end
	return
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
			#
			for k in new_inds
				matching = findfirst(x->x==k,old_inds)
				if matching != nothing
					tensor.data[i][j] = TTNKit.replaceinds(tensor.data[i][j],old_inds[matching],k)
				end
			end
			#=
			do_inds_match = [TTNKit.tags(old_inds[i1])==TTNKit.tags(new_inds[i1]) for i1 in 1:length(old_inds)]
			if all(do_inds_match)
				tensor.data[i][j] = TTNKit.replaceinds(tensor.data[i][j],old_inds,new_inds)
			else
				nonmatching_sites = findall(x->x==false,do_inds_match)
				nonmatching_old = [TTNKit.tags(old_inds[i1]) for i1 in nonmatching_sites]
				nonmatching_new = [TTNKit.tags(new_inds[i1]) for i1 in nonmatching_sites]
				println("Index Tags Don't Match")
				display(nonmatching_old)
				display(nonmatching_new)
			end
			=#
		end
	end
	return tensor
end


function get_particles_needed(num_layers;kwargs...)
	if num_layers % 2 != 0.0
		edge_length = Int(sqrt(2^(num_layers-1)))
	else
		edge_length = Int(sqrt(2^num_layers))
	end
	phi = get(kwargs, :phi, 1/edge_length)
	nu = get(kwargs, :nu, 0.0)
	parts_needed = Int(ceil(phi * nu * (2^num_layers)))
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
    
    ttn2 = rewrite_inds(copy(old_ttn2),ttn1)

    elT = promote_type(eltype(ttn1), eltype(ttn2))
    # check in case if symmetric the Top node for qn correspondence
    if !(TTNKit.sectortype(net) == Int64)
        fl1 = TTNKit.flux(ttn1[TTNKit.number_of_layers(net), 1])
        fl2 = TTNKit.flux(ttn2[TTNKit.number_of_layers(net), 1])
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
            #
            if prod(size(rpre1)) > 2^5
            	println("Stop here: ",typeof(res))
            	return res
            end
            #
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

function get_phase_diag_MISF(ts_start,ts_end,ts_count,chem_start,chem_end,chem_count,edge_length,particles_count,us,nu; kwargs...)
	tss = [ts_start + (i-1)*(ts_end-ts_start)/ts_count for i in 1:ts_count+1]
	mus = [chem_start + (i-1)*(chem_end-chem_start)/chem_count for i in 1:chem_count+1]
	mis_m = []
	mis_t = []
	sfs_m = []
	sfs_t = []
	for i in 1:length(tss)
		ts = ts[i]
		for j in 1:length(mus)
			chem = mus[i]
			ttn, harphof_ham, hh_sp = build_full_harperhofstadter(edge_length,particles_count,us,ts,nu; chem_strength=chem, kwargs...)
			rez = get_ydir_greenfunc(edge_length,ttn; if_plot=false)
			if all(rez[2][:,2:end] .< 10^-20)
				append!(mis_m,[mu])
				append!(mis_t,[ts])
			else
				append!(sfs_m,[mu])
				append!(sfs_t,[ts])
			end
		end
	end
	fig = figure()
	scatter(mis_t,mis_m,c="b",label="MI")
	scatter(sfs_t,sfs_m,c="g",label="SF")
	xlabel("hopping strength")
	ylabel("chemical strength")
	legend()
	return mis_m,mis_t,sfs_m,sfs_t
end

function get_mu_transitions(ts,num_layers,n=1; kwargs...)
	if get(kwargs, :max_occ, 4) < n + 1
		println("Max Occupation not larger enough for $n particle filling")
		return
	end
	if_mm = get(kwargs, :mm, true)
	if_mp = get(kwargs, :mp, true)
	if if_mm && if_mp
		results = [0.0,0.0]
	else
		results = [0.0]
	end
	
	full_ttn,full_ham,hh_sp_full = build_full_harperhofstadter(num_layers,n*(2^num_layers),ts,1/2; if_chem=false, no_magF=true, if_sweep=true, kwargs...)
	energy_full = hh_sp_full.current_energy
	
	#Threads.@threads for i in 1:2
	#	if i == 1
			if if_mm
				minus_ttn,minus_ham,hh_sp_minus = build_full_harperhofstadter(num_layers,n*(2^num_layers)-1,ts,1/2; if_chem=false, no_magF=true, if_sweep=true, kwargs...)
				energy_minus = hh_sp_minus.current_energy
				mu_minus = energy_full - energy_minus
				results[1] = mu_minus
			end
	#	else
			if if_mp
				plus_ttn,plus_ham,hh_sp_plus = build_full_harperhofstadter(num_layers,n*(2^num_layers)+1,ts,1/2; if_chem=false, no_magF=true, if_sweep=true, kwargs...)
				energy_plus = hh_sp_plus.current_energy
				mu_plus = energy_plus - energy_full
				if length(results) > 1
					results[2] = mu_plus
				else
					results[1] = mu_plus
				end
			end	
	#	end
	#end	
	return results
end

function get_mu_trans_range_t(ts_start,ts_end,ts_count,num_layers,n=1; kwargs...)
	tss = [ts_start + i*(ts_end-ts_start)/ts_count for i in 0:ts_count]
	if_mm = get(kwargs, :mm, true)
	if_mp = get(kwargs, :mp, true)
	mms = [0.0 for i in 1:length(tss)]
	mps = [0.0 for i in 1:length(tss)]
	#Threads.@threads for i in 1:length(tss)
	for i in 1:length(tss)
		println(round(100*i/length(tss),digits=2),"%")
		ts = tss[i]
		results = get_mu_transitions(ts,num_layers,n; mm=true, mp=true, kwargs...)
		if if_mm
			mms[i] = results[1]
			if if_mp
				mps[i] = results[2]
			end
		elseif if_mp
			mps[i] = results[1]
		end
	end
	if get(kwargs, :if_plot, true)
		fig = figure()
		if if_mm
			plot(tss,mms,"-p",label="Minus")
		end
		if if_mp
			plot(tss,mps,"-p",label="Plus")
		end
		xlabel("Hopping Strength, t")
		ylabel("Chemical Potential, mu")
		title_string = "Mott Transition, " * get(kwargs, :plot_title, "")
		title(title_string)
	end
	if length(mms) < 1
		return mps,tss
	elseif length(mps) < 1
		return mms,tss
	else
		return mms,mps,tss
	end
end

function get_mu_trans_range_N(ts,edge_start,edge_count,n=1; kwargs...)
	edges = [edge_start*(2^i) for i in 0:edge_count-1]
	if_mm = get(kwargs, :mm, true)
	if_mp = get(kwargs, :mp, true)
	mms = []
	mps = []
	for j in 1:length(edges)
		results = get_mu_transitions(ts,edges[j],n; mm=true, mp=true, kwargs...)
		if if_mm
			append!(mms,[results[1]])
			if if_mp
				append!(mps,[results[2]])
			end
		elseif if_mp
			append!(mps,[results[1]])
		end
	end
	if get(kwargs, :if_plot, true)
		fig = figure()
		if if_mm
			plot(1 ./ edges,mms,"-p",label="Minus")
		end
		if if_mp
			plot(1 ./ edges,mps,"-p",label="Plus")
		end
		xlabel("1/Nsites")
		ylabel("Chemical Potential, mu")
		title_string = "Mu Transitions, " * get(kwargs, :plot_title, "")
		title(title_string)
	end
	if length(mms) < 1
		return mps,tss
	elseif length(mps) < 1
		return mms,tss
	else
		return mms,mps,tss
	end
end

function get_linfit(x,y)
	model(t,p) = p[1] .+ p[2] .* t
	p0 = [0.5,1.0]
	fit = curve_fit(model,x,y,p0)
	#fig = figure()
	#plot(x,model(x,fit.param))
	#plot(x,y)
	return fit.param[1],LsqFit.standard_errors(fit)[1]
end

function plot_thermo_trans(mm_limits,mm_errors,mp_limits,mp_errors,t_values)
	fig = figure()
	errorbar(t_values,mm_limits,yerr=[mm_errors,mm_errors],label="Minus")
	errorbar(t_values,mp_limits,yerr=[mp_errors,mp_errors],label="Plus")
	xlabel("Hopping Strength")
	ylabel("Chemical Strength")
	title("Thermodynamic Limit of Mott Transition")
	legend()
	return
end

function get_thermodynamic_transitions(mm_data,mp_data,t_values,which_counts; kwargs...)
	number_sizes = length(which_counts)
	mm_limits = [0.0 for i in 1:length(t_values)]
	mm_errors = [0.0 for i in 1:length(t_values)]
	mp_limits = [0.0 for i in 1:length(t_values)]
	mp_errors = [0.0 for i in 1:length(t_values)]
	for i in 1:length(t_values)
		ts = t_values[i]
		mm_vals = [mm_data[j][i] for j in 1:number_sizes]
		mp_vals = [mp_data[j][i] for j in 1:number_sizes]
		mm_results = get_linfit(1 ./ which_counts,mm_vals)
		mp_results = get_linfit(1 ./ which_counts,mp_vals)
		mm_limits[i] = mm_results[1]
		mm_errors[i] = mm_results[2]
		mp_limits[i] = mp_results[1]
		mp_errors[i] = mp_results[2]
	end
	
	if get(kwargs, :if_plot, true)
		plot_thermo_trans(mm_limits,mm_errors,mp_limits,mp_errors,t_values)
	end
	return mm_limits,mm_errors,mp_limits,mp_errors,t_values
end

function get_mag_string(no_magF)
	if no_magF
		return "No Mag"
	else
		return "Mag On"
	end	
end

function check_overlaps(all_ttns)
	count = length(all_ttns)
	overlaps = zeros(count,count)
	for i in 1:count
		for j in 1:count
			if i == j
				overlaps[i,j] = 1.0
				continue
			end
			overlaps[i,j] = localinner(all_ttns[i],all_ttns[j])
		end
	end
	fig = figure()
	imshow(overlaps)
	colorbar()
	return overlaps
end

function get_part_counts_range_fillings(site_count,phi)
	starting_count = 0.25 * phi * site_count
	if starting_count <= 1.0
		println("Lattice too small")
		return
	else
		starting_count = Int(starting_count)
	end
	ending_count = Int(0.75 * phi * site_count)
	num_points = ending_count - starting_count
	part_counts = [starting_count + i for i in 0:num_points]
	nu_values = [part_counts[i]/(phi*site_count) for i in 1:length(part_counts)]
	return part_counts,nu_values
end
#=

#final_time = 0.1
if_per = false
mag_off = true
evolve = true
chemical = false
mu = 0.5
#max_occupation = 3
bc_string = get_periodic_title_string(if_per)
mag_string = get_mag_string(mag_off)
expan = TTNKit.DefaultExpander(0.5)
#us = 1.0
ts = 0.01
nu = 1/2
layers = 6
tot_sites = 2^layers
edge_sites = Int(sqrt(2^layers))
alpha = 1/edge_sites
num_particles = Int(edge_sites/2)#get_particles_needed(layers; nu=nu)#tot_sites - 
mdim = 70
nswps = 3

println("Using $num_particles particles on $tot_sites sites")


og_ttn, hamilt, dm_sp = build_full_harperhofstadter(layers,num_particles,ts,nu; max_dim=mdim, num_sweeps=nswps,phi=alpha, if_periodic=if_per,max_occ=1,if_sweep=evolve,sweep_type="dmrg",expander=expan,if_chem=chemical,chem_strength=mu,no_magF=mag_off,output_level=0)

rez1 = get_occupancy(dm_sp.ttn)
rez2 = get_current_yfunc(dm_sp.ttn)
rez3 = get_ydir_greenfunc(dm_sp.ttn)
rez4 = get_xdir_greenfunc(dm_sp.ttn)
=#
#=
og_ttn_pin, hamilt_pin, dm_sp_pin = build_full_harperhofstadter(layers,num_particles,ts,nu; max_dim=mdim, num_sweeps=nswps,phi=alpha, if_periodic=if_per,max_occ=1,if_pinning_pot=true,if_sweep=evolve,sweep_type="dmrg",expander=expan,if_chem=chemical,chem_strength=mu,no_magF=mag_off,output_level=0)

rez4 = get_occupancy(dm_sp_pin.ttn; plot_title="Pinning")
rez5 = get_current_yfunc(dm_sp_pin.ttn; plot_title="Pinning")
rez6 = get_ydir_greenfunc(dm_sp_pin.ttn; plot_title="Pinning")

fig = figure()
imshow(rez1 - rez4)
colorbar()
title("Difference btw Pinning Potential")
=#

#=rez6 = get_xdir_greenfunc(Int(sqrt(2^layers)),dm_sp.ttn; plot_title="$bc_string")
rez5 = get_xdir_greenfunc(Int(sqrt(2^layers)),dm_sp.ttn; direction="reverse",plot_title="$bc_string")
rez4 = get_current_xfunc(Int(sqrt(2^layers)),dm_sp.ttn; plot_title="$bc_string")

=#
#end
#
#all_paths = get_all_sites_paths_and_plot(dm_sp.ttn,edge_sites; likely_path=true)

























"fin"

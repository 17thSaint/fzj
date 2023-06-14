using PyPlot
include("../review-practice-codes/ttn.jl")

function long_range_scaling(x_final,virt_edge_length,initial_strength; kwargs...)
	if_plot = get(kwargs, :if_plot, false)
	if_save_data = get(kwargs, :if_save_data, false)
	if_save_fig = get(kwargs, :if_save_fig, false)
	if_hard_cutoff = get(kwargs, :cliff, false)
	if_rounding = get(kwargs, :rounding, true)
	if if_hard_cutoff
		if_rounding = false
	end
	final_minimum = get(kwargs, :limit, 10^-3)
	trunc_minimum = get(kwargs, :trunc_min, 10^-6)
	trunc = get(kwargs, :trunc, trunc_minimum*initial_strength)
	scaling_func = get(kwargs, :scaling, "flat")
	
	strengths = zeros(virt_edge_length)
	
	if scaling_func == "flat"
		strengths[1:x_final+1] .= initial_strength
	elseif scaling_func == "exp"
		strengths = map(1:virt_edge_length) do x
			initial_strength * exp(-log(1/final_minimum)*(x-1)/x_final)	
		end
	elseif scaling_func == "power"
		println("Not done power")
	end
	
	if if_hard_cutoff
		strengths[x_final + 2:end] .= 0.0
	elseif if_rounding
		final_index = findfirst(x -> x .<= trunc,strengths)
		if !isnothing(final_index)
			strengths[final_index:end] .= 0.0
		end
	end
	
	if_plot ? plot_long_range_scaling(strengths,virt_edge_length; kwargs...) : nothing
	if_save_data ? save_long_range_scaling(strengths,virt_edge_length; kwargs...) : nothing

	return strengths
end

function save_long_range_scaling(strengths,virt_edge_length; kwargs...)
	filename = get(kwargs, :name, "scaling-strength")
	location = get(kwargs, :location, pwd())
	xcoord = [i for i in 0:virt_edge_length-1]
	scaling_data = Dict([("strengths",strengths),("xcoord",xcoord)])
	write_data_hdf5(filename,scaling_data,location=location; kwargs...)
	return
end

function plot_long_range_scaling(strengths,virt_edge_length; kwargs...)
	fig = figure()
	xcoord = range(0,virt_edge_length-1,virt_edge_length)
	plot(xcoord,strengths,"-p")
	xlabel("Distance of Interaction")
	ylabel("Strength")
	title("Interaction Strength Distance Scaling")
	
	if get(kwargs, :if_save_fig, false)
		filename = get(kwargs, :name, "scaling-strength")
		save_figure(filename; kwargs...)
	end
	return
end

function build_HH_net(num_layers; kwargs...)
	conserve_qns = get(kwargs, :syms, true)
	max_occ = 1
	
	net = TTNKit.BinaryRectangularNetwork(num_layers, TTNKit.ITensorNode, "Boson";conserve_qns=conserve_qns,dim=max_occ+1)
	
	return net
end

function get_interaction_coords(given_site,inter_dist,lat) # written by ChatGPT 12.06.2023 then vastly edited 13.06.2023 by me
	virtual, physical = given_site
	coordinates = []
    
	virt_edge_length, phys_edge_length = size(lat)
	#if typeof(given_site) == Int64
	#	given_site = TTNKit.coordinate(lat,given_site)
	#end

	for shift in [-inter_dist % virt_edge_length,inter_dist % virt_edge_length]
		new_virtual = virtual + shift
				
		# Apply periodic boundary conditions along the virtual-axis
		if new_virtual < 1
			new_virtual += virt_edge_length
		elseif new_virtual > virt_edge_length
			new_virtual -= virt_edge_length
		end

               	# Check if new coordinates are within lattice dimensions
		if 1 <= new_virtual <= virt_edge_length && new_virtual != virtual
 			append!(coordinates, [[new_virtual,physical]])
 		else
 			println("Still Outside Lattice")
		end

	end
	return coordinates
end

function long_range_HH_ham(net,t_strength,phi; kwargs...)
	resulting_ham = []
	phys_edge_length,virt_edge_length = get_lattice_dims(net)
	
	u_strength = get(kwargs, :u_strength, 1.0)
	scaling_distance = get(kwargs, :scaling_dist, 0) #[[u_strength/2]; [0 for i in 1:virt_edge_length-1]]
	long_range_strengths = long_range_scaling(scaling_distance,virt_edge_length,u_strength; kwargs...)
	
	if_interaction = !all(long_range_strengths.==0)
	if_periodic = get(kwargs, :if_periodic, true)
	if_hopping = get(kwargs, :if_hopping, true)
	if_chem = get(kwargs, :if_chem, false)
	if_pinning_pot = get(kwargs, :if_pinning_pot, false)
	vpinning = get(kwargs, :vpinning, 2.5)
	no_magF = get(kwargs, :no_magF, false)
	chem_strength = get(kwargs, :chem_strength, 0.0)
	
	lat = TTNKit.physical_lattice(net)
	
	if if_hopping
		hopping = TTNKit.OpSum()
		for (s1,s2) in TTNKit.nearest_neighbours(lat,collect(1:TTNKit.number_of_sites(lat)); periodic=if_periodic)
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)
						
			coeff = get_inter_coeff(s1_coord,s2_coord,t_strength,phi,phys_edge_length,virt_edge_length; kwargs...)
			hopping += (coeff,"Adag",s1_coord,"A",s2_coord)
			hopping += (conj(coeff),"Adag",s2_coord,"A",s1_coord)
		end
		append!(resulting_ham,[hopping])
	end
	
	if if_interaction
		interaction = TTNKit.OpSum()
		for i in 1:length(long_range_strengths)
			if long_range_strengths[i] == 0
				continue
			else
				local_strength = long_range_strengths[i]
				if i == 1
					for j in TTNKit.eachindex(lat)
						interaction += (local_strength/2,"Adag * Adag * A * A",TTNKit.coordinate(lat,j))
					end
					append!(resulting_ham,[interaction])
				else
					for j in TTNKit.eachindex(lat)
						interaction_sites = get_interaction_coords(TTNKit.coordinate(lat,j),i-1,lat)
						for k in interaction_sites
							interaction += (local_strength/2,"Adag * A",TTNKit.coordinate(lat,j),"Adag * A",Tuple(k))
							interaction -= (local_strength/2,"Adag * A",TTNKit.coordinate(lat,j))
						end
					end
					append!(resulting_ham,[interaction])
				end
			end
		end
	end
	
	if if_chem && chem_strength != 0.0
		chem = TTNKit.OpSum()
		for i in TTNKit.eachindex(lat)
			chem -= (chem_strength,"N",TTNKit.coordinate(lat,i))
		end
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
	
	if length(resulting_ham) > 1
		return sum(resulting_ham)
	else
		return resulting_ham[1]
	end
end

function get_densdens_corrs(ttn,distances; kwargs...)
	phys_edge_length,virt_edge_length = get_lattice_dims(ttn)
	direction = get(kwargs, :direction, "virt")
	dim_dict = Dict([("virt",virt_edge_length),("phys",phys_edge_length)])
	chosen_dim = dim_dict[direction]
	other_direction = collect(keys(dim_dict))[findfirst(x -> x != direction,collect(keys(dim_dict)))]
	other_dim = dim_dict[other_direction]
	lat = TTNKit.physical_lattice(TTNKit.network(ttn))
	densdens_corr = zeros(length(distances),other_dim)
	corr_errors = zeros(length(distances),other_dim)
	for j in 1:length(distances)
		distance = distances[j]
		for k in 1:other_dim
			all_values = zeros(chosen_dim)
			for i in 1:chosen_dim
				next_site = i + (distance % chosen_dim)
				if next_site > chosen_dim
					next_site = next_site % chosen_dim
					if next_site == 0
						next_site = i
					end
				end
				if direction == "phys"
					pos1 = TTNKit.linear_ind(lat,(k,i))
					pos2 = TTNKit.linear_ind(lat,(k,next_site))
				elseif direction == "virt"
					pos1 = TTNKit.linear_ind(lat,(i,k))
					pos2 = TTNKit.linear_ind(lat,(next_site,k))
				else
					println("Bad direction")
					return
				end
				value = TTNKit.correlation(ttn,"N","N",pos1,pos2)
				all_values[i] = real(value)
			end
			densdens_corr[j,k] = mean(all_values)
			corr_errors[j,k] = std(all_values)
		end
	end
	if get(kwargs, :if_plot, true)
		title_string = "DensDens Corr, " * get(kwargs, :plot_title, "$direction Edge Count = $chosen_dim")
		fig = figure()
		for i in 1:length(distances)
			plot([j for j in 1:other_dim],densdens_corr[i,:],"-p",label="$(distances[i])")
		end
		yscale("log")
		title(title_string)
		xlabel("$other_direction Axis")
		ylabel("Corr")
		legend()
	end
	return densdens_corr,corr_errors
end


#=
function check_duplicates(file_name)
	file_string,file_type = split(file_name,".")
	if file_name in readdir()
		println("Found Duplicate File, renaming")
	end
	while file_name in readdir()
		string_elems = split(file_string,"-")
		if "mk" in string_elems
			mk_number_loc = findfirst(x -> x == "mk",string_elems) + 1
			string_elems[mk_number_loc] = string(parse(Int,string_elems[mk_number_loc]) + 1)
		else
			append!(string_elems,["mk","2"])
		end
		file_string = join(string_elems,"-")
		file_name = file_string * "." * file_type
	end
	return file_name
end

function save_figure(file_name,file_location=pwd())
	current_location = pwd()
	cd(file_location)
	file_name = check_duplicates(file_name * ".png")
	savefig(file_name,format="png")
	cd(current_location)
	return
end
=#

function get_allAVG_densdenscorr(ttn,distances; kwargs...)
	all_corrs,all_errors = get_densdens_corrs(ttn,distances; kwargs...,if_plot=false)
	avg_corrs = [mean(all_corrs[i,:]) for i in 1:length(distances)]
	avg_errors = [mean(all_errors[i,:]) for i in 1:length(distances)]
	
	if_save_fig = get(kwargs, :if_save_fig, false)
	if_save_data = get(kwargs, :if_save_data, false)
	
	if_plot = get(kwargs, :if_plot, true)
	if if_plot
		plot_title = get(kwargs, :plot_title, "")
		title_string = "AVG DensDens Corr, " * plot_title
		plot_allAVG_densdenscorr(distances,avg_corrs,avg_errors,title_string; kwargs...)
	end
	
	if_save_data ? save_allAVG_densdenscorr(distances,avg_corrs,avg_errors; kwargs...) : nothing
	
	return avg_corrs,avg_errors,distances
end

function plot_allAVG_densdenscorr(distances,avg_corrs,avg_errors,title_string; kwargs...)
	fig = figure()
	errorbar(distances,avg_corrs,yerr=[avg_errors,avg_errors])
	yscale("log")
	title(title_string)
	xlabel("Distance")
	ylabel("AVG Corr")
	
	if get(kwargs, :if_save_fig, false)
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "densdens")
		save_figure(filename;kwargs...)
	end
	return
end

function save_allAVG_densdenscorr(distances,avg_corrs,avg_errors; kwargs...)
	location = get(kwargs, :location, pwd())
	filename = get(kwargs, :name, "densdens")
	avg_data_dict = Dict([("dists",distances),("vals",avg_corrs),("errors",avg_errors)])
	write_data_hdf5(filename,avg_data_dict,location)
end

function get_mdim(num_layers,shift=(false,0.5))
	if num_layers <= 4
		maxdims = 50
	elseif 5 <= num_layers <= 6
		 maxdims = 100
	elseif 7 <= num_layers <= 8
		maxdims = 150
	end
	if shift[1]
		maxdims *= 1 + shift[2]
	end
	return Int(round(maxdims,digits=0))
end

#
if_per = false
mag_off = true
evolve = true
chemical = false
mu = 0.5
#max_occupation = 3
expan = TTNKit.DefaultExpander(0.5)
ts = 0.001
nu = 1/2
layers = 4#parse(Int,ARGS[1])
tot_sites = 2^layers
if layers % 2 == 0
	edge_sites = Int(sqrt(2^layers))
	num_particles = Int(edge_sites/2)
else
	num_particles = Int(sqrt(2^(layers+1))/2)
end
if !mag_off
	alpha = num_particles/(mu * (tot_sites))
else
	alpha = 0.0
end
mdim = get_mdim(layers,(true,1))
nswps = 3
println("Using $num_particles particles on $tot_sites sites")

loc = pwd()

if_cliff = true
sc_type = "flat"
limit = 0.5
dists = [i for i in 1:2*edge_sites]

net = build_HH_net(layers; syms=true)
all_ttns = []
#for i in 1:3
	longrange_dist = 1#parse(Int, ARGS[1])
	title_string = "LR $sc_type = $longrange_dist, md = $mdim"

	ham = long_range_HH_ham(net,ts,alpha; scaling=sc_type,limit=limit,scaling_dist=longrange_dist,cliff=if_cliff,if_periodic=if_per,if_chem=chemical,no_magF=mag_off)

	og_ttn, hamilt, dm_sp = build_full_harperhofstadter(layers,num_particles,ts,nu; ttn_net=net,ham_op=ham,max_dim=mdim, num_sweeps=nswps,phi=alpha, if_periodic=if_per,max_occ=1,if_sweep=evolve,sweep_type="dmrg",expander=expan,if_chem=chemical,chem_strength=mu,no_magF=mag_off,output_level=0)
	#append!(all_ttns,[dm_sp.ttn])
	#rez = get_densdens_corrs(dm_sp.ttn,dists;plot_title=title_string)
	#rez2 = get_densdens_corrs(dm_sp.ttn,dists;plot_title=title_string*" Phys",direction="phys")
	datafile_name = "densdens-test-md-$mdim-n-$tot_sites-lr-$longrange_dist"
	rez3 = get_allAVG_densdenscorr(dm_sp.ttn,dists;if_plot=true,plot_title=title_string, if_save_fig=true, if_save_data=true, name=datafile_name,location=loc)
	
	#rez1 = get_occupancy(dm_sp.ttn; plot_title=title_string, if_save=true, name="occ-test-md-$mdim-lr-$longrange_dist.png")
	#
	#rez2 = get_current_yfunc(dm_sp.ttn)
	#rez3_fqh = get_ydir_greenfunc(dm_sp.ttn)
	#rez3_sf = get_ydir_greenfunc(dm_sp.ttn; plot_title=title_string)
	#rez4 = get_xdir_greenfunc(dm_sp.ttn; plot_title=title_string)
	#
#end
#




























"fin"

#using PyPlot
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
		strengths[1] = initial_strength
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
	
	if_plot || if_save_fig ? plot_long_range_scaling(strengths,virt_edge_length; kwargs...) : nothing
	if_save_data ? save_long_range_scaling(strengths,virt_edge_length; kwargs...) : nothing

	return strengths
end

function save_long_range_scaling(strengths,virt_edge_length; kwargs...)
	filename = get(kwargs, :name, "scaling-strength")
	location = get(kwargs, :location, pwd())
	xcoord = [i for i in 0:virt_edge_length-1]
	scaling_data = Dict([("strengths",strengths),("xcoord",xcoord)])
	write_data_jld2(filename,scaling_data,location=location; kwargs...)
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

function get_allAVG_densdenscorr(ttn,distances; kwargs...)
	all_corrs,all_errors = get_densdens_corrs(ttn,distances; kwargs...,if_plot=false)
	avg_corrs = [mean(all_corrs[i,:]) for i in 1:length(distances)]
	avg_errors = [mean(all_errors[i,:]) for i in 1:length(distances)]
	
	if_save_fig = get(kwargs, :if_save_fig, false)
	if_save_data = get(kwargs, :if_save_data, false)
	if_plot = get(kwargs, :if_plot, true)

	if_save_data ? save_allAVG_densdenscorr(distances,avg_corrs,avg_errors; kwargs...) : nothing
	if_plot || if_save_fig ? plot_allAVG_densdenscorr(distances,avg_corrs,avg_errors; kwargs...) : nothing
	
	return avg_corrs,avg_errors,distances
end

function plot_allAVG_densdenscorr(distances,avg_corrs,avg_errors; kwargs...)
	data_dict = get(kwargs, :data_dict, nothing)
	if !isnothing(data_dict)
		distances,avg_corrs,avg_errors = data_dict["dists"], data_dict["vals"], data_dict["errors"]
	end
	plot_title = get(kwargs, :plot_title, "")
	title_string = "AVG DensDens Corr, " * plot_title
	fig = figure()
	#errorbar(distances,avg_corrs,yerr=[avg_errors,avg_errors])
	plot(distances,avg_corrs,"-p")
	yscale("log")
	title(title_string)
	xlabel("Distance")
	ylabel("AVG Corr")
	
	if get(kwargs, :if_save_fig, false)
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "densdens")
		filename = check_plot_label(filename,"densdens")
		save_figure(filename;kwargs...)
	end
	return
end

function save_allAVG_densdenscorr(distances,avg_corrs,avg_errors; kwargs...)
	location = get(kwargs, :location, pwd())
	filename = get(kwargs, :name, "densdens")
	metadata = get(kwargs, :metadata, nothing)
	avg_data_dict = Dict([("dists",distances),("vals",avg_corrs),("errors",avg_errors)])
	write_data_jld2(filename,avg_data_dict,location,metadata)
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

function bulk_density(ttn,bulk_width=1)
	occ_mat = get_occupancy(ttn;if_plot=false)
	num_particles = sum(occ_mat)
	bulk_occ_mat = occ_mat[1+bulk_width:end-bulk_width,1+bulk_width:end-bulk_width]
	bulk_density = sum(bulk_occ_mat)/prod(size(bulk_occ_mat))
	return bulk_density
end

function deriv_bulk_dens(ttn1,ttn2,alpha_change,bulk_width=1)
	bulk_dens_1 = bulk_density(ttn1,bulk_width)
	bulk_dens_2 = bulk_density(ttn2,bulk_width)
	deriv = (bulk_dens_1 - bulk_dens_2)/alpha_change
	return deriv
end
#=

#params_dict = Dict([("layers",5),("mdim",70),("nn_strength",0.0),("alpha",0.1),("mag_off",false),("lr",0)])
# usually in params: mag_off, layers, mdim, longrange_dist
params_dict = make_args_dict(ARGS)
if_change = get(params_dict, "if_change", false)
change = get(params_dict, "change", 0.0001)
limit = get(params_dict, "nn_strength", 1.0)
layer_count = get(params_dict, "layers", 4)
mag_off = get(params_dict, "mag_off", true)
mdim = get(params_dict, "mdim", get_mdim(layer_count,(false,1)))
longrange_dist = get(params_dict, "lr", 0)
if if_change in ["pos","neg"]
	if if_change == "pos"
		alpha = get(params_dict, "alpha", nothing) + change
		params_dict["alpha"] = round(alpha,digits=4)
	elseif if_change == "neg"
		alpha = get(params_dict, "alpha", nothing) - change
		params_dict["alpha"] = round(alpha,digits=4)
	end
else
	alpha = get(params_dict, "alpha", nothing)
end
if layer_count % 2 == 0
	edge_sites = Int(sqrt(2^layer_count))
	num_particles = get(params_dict, "particles", Int(edge_sites/2))
else
	edge_sites = Int(sqrt(2^(layer_count+1)))
	num_particles = get(params_dict, "particles", Int(sqrt(2^(layer_count+1))/2))
end

#
sweep_type = "dmrg"
max_occ = 1
if_per = false
evolve = true
chemical = false
mu = 0.5
#max_occupation = 3
expan = TTNKit.DefaultExpander(0.5)
ts = 0.001
nu = 1/2
tot_sites = 2^layer_count

if isnothing(alpha)
	if !mag_off
		alpha = num_particles/(mu * (tot_sites))
	else
		alpha = 0.0
	end
end
nswps = 3
#

plotting = false
save_plot = false
save_data = true

loc = "/local/geraghty/cluster-data"
if_cliff = true
sc_type = "exp"
dists = [i for i in 1:2*edge_sites]
lr_scaling = long_range_scaling(longrange_dist,edge_sites,1.0; cliff=if_cliff,limit=limit,scaling=sc_type,if_plot=false)

#
metadata_dict = Dict([("if_per",if_per),("mag_off",mag_off),("chemical",chemical),("mu",mu),("ts",ts),("nu",nu),("layers",layer_count),("particles",num_particles),("alpha",alpha),("mdim",mdim),("nswps",nswps),("if_cliff",if_cliff),("sc_type",sc_type),("longrange_dist",longrange_dist),("max_occ",max_occ),("sweep_type",sweep_type),("limit",limit),("lr_scaling",lr_scaling),("if_change",if_change),("change",change)])

if length(keys(params_dict)) == 0
	datafile_name = "layers-$layer_count-particles-$num_particles-mdim-$mdim-mag-$(!mag_off)-lr-$longrange_dist"
else
	datafile_name = make_parameters_filename(params_dict)
end

#
println(datafile_name)
title_string = "Np = $num_particles, LR = $longrange_dist at $limit"
println("Starting Script using $num_particles particles on $tot_sites sites with $(!mag_off) Mag Field, Bond Dim = $mdim, and Long Range Dist = $longrange_dist")
starting = time()
net = build_HH_net(layer_count; syms=true)
ham = long_range_HH_ham(net,ts,alpha; scaling=sc_type,limit=limit,scaling_dist=longrange_dist,cliff=if_cliff,if_periodic=if_per,if_chem=chemical,no_magF=mag_off)
og_ttn, hamilt, dm_sp = build_full_harperhofstadter(layer_count,num_particles,ts,nu; ttn_net=net,ham_op=ham,if_save_data=save_data,name="ttn-"*datafile_name,location=loc,metadata=metadata_dict,max_dim=mdim, num_sweeps=nswps,phi=alpha, if_periodic=if_per,max_occ=max_occ,if_sweep=evolve,sweep_type=sweep_type,expander=expan,if_chem=chemical,chem_strength=mu,no_magF=mag_off,output_level=0)
total_time = time() - starting
println("Running time = $total_time")
=#

























"fin"

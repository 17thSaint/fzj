using Pkg
Pkg.activate(".")
include("../review-practice-codes/ttn.jl")
using Profile,MKL

function spin_matrix_element(m1,m2,spin,direction::String)
	if direction == "X"
		return ((m1 == m2+1) + (m1+1 == m2)) * 0.5 * sqrt(spin*(spin+1) - m1*m2)
	elseif direction == "Y"
		return ((m1 == m2+1) - (m1+1 == m2)) * 0.5 * sqrt(spin*(spin+1) - m1*m2)
	elseif direction == "Z"
		return (m1 == m2) * m1 * -1
	end
end
	
function spin_matrix_full(spin,direction::String)
	spin_mat = zeros(Int(2*spin+1),Int(2*spin+1))
	for i in -spin:spin
		for j in -spin:spin
			spin_mat[Int(i+spin+1),Int(j+spin+1)] = spin_matrix_element(j,i,spin,direction)
		end
	end
	return spin_mat
end

function get_magnetization(wavefunc,spin_value::Float64,direction::String; kwargs...)
	phys_dim_length = typeof(wavefunc) == MPS ? length(wavefunc) : get_lattice_dims(wavefunc)[1]
	magnetization = [0.0*im for i in 1:phys_dim_length]
	sites = [i for i in 1:length(magnetization)]
	for j in 1:Int(2*spin_value)
		for k in 1:Int(2*spin_value)
			spin_matrix_value = spin_matrix_element(k-spin_value-1,j-spin_value-1,spin_value,direction)
			#=
			if spin_matrix_value != dict_spin_version
				println("Different values: ",dict_spin_version,", ",spin_matrix_value,". Probably should be ",spin_matrix_element(j,k,spin_value,direction))
				spin_matrix_value = dict_spin_version
			end
			=#
			if spin_matrix_value != 0.0
				if typeof(wavefunc) != MPS
					part = round(TTNKit.correlation(wavefunc,"Adag","A",j,k) .* spin_matrix_value,digits=8)
				else
					part = round.(expect(wavefunc,"Cr$(j) * Anh$(k)") .* spin_matrix_value,digits=8)
				end
				if direction == "Y"
					magnetization .+= part / im
				else
					magnetization .+= part
				end
			end
		end
	end
	
	if_plot = get(kwargs, :if_plot, true)
	if_plot ? plot_magnetization(magnetization,sites,direction; kwargs...) : nothing
	
	return magnetization,sites
end

function plot_magnetization(mags,sites,direction; kwargs...)
	title_string = "Magnetization $direction, " * get(kwargs, :plot_title, "")
	plot_label = get(kwargs, :plot_label, "")
	isempty(plot_label) ? fig = figure() : nothing
	plot(sites,mags,"-p",label=plot_label)
	legend()
	xlabel("Site Number")
	ylabel("M" * direction)
	title(title_string)
end

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
	elseif scaling_func == "lr_flat"
		strengths[1] = initial_strength
		strengths[2:x_final+1] .= final_minimum
	elseif scaling_func == "rydberg"
		blockade_radius = initial_strength
		strengths = map(0:virt_edge_length-1) do x
			1.0 * (blockade_radius^6) / (blockade_radius^6 + x^6)
		end
		x_final < length(strengths) ? strengths[x_final+2:end] .= 0.0 : nothing
	end
	
	if if_hard_cutoff
		strengths[x_final + 2:end] .= 0.0
	elseif if_rounding
		final_index = findfirst(x -> abs(x) .<= trunc,strengths)
		if !isnothing(final_index)
			strengths[final_index:end] .= 0.0
		end
	end
	
	if_plot || if_save_fig ? plot_long_range_scaling(strengths,virt_edge_length; kwargs...) : nothing
	#if_save_data ? save_long_range_scaling(strengths,virt_edge_length; kwargs...) : nothing

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
	if_fermion = get(kwargs, :if_fermion, false)
	particle_type = if_fermion ? "Fermion" : "Boson"
	max_occ = get(kwargs,:max_occ,1)
	
	net = if_fermion ? TTNKit.BinaryRectangularNetwork(num_layers, TTNKit.ITensorNode, particle_type;conserve_nf=conserve_qns,conserve_nfparity=false) : TTNKit.BinaryRectangularNetwork(num_layers, TTNKit.ITensorNode, particle_type;conserve_qns=conserve_qns,dim=max_occ+1)
	
	return net
end

# isotropic not implemented for cylinder
function get_interaction_coords(given_site,inter_dist,lat,if_periodic_virt,if_anis) # written by ChatGPT 12.06.2023 then vastly edited 13.06.2023 by me
	virtual, physical = given_site
	coordinates = []
    
	virt_edge_length, phys_edge_length = size(lat)
	#if typeof(given_site) == Int64
	#	given_site = TTNKit.coordinate(lat,given_site)
	#end

	for shift in [-inter_dist % virt_edge_length,inter_dist % virt_edge_length]
		new_virtual = virtual + shift
				
		# Apply periodic boundary conditions along the virtual-axis
		if if_periodic_virt
			if new_virtual < 1
				new_virtual += virt_edge_length
			elseif new_virtual > virt_edge_length
				new_virtual -= virt_edge_length
			end
		end

        # Check if new coordinates are within lattice dimensions
		if 1 <= new_virtual <= virt_edge_length && new_virtual != virtual
 			append!(coordinates, [[physical,new_virtual]])
 		#else
 			#println("Still Outside Lattice")
		end

	end

	if !if_anis
		for shift in [-inter_dist % phys_edge_length,inter_dist % phys_edge_length]
			new_physical = physical + shift
			if if_periodic_virt
				if new_physical < 1
					new_physical += phys_edge_length
				elseif new_physical > phys_edge_length
					new_physical -= phys_edge_length
				end
			end

			if 1 <= new_physical <= phys_edge_length && new_physical != physical
				append!(coordinates, [[new_physical,virtual]])
			end
		end
	end

	return unique(coordinates)
end

function long_range_HH_ham(net,t_strength,phi; kwargs...)
	resulting_ham = []
	phys_edge_length,virt_edge_length = get_lattice_dims(net)
	println("Phys = ",phys_edge_length,", Virt = ",virt_edge_length)
	
	scaling_distance = get(kwargs, :scaling_dist, 0)
	
	restricted_size = get(kwargs, :restricted_size, [phys_edge_length,virt_edge_length])
	if_periodic_virt = get(kwargs, :if_periodic_virt, false)
	if_periodic_phys = get(kwargs, :if_periodic_phys, false)
	if_per = [if_periodic_phys,if_periodic_virt]
	if_hopping = get(kwargs, :if_hopping, true)
	if_nn_int = get(kwargs, :if_nn_int, false)
	onsite_strength = get(kwargs, :onsite_strength, 0.0)
	if_pinning_pot = get(kwargs, :if_pinning_pot, false)
	vpinning = get(kwargs, :vpinning, 2.5)
	no_magF = get(kwargs, :no_magF, false)
	chem_strength = get(kwargs, :chem_strength, 0.0)
	centralflux_strength = get(kwargs, :centralflux_strength, 0.0)
	twist_angle = get(kwargs, :twist_angle, 0.0)
	#hopping_anisotropy = get(kwargs, :hopping_anisotropy, 1.0) t_phys / t_synth = anisotropy
	
	long_range_strengths = long_range_scaling(scaling_distance,virt_edge_length,onsite_strength; kwargs...)
	#long_range_strengths[1] = 0.0
	display(long_range_strengths)
	if_interaction = !all(long_range_strengths.==0)
	
	lat = TTNKit.physical_lattice(net)
	
	hopping_old = get(kwargs, :hopping_old, false)
	if if_hopping && hopping_old
		if_periodic_phys ? nothing : centralflux_strength = 0.0
		hopping = TTNKit.OpSum()
		#
		for (s1,s2) in TTNKit.nearest_neighbours(lat,collect(1:TTNKit.number_of_sites(lat)))
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)

			if any([s1_coord[1],s2_coord[1]] .> restricted_size[1])
				continue
			end
			if any([s1_coord[2],s2_coord[2]] .> restricted_size[2])
				continue
			end
						
			coeff = get_inter_coeff(s1_coord,s2_coord,t_strength,phi,phys_edge_length,virt_edge_length; kwargs...)
			
			if s1_coord[2] == s2_coord[2]
				if s1_coord[1] > s2_coord[1]
					coeff *= exp(im*2*pi*centralflux_strength/size(lat)[1])
				else
					coeff *= exp(-im*2*pi*centralflux_strength/size(lat)[1])
				end
			end
			
			hopping += (coeff,"Adag",s1_coord,"A",s2_coord)
			hopping += (conj(coeff),"Adag",s2_coord,"A",s1_coord)
			
		end
		#
		if if_periodic_virt
			for i in 1:restricted_size[1]
				s2_coord = (i,1)
				s1_coord = (i,restricted_size[2])
				coeff = get_inter_coeff(s1_coord,s2_coord,t_strength,phi,phys_edge_length,virt_edge_length; kwargs...)
				hopping += (coeff,"Adag",s1_coord,"A",s2_coord)
				hopping += (conj(coeff),"Adag",s2_coord,"A",s1_coord)
			end
		end

		if if_periodic_phys
			for i in 1:restricted_size[2]
				s2_coord = (1,i)
				s1_coord = (restricted_size[1],i)
				coeff = get_inter_coeff(s1_coord,s2_coord,t_strength,phi,phys_edge_length,virt_edge_length; kwargs...)
				#coeff *= exp(im*2*pi*centralflux_strength/size(lat)[1])
				#coeff *= exp(im*twist_angle*2*pi)
				hopping += (coeff,"Adag",s1_coord,"A",s2_coord)
				hopping += (conj(coeff),"Adag",s2_coord,"A",s1_coord)
			end
		end

		append!(resulting_ham,[hopping])
	else

	#if if_hopping
		hopping = TTNKit.OpSum()
		for s_phys in 1:restricted_size[1]
			for s_synth in 1:restricted_size[2]
				starting_site = [s_phys,s_synth]
				for which_axis in [1,2]
						ending_site = starting_site .+ ((which_axis == 1,which_axis == 2))

						# enforce boundary conditions
						if ending_site[which_axis] > restricted_size[which_axis]
							if if_per[which_axis]
								ending_site[which_axis] = 1
							else
								continue
							end
						end

						if ending_site[which_axis] < 1
							if if_per[which_axis]
								ending_site[which_axis] = restricted_size[which_axis]
							else
								continue
							end
						end

						coeff = round(get_inter_coeff(starting_site,ending_site,t_strength,phi,phys_edge_length,virt_edge_length; kwargs...),digits=10)
						#dir == -1 ? coeff = conj(coeff) : nothing
						hopping += (coeff,"Adag",Tuple(starting_site),"A",Tuple(ending_site))
						hopping += (conj(coeff),"Adag",Tuple(ending_site),"A",Tuple(starting_site))
				end
			end
		end
		append!(resulting_ham,[hopping])
	end
	
	if if_interaction
		if kwargs[:scaling] == "rydberg"
			if_anis = false
		else
			if_anis = true
		end
		interaction = TTNKit.OpSum()
		for (idx,stren) in enumerate(long_range_strengths)
			if stren == 0.0
				continue
			else
				if idx == 1
					#=for j in TTNKit.eachindex(lat)
						s_coord = TTNKit.coordinate(lat,j)
						if s_coord[1] > restricted_size[1] || s_coord[2] > restricted_size[2]
							continue
						end
						interaction += (stren,"N * N",s_coord)
						interaction -= (stren,"N",s_coord)
					end=#
					continue
				else
					for j in TTNKit.eachindex(lat)
						s_coord = TTNKit.coordinate(lat,j)
						if s_coord[1] > phys_edge_length || s_coord[2] > virt_edge_length
							continue
						end
						interaction_sites = get_interaction_coords(s_coord,idx-1,lat,if_periodic_virt,if_anis)
						
						for k in interaction_sites
							if k[1] > restricted_size[1] || k[2] > restricted_size[2]
								continue
							end
							interaction += (stren,"Adag * A",s_coord,"Adag * A",Tuple(k))
						end
					end
				end
			end
		end
		append!(resulting_ham,[interaction])
	end
	
	if restricted_size != [phys_edge_length,virt_edge_length]
		restrict_size = TTNKit.OpSum()
		for i in restricted_size[1]+1:phys_edge_length
			for j in 1:virt_edge_length
				restrict_size += (1e10,"N",(i,j))
			end
		end
		for i in restricted_size[2]+1:virt_edge_length
			for j in 1:restricted_size[1]
				restrict_size += (1e10,"N",(j,i))
			end
		end
		append!(resulting_ham,[restrict_size])
	end

	if chem_strength != 0.0
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

function long_range_HH_ham(metadata::Dict)
	net = metadata["net"]
	t_strength = metadata["t_strength"]
	phi = metadata["phi"]
	model_paras = dict_to_symbols(metadata)
	return long_range_HH_ham(net,t_strength,phi; model_paras...)
end

function get_densdens_corrs(ttn::TTNKit.TreeTensorNetwork,distances; kwargs...)
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

function radial_box_dist(ttn)
	occ_mat = get_occupancy(ttn;if_plot=false)
	num_particles = sum(occ_mat)
	range_rad = Int(maximum(size(occ_mat))/2)
	rads = [i for i in 1:range_rad]
	neg_rads = sort(-1 .* rads)
	full_rads = [neg_rads; rads]
	vals_top = zeros(range_rad)
	vals_bot = zeros(range_rad)
	for i in 1:range_rad
		change = i - 1
		if i == 1
			vals_top[i] = sum(occ_mat[range_rad,range_rad:range_rad + 1])
			vals_bot[i] = sum(occ_mat[range_rad + 1,range_rad:range_rad + 1])
		else
			top_vert_left = sum(occ_mat[range_rad - change:range_rad,range_rad - change])
			top_vert_right = sum(occ_mat[range_rad - change:range_rad,range_rad + change + 1])
			top_horiz = sum(occ_mat[range_rad - change,range_rad - change + 1:range_rad + change])
			vals_top[i] = top_vert_left + top_vert_right + top_horiz
			
			vals_bot[i] = sum(occ_mat[range_rad + 1:range_rad + change + 1,range_rad - change]) + sum(occ_mat[range_rad + 1:range_rad + change + 1,range_rad + change + 1]) + sum(occ_mat[range_rad + change + 1,range_rad - change + 1:range_rad + change])
		end
	end
	combined_vals = [reverse(vals_top); vals_bot] ./ num_particles
	return full_rads,combined_vals
end

function bulk_density(ttn::TreeTensorNetwork,bulk_width_phys=1,bulk_width_virt=1; kwargs...)
	if isnothing(ttn)
		occ_mat = get(kwargs, :occ_mat, nothing)
	else
		occ_mat = get_occupancy(ttn;if_plot=false)
	end
	size(occ_mat)[1] == size(occ_mat)[2] ? bulk_width_virt = bulk_width_phys : nothing
	num_particles = sum(occ_mat)
	bulk_occ_mat = occ_mat[1+bulk_width_phys:end-bulk_width_phys,1+bulk_width_virt:end-bulk_width_virt]
	bulk_density = sum(bulk_occ_mat)/prod(size(bulk_occ_mat))
	return bulk_density
end

function deriv_bulk_dens(ttn1,ttn2,alpha_change,bulk_width_phys=1,bulk_width_virt=1; kwargs...)
	if isnothing(ttn1) && isnothing(ttn2)
		occ_mat1 = get(kwargs, :occ_mat1, nothing)
		occ_mat2 = get(kwargs, :occ_mat2, nothing)
		bulk_dens_1 = bulk_density(nothing,bulk_width_phys,bulk_width_virt; occ_mat=occ_mat1)
		bulk_dens_2 = bulk_density(nothing,bulk_width_phys,bulk_width_virt; occ_mat=occ_mat2)
	else
		bulk_dens_1 = bulk_density(ttn1,bulk_width_phys,bulk_width_virt)
		bulk_dens_2 = bulk_density(ttn2,bulk_width_phys,bulk_width_virt)
	end
	deriv = (bulk_dens_1 - bulk_dens_2)/alpha_change
	return deriv
end

function get_position_pairs(phys_length,virt_length)
	position_pairs = []
	reverse_position_pairs = []
	for i in 1:phys_length
		for j in 1:virt_length
			left_site = i + phys_length*(j-1)
			right_site = left_site+1
			if right_site > phys_length*j
				right_site = phys_length*(j-1) + 1
			end
			append!(position_pairs,[(left_site,right_site)])
			append!(reverse_position_pairs,[(right_site,left_site)])
		end
	end
	return sort(position_pairs)
end

function reorder_vector_to_matrix(left_moving,right_moving)
	side_length = Int(sqrt(length(left_moving)))
	left_moving_mat = zeros(side_length,side_length) .* im
	right_moving_mat = zeros(side_length,side_length) .* im
	for i in 1:side_length
			left_moving_mat[i,:] = left_moving[side_length*(i-1)+1:side_length*i]
			right_moving_mat[i,:] = right_moving[side_length*(i-1)+1:side_length*i]
	end
	return right_moving_mat .+ left_moving_mat
end

function ttn_current_site(psi::TreeTensorNetwork,virt_site; kwargs...)
	centralflux_strength = get(kwargs, :centralflux_strength, 0.0)
	phys_length,virt_length = get_lattice_dims(psi)
	phys_left = (virt_site-1)*phys_length + Int(phys_length/2)
	coeff = im*2*pi*exp(im*2*pi*centralflux_strength/phys_length)
	left_moving = conj(coeff)*TTNKit.correlation(psi,"Adag","A",phys_left,phys_left+1)
	right_moving = coeff*TTNKit.correlation(psi,"Adag","A",phys_left+1,phys_left)
	return right_moving + left_moving
end

function ttn_current(psi::TreeTensorNetwork; kwargs...)
	centralflux_strength = get(kwargs, :centralflux_strength, 0.0)
	phys_length,virt_length = get_lattice_dims(psi)
	position_pairs = get_position_pairs(phys_length,virt_length)
	right_moving = [0.0*im for i in 1:length(position_pairs)]
	left_moving = [0.0*im for i in 1:length(position_pairs)]
	coeff = im*2*pi*exp(im*2*pi*centralflux_strength/phys_length)
	for i in 1:length(position_pairs)
		left_site,right_site = position_pairs[i]
		left_moving[i] = TTNKit.correlation(psi,"Adag","A",left_site,right_site)
		right_moving[i] = TTNKit.correlation(psi,"Adag","A",right_site,left_site)
	end
	total_current = coeff*sum(right_moving) + conj(coeff)*sum(left_moving)
	return total_current,reorder_vector_to_matrix(conj(coeff) .* left_moving,coeff .* right_moving)
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

    return sqrt(dx^2 + dy^2),(dx,dy)
end

function physical_distance_correlation(psi::TreeTensorNetwork; kwargs...)
	if_plot = get(kwargs, :if_plot, false)
	if_periodic_phys = get(kwargs, :if_periodic_phys, true)
	if_periodic_virt = get(kwargs, :if_periodic_virt, false)
	densmat = get(kwargs, :densmat, nothing)

	lat = TTNKit.physical_lattice(TTNKit.network(psi))

	phys_length,virt_length = get_lattice_dims(psi)
	all_corrs = [[] for i in 1:virt_length]
	dists = [[] for i in 1:virt_length]
	for s=1:virt_length
		for j=1:phys_length
			top = phys_length#mod1(j+Int(phys_length/2),phys_length)
			#println("Doing site ",(j,s))
			for jj=1:top#phys_length
				if isnothing(densmat)
					corr_val = TTNKit.correlation(psi,"Adag","A",(j,s),(jj,s))
					corr_val /= sqrt(TTNKit.expect(psi,"N",(j,s)) * TTNKit.expect(psi,"N",(jj,s)))
					#corr_val += conj(corr_val)
				else
					corr_val = densmat[TTNKit.linear_ind(lat,(j,s)),TTNKit.linear_ind(lat,(jj,s))]
					normalization = sqrt(densmat[TTNKit.linear_ind(lat,(j,s)),TTNKit.linear_ind(lat,(j,s))] * densmat[TTNKit.linear_ind(lat,(jj,s)),TTNKit.linear_ind(lat,(jj,s))]) 
					corr_val /= normalization
					#corr_val += conj(corr_val)
				end
				dist_btw = abs(j-jj)#find_dist((s,j),(s,jj),(virt_length,phys_length),(if_periodic_virt,if_periodic_phys))
				if dist_btw in dists[s]
					append!(all_corrs[s][findfirst(x -> x == dist_btw,dists[s])],abs(corr_val))
				else
					append!(dists[s],[dist_btw])
					sort!(dists[s])
					insert!(all_corrs[s],findfirst(x -> x == dist_btw,dists[s]),[abs(corr_val)])
				end
			end
		end
	end

	corrs = [[] for i in 1:virt_length]
	for i in 1:virt_length
		corrs[i] = ([mean(all_corrs[i][j]) for j in 1:length(all_corrs[i])])
	end

	if if_periodic_phys
		middle = Int(phys_length/2) + 1
		corr_lengths = correlation_length(dists[1][2:middle],[corrs[i][2:middle] for i in 1:length(corrs)]; if_plot=if_plot)
	else
		middle = phys_length - 1
		corr_lengths = correlation_length(dists[1][2:middle],[corrs[i][2:middle] for i in 1:length(corrs)]; if_plot=if_plot)
	end

	if if_plot
		#fig = figure()
		title_string = get(kwargs, :plot_title, "")
		for i in 1:virt_length
			plot(dists[i],abs.(corrs[i]),"p",label="$i")
			yscale("log")
		end
		legend()
		xlabel("Distance")
		ylabel("Correlation")
		title(title_string)
	end

	return dists,corrs,corr_lengths
end

function correlation_length(dists,phys_correlations; kwargs...)
	if_plot = get(kwargs, :if_plot, false)

	exp_fit(x,p) = p[1].* exp.(-x ./ p[2]) .+ p[3]

	all_fits = [curve_fit(exp_fit,dists,phys_correlations[i],[1.0,1.0,0.0]) for i in 1:length(phys_correlations)]
	corr_lengths = [all_fits[i].param[2] for i in 1:length(all_fits)]

	if if_plot
		#fig = figure()
		for i in 1:length(phys_correlations)
			#scatter(dists,phys_correlations[i],label="Site $i")
			plot(dists,exp_fit(dists,all_fits[i].param),"-",label="$(round(corr_lengths[i],digits=4))")
		end
		xlabel("Distance")
		ylabel("Correlation")
		title("Correlation Lengths")
		legend()
	end
	return corr_lengths
end

function distance_correlation(psi::TreeTensorNetwork; kwargs...)
	if_plot = get(kwargs, :if_plot, true)
	if_periodic_phys = get(kwargs, :if_periodic_phys, true)
	if_periodic_virt = get(kwargs, :if_periodic_virt, false)
	densmat = get(kwargs, :densmat, nothing)

	lat = TTNKit.physical_lattice(TTNKit.network(psi))

	phys_length,virt_length = get_lattice_dims(psi)
	all_corrs = []
	dists = []
	for s=1:virt_length, j=1:phys_length
		println(round(100*s/(virt_length),digits=2),"%")
		for ss=1:virt_length, jj=1:phys_length
			if isnothing(densmat)
				corr_val = TTNKit.correlation(psi,"Adag","A",(s,j),(ss,jj)) + TTNKit.correlation(psi,"Adag","A",(ss,jj),(s,j))
			else
				corr_val = densmat[TTNKit.linear_ind(lat,(s,j)),TTNKit.linear_ind(lat,(ss,jj))]
				corr_val += conj(corr_val)
			end
			dist_btw = find_dist((s,j),(ss,jj),(virt_length,phys_length),(if_periodic_virt,if_periodic_phys))
			if dist_btw in dists
				append!(all_corrs[findfirst(x -> x == dist_btw,dists)],corr_val)
			else
				append!(dists,[dist_btw])
				sort!(dists)
				insert!(all_corrs,findfirst(x -> x == dist_btw,dists),[corr_val])
			end
		end
	end

	corrs = ([sum(all_corrs[i]) for i in 1:length(all_corrs)])
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

function momentum_occupation(psi::TreeTensorNetwork,p_count::Int,p_end::Real,direction="phys"; kwargs...)
	if_neg = get(kwargs, :if_neg, true)
	if_save_data = get(kwargs, :if_save_data, false)
	if_plot = p_count != 1 ? get(kwargs, :if_plot, false) : false
	p_start = get(kwargs, :p_start, 0.0)
	other_p = get(kwargs, :other_p, 0.0)
	densmat = get(kwargs, :densmat, nothing)
	if_fermion = get(kwargs, :if_fermion, false)
	creation = if_fermion ? "Cdag" : "Adag"
	annihilation = if_fermion ? "C" : "A"

	lat = TTNKit.physical_lattice(TTNKit.network(psi))

	if if_neg
		p_start = -p_end
	end

	phys_length,virt_length = get_lattice_dims(psi)

	momenta = range(p_start,stop=p_end,length=p_count)
	mom_occs = zeros(p_count) .* im
	
	for s=1:virt_length, j=1:phys_length
		for ss=1:virt_length, jj=1:phys_length
			if isnothing(densmat)
				corr_val = TTNKit.correlation(psi,creation,annihilation,(s,j),(ss,jj))
				corr_val += conj(corr_val)
			else
				corr_val = densmat[TTNKit.linear_ind(lat,(s,j)),TTNKit.linear_ind(lat,(ss,jj))]
				corr_val += conj(corr_val)
			end
			for (i,p) in enumerate(momenta)
				pvec = direction == "phys" ? [other_p,p] : [p,other_p]
				mom_occs[i] += corr_val * exp(im*pi*dot(pvec,[s-ss,j-jj]))
			end
		end
	end

	if_plot ? plot_momentum_occupation(momenta,abs2.(mom_occs); kwargs...) : nothing

	return momenta,mom_occs
end

function momentum_occupation(psi::TreeTensorNetwork,p_count::Int,p_end::Real; kwargs...)
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

	if_plot ? plot_momentum_occupation(momenta,real.(mom_occs); kwargs...) : nothing

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

function is_within_lattice(point::Tuple{Int,Int}, lattice_size::Tuple{Int,Int})
    return 1 <= point[1] <= lattice_size[1] && 1 <= point[2] <= lattice_size[2]
end

function loop_sites(starting_site,which_quadrant,phys_length,virt_length; kwargs...)
	if_periodic_phys = get(kwargs, :if_periodic_phys, false)
	if_periodic_virt = get(kwargs, :if_periodic_virt, false)
	loop_length = get(kwargs, :loop_length, 1)

    # Determine the direction of the loop based on the quadrant
    dx = (which_quadrant in [1,4]) ? loop_length : -loop_length
    dy = (which_quadrant in [1,2]) ? loop_length : -loop_length

    # The second point is one step to the right or left
    first_x_move = if_periodic_phys ? (mod1(starting_site[1] + dx,phys_length), starting_site[2]) : (starting_site[1] + dx, starting_site[2])

    # The third point is one step up or down from the second point
    first_y_move = if_periodic_virt ? (mod1(first_x_move[1],virt_length), first_x_move[2] + dy) : (first_x_move[1], first_x_move[2] + dy)

    # The fourth point is one step to the left or right from the third point
    second_x_move = if_periodic_phys ? (mod1(first_y_move[1] - dx,phys_length), first_y_move[2]) : (first_y_move[1] - dx, first_y_move[2])

	resulting_loop = [starting_site,first_x_move,first_y_move,second_x_move]

	# check that the loop is inside the lattice
	if any([!is_within_lattice(p,(phys_length,virt_length)) for p in resulting_loop])
		#println("Loop is outside lattice")
		return nothing
	else
    	# Return the four pointsYour location
    	return resulting_loop
	end
end

function closed_loop(psi::TreeTensorNetwork, starting_site; kwargs...)
	phys_length,virt_length = get_lattice_dims(psi)
	which_direction = get(kwargs, :direction, 1)
	loop_length = get(kwargs, :loop_length, 1)
	if_fermion = get(kwargs, :if_fermion, false)
	creation = if_fermion ? "Cdag" : "Adag"
	annihilation = if_fermion ? "C" : "A"

	# find loop to use by trying all directions
	sites_to_loop = loop_sites(starting_site,which_direction,phys_length,virt_length; kwargs...)
	while isnothing(sites_to_loop) && which_direction < 4
		which_direction += 1
		sites_to_loop = loop_sites(starting_site,which_direction,phys_length,virt_length; kwargs...)
	end
	if isnothing(sites_to_loop)
		println("No loop found")
		return nothing
	end
		
	calced_values = zeros(length(sites_to_loop)) .* im
	for (idx,s) in enumerate(sites_to_loop)
		next_site = idx == length(sites_to_loop) ? sites_to_loop[1] : sites_to_loop[idx+1]
		calced_values[idx] = TTNKit.correlation(psi,creation,annihilation,next_site,s)
	end

	return angle(prod(calced_values)),calced_values,sites_to_loop
end

function cdw_structure_factor(rho,qvec::Tuple,psi::TreeTensorNetwork; kwargs...)
	if_periodic_phys = get(kwargs, :if_periodic_phys, false)
	if_periodic_synth = get(kwargs, :if_periodic_synth, false)

	lat = TTNKit.physical_lattice(TTNKit.network(psi))
	phys_len,synth_len = size(lat)[1],size(lat)[2]

	occs = get_occupancy(psi; if_plot=false,densmat=rho)

	struc_fact = 0.0
	for j in 1:phys_len
		for s in 1:synth_len
			p1 = (j,s)
			p1_linear = TTNKit.linear_ind(lat,p1)
			for jj in 1:phys_len
				for ss in 1:synth_len
					p2 = (jj,ss)
					p2_linear = TTNKit.linear_ind(lat,p2)
					dist = find_dist(p1, p2, (phys_len,synth_len), (if_periodic_phys,if_periodic_synth))[2]
					struc_fact += occs[p1[1],p1[2]] * occs[p2[1],p2[2]] * exp(im * dot(qvec,dist))
				end
			end
		end
	end
	return struc_fact / sum(occs)
end

function cdw_struct_full(rho,psi::TreeTensorNetwork,howmany=100,qmax=3.0; kwargs...)
	if_plot = get(kwargs, :if_plot, true)

	qs = range(-qmax,stop=qmax,length=howmany)
	struct_factor = zeros(ComplexF64,howmany,howmany)
	for (i,qx) in enumerate(qs)
		for (j,qy) in enumerate(qs)
			struct_factor[i,j] = cdw_structure_factor(rho,(qx,qy),psi; kwargs...)
		end
	end

	if if_plot
		fig = figure()
		imshow(abs.(struct_factor),extent=[qs[1],qs[end],qs[1],qs[end]])
		colorbar()
		xlabel("qx")
		ylabel("qy")
		title("CDW Structure Factor")
	end

	return struct_factor,qs
end

function distance_correlation(rho::Matrix,wavefunc::TreeTensorNetwork,Lx::Int64,Ly::Int64,direction::String="x")
    #=if layers % 2 == 0.0
		Lx = Int(sqrt(2^layers))
		Ly = Int(sqrt(2^layers))
	else
		Lx = Int(sqrt(2^(layers+1)))
		Ly = Int(sqrt(2^(layers-1)))
	end=#
	lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))

    if direction == "x"
        len = Lx
        other_len = Ly
    else
        len = Ly
        other_len = Lx
    end
    dist_corrs = zeros(Float64,len-1)
    corr_counts = zeros(Int64,len-1)

    for x1 in 1:len
        for y1 in 1:other_len
            s1 = direction == "x" ? TTNKit.linear_ind(lat,(x1,y1)) : TTNKit.linear_ind(lat,(y1,x1))
            for x2 in 1:len-1
                if x1 == x2
                    continue
                end
                s2 = direction == "x" ? TTNKit.linear_ind(lat,(x2,y1)) : TTNKit.linear_ind(lat,(y1,x2))
                dist_corrs[Int(abs(x1-x2))] += abs(rho[s1,s2])
                corr_counts[Int(abs(x1-x2))] += 1
            end
        end
    end

    dist_corrs ./= corr_counts

    return dist_corrs
end

function rydberg_2pcorr(rho::Matrix; kwargs...)
	if_plot = get(kwargs, :if_plot, true)
	occs = get_occupancy(rho; if_plot=false)

	dist_corrs::Dict{Float64,Float64} = Dict()
	dist_counts::Dict{Float64,Int64} = Dict()
	for x1 in 1:size(occs,1)
		for y1 in 1:size(occs,2)
			for x2 in 1:size(occs,1)
				for y2 in 1:size(occs,1)
					dist_btw = round(sqrt((x2 - x1)^2 + (y2-y1)^2),digits=4)
					if dist_btw in keys(dist_corrs)
						dist_corrs[dist_btw] += 1.0
						dist_counts[dist_btw] += 1
					else
						dist_corrs[dist_btw] = 1.0
						dist_counts[dist_btw] = 1
					end
					if (x1,y1) == (x2,y2)
						dist_corrs[0.0] -= occs[x1,y1] / (occs[x1,y1]*occs[x2,y2])
					end
				end
			end
		end
	end

	for k in keys(dist_corrs)
		dist_corrs[k] /= dist_counts[k]
	end

	if if_plot
		plot_title = get(kwargs,:plot_title,"")
		fig = figure()
		scatter(collect(keys(dist_corrs)),collect(values(dist_corrs)))
		xlabel("Distance")
		ylabel("Two Particle Correlation")
		title(plot_title)
	end

	return dist_corrs
end

function rydberg_2pcorr(wavefunc::TreeTensorNetwork; kwargs...)
	if_plot = get(kwargs, :if_plot, true)

	site_count = Int(2^TTNKit.number_of_layers(wavefunc))
	coords = TTNKit.physical_coordinates(TTNKit.network(wavefunc))

	onsite_occs = abs.(TTNKit.expect(wavefunc,"N"))

	dist_corrs::Dict{Float64,Vector{Float64}} = Dict()
	#dist_counts::Dict{Float64,Int64} = Dict()
	for x1 in 1:Int(sqrt(site_count))
		for y1 in 1:Int(sqrt(site_count))
			s1 = findfirst(i -> (x1,y1) == i,coords)
			for x2 in 1:Int(sqrt(site_count))
				for y2 in 1:Int(sqrt(site_count))
					s2 = findfirst(i -> (x2,y2) == i,coords)

					dist_btw = round(sqrt((x2 - x1)^2 + (y2-y1)^2),digits=4)
					if dist_btw in keys(dist_corrs)
						append!(dist_corrs[dist_btw],[abs(TTNKit.correlation(wavefunc,"N","N",s1,s2))])
					else
						dist_corrs[dist_btw] = [abs(TTNKit.correlation(wavefunc,"N","N",s1,s2))]
					end
					
					if (x1,y1) == (x2,y2)
						dist_corrs[0.0][end] -= onsite_occs[x1,y1]
					end

					dist_corrs[dist_btw][end] /= onsite_occs[x1,y1] * onsite_occs[x2,y2]
				end
			end
		end
	end

	avg_dist_corrs::Dict{Float64,Float64} = Dict()
	for (k,v) in dist_corrs
		avg_dist_corrs[k] = mean(v)
	end


	if if_plot
		plot_title = get(kwargs,:plot_title,"")
		fig = figure()
		scatter(collect(keys(avg_dist_corrs)),collect(values(avg_dist_corrs)))
		xlabel("Distance")
		ylabel("Two Particle Correlation")
		title(plot_title)
	end

	return avg_dist_corrs
end

function average_close_keys(dict, bin_width)
	# Create bins
	bins = Dict()

	for (key, value) in dict
		# Find the bin this key belongs to
		bin_key = round(key / bin_width) * bin_width

		# If the bin doesn't exist, create it
		if !haskey(bins, bin_key)
			bins[bin_key] = []
		end

		# Add the value to the bin
		push!(bins[bin_key], value)
	end

	# Now average the values in each bin
	averaged_dict = Dict()
	for (key, values) in bins
		averaged_dict[key] = mean(values)
	end

	return averaged_dict
end

function check_fluxes(alpha,Lx::Int64,Ly::Int64,if_periodic_x::Bool,if_periodic_y::Bool)
    if alpha == 0.0
        return nothing
    end
    x_shift,y_shift = !if_periodic_x, !if_periodic_y
    num_fluxes = round(alpha*(Lx - x_shift) * (Ly - y_shift),digits=5)
    println("Number of Fluxes = ",num_fluxes," for Lx = ",Lx," and Ly = ",Ly)
    if !isinteger(num_fluxes)
        error("Number of fluxes is not an integer")
    end

    if if_periodic_x && !isinteger(num_fluxes/Lx)
        error("Number of fluxes is not an integer multiple of Lx")
    end

    if if_periodic_y && !isinteger(num_fluxes/Ly)
        error("Number of fluxes is not an integer multiple of Ly")
    end

    return nothing
end


# Usage:
# averaged_dict = average_close_keys(your_dict, 0.1)

#= Momentum occupation testing
lnet = build_HH_net(4; syms=true)
states = fill("0", 16)
old_ttn = TTNKit.ProductTreeTensorNetwork(lnet,states)
ttn = initialize_ttn(old_ttn,50,1)
freeboson_ham = long_range_HH_ham(lnet,1.0,0.0)
old, hamilt, dm_sp = find_ground_state(4,1,1.0; ham_op=freeboson_ham,ttn_net=lnet,seed_ttn=ttn,sweep_type="dmrg",output_level=1,mdim=50,num_sweeps=10,if_save_data=false)
fb_gs = dm_sp.ttn
fb_occ_mat = get_occupancy(fb_gs)
=#
	


#
if false

#nnst = 0.0
#layers = 6
#=layers = 4
num_parts = 2
ref_dict = Dict([("layers",layers),("particles",num_parts)])
loc = get_folder_location("cluster-data/synth-dims")
all_files = find_data_file(ref_dict,"ttn",loc)
alphas = zeros(length(all_files))
for (idx,file) in enumerate(all_files)
	alphas[idx] = get_params_dict_from_filename(file)["alpha"]
end=#

#layers = 6
#lr = 7
#anises = [0.01,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.6,0.8,0.9,1.1,1.3,1.5,1.7,1.9,2.0,2.5,3.0,3.5,4.0,6.0,8.0,9.0,10.0,15.0,20.0,25.0,30.0,40.0,50.0,70.0,90.0,100.0,1000.0,10000.0]
#anises = range(1.0,5.0,length=5)
#strens = range(0.0,5.0,length=10)
#alphas = [4/(0.5*64)]#range(4/(0.2*64),4/(0.8*64),length=20)
#strens = range(0.1,0.5,length=3)
#for (idx,anis) in enumerate(anises)
#for (idx,stren) in enumerate(strens)
	params_dict = Dict([("hopping_anisotropy",1.0),("make_smaller_lattice",[8,8]),("nrgtol",5e-5),("particles",4),("layers",6),("mdim",100),("if_save_data",true),("filling",0.5),("onsite_strength",10.0),("lr","all"),("if_periodic_phys",false),("if_periodic_virt",false)])
	# usually in params: mag_off, layers, mdim, longrange_dist
	#params_dict = make_args_dict(ARGS)
	open_cores = get(params_dict, "open_cores", 5)
	if typeof(open_cores) != String
		BLAS.set_num_threads(open_cores)	
		display(BLAS.get_config())
	end
	#true
	nrgtol = get(params_dict, "nrgtol", 1E-4)
	cutoff = get(params_dict, "cutoff", 0.0)
	if_NN = get(params_dict, "if_nn_int", false)
	if_pinning = get(params_dict, "if_pinning", false)
	if_gpu = get(params_dict, "if_gpu", false)
	if_change = get(params_dict, "if_change", false)
	if_densmat = get(params_dict, :if_densmat, true)
	change = get(params_dict, "change", 0.0001)
	onsite_strength = get(params_dict, "onsite_strength", 0.0)
	anis = get(params_dict, "hopping_anisotropy", 1.0)
	layer_count = Int(get(params_dict, "layers", 4))
	mag_off = get(params_dict, "mag_off", true)
	mdim = get(params_dict, "mdim", 300)
	longrange_dist = get(params_dict, "lr", 0)
	centralflux_strength = get(params_dict, "centralflux_strength", 0.0)
	twist_angle = get(params_dict, "twist_angle", 0.0)
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
	#
	if layer_count % 2 == 0
		edge_sites = Int(sqrt(2^layer_count))
		phys_edge_length,synth_edge_length = edge_sites,edge_sites
		num_particles = get(params_dict, "particles", Int(edge_sites/2))
	else
		edge_sites = Int(sqrt(2^(layer_count+1)))
		phys_edge_length,synth_edge_length = edge_sites,Int(edge_sites/2)
		num_particles = get(params_dict, "particles", Int(sqrt(2^(layer_count+1))/2))
	end

	make_smaller_lattice = get(params_dict, "make_smaller_lattice", [phys_edge_length,synth_edge_length])
	if make_smaller_lattice != [phys_edge_length,synth_edge_length]
		phys_edge_length,synth_edge_length = make_smaller_lattice
		edge_sites = phys_edge_length
	end

	if_per_phys = get(params_dict, "if_periodic_phys", true)
	if_per_virt = get(params_dict, "if_periodic_virt", false)
	if isnothing(alpha)
		filling = get(params_dict, "filling", 1.0)
		phys_shift,synth_shift = !if_per_phys,!if_per_virt
		alpha = num_particles/(filling*(phys_edge_length - phys_shift)*(synth_edge_length - synth_shift))
		mag_off = false
	else
		mag_off = alpha == 0.0
	end
	check_fluxes(alpha,phys_edge_length,synth_edge_length,if_per_phys,if_per_virt)
	#

	#
	#nu = 1.0
	sweep_type = "dmrg"
	max_occ = get(params_dict, "max_occ", 1)
	
	evolve = true
	#max_occupation = 3
	expan = TTNKit.DefaultExpander(1.0)#TTNKit.NoExpander()
	herenoise = [0.0]
	ts = 1.0
	mu = get(params_dict, "chem_strength", 0.0)
	tot_sites = 2^layer_count
	syms = get(params_dict, "syms", true)
	#=
	nu = 1.0
	if isnothing(alpha)
		if !mag_off
			alpha = num_particles/(nu * (tot_sites))
		else
			alpha = 0.0
		end
	end
	=#
	nswps = 100
	#alpha = 7/64
	#num_particles = Int(sqrt(2^layer_count)/2)#Int(alpha * tot_sites * nu)
	

	plotting = false
	save_plot = false
	save_data = get(params_dict, "if_save_data", true)
	if_cluster = any([occursin("local",pwd()),occursin("Local",pwd()),occursin("geraghty",pwd())])
	if_continuous_saving = get(params_dict,"if_continuous_saving",if_cluster || layer_count >= 7)
	save_data ? nothing : if_continuous_saving = false

	
	if_cliff = false
	trunc = get(params_dict,"trunc",1e-3)
	sc_type = "flat"
	dists = [i for i in 1:2*edge_sites]
	longrange_dist == "all" ? longrange_dist = edge_sites-1 : nothing
	lr_scaling = long_range_scaling(longrange_dist,edge_sites,onsite_strength; cliff=if_cliff,scaling=sc_type,if_plot=false,trunc=trunc)

	if if_per_phys && if_per_virt
		dataloc = get_folder_location("cluster-data/synth-dims/torus")
	elseif if_per_phys || if_per_virt
		dataloc = get_folder_location("cluster-data/synth-dims")
	elseif !if_per_phys && !if_per_virt
		dataloc = get_folder_location("cluster-data/synth-dims/obc")
	end
	if sc_type == "rydberg"
		dataloc = get_folder_location("cluster-data/synth-dims/rydberg")
	end
	loc = get(params_dict, "dataloc", dataloc)
	
	#counting = 50
	#strens = range(num_particles/(0.2*tot_sites),num_particles/(2.0*tot_sites),length=counting) #range(0.02,0.25,length=counting)
	#sforderparams = zeros(counting)
	#=alpha_start = 0.0525
	alpha_end = 0.0725
	alpha_count = 5
	alphas = [alpha_start + (i-1)*(alpha_end-alpha_start)/(alpha_count-1) for i in 1:alpha_count] .- change/2
	alphas = [alphas; alphas .+ change]
	#
	alpha_center = mag_off ? 0.0 : 1 * num_particles / tot_sites
	wavefuncs = []
	currents = []
	nrgs = []
	#display(alphas)
	counting = 50
	centralflux_strength = 0.0
	#parts = [i for i in 1:Int(tot_sites/2)]
	#fillings = range(0.2,3.0,length=counting)
	strens = range(num_particles/(0.2*tot_sites),num_particles/(3.0*tot_sites),length=counting) #range(0.02,0.25,length=counting)
	centermoms = [0.0 for i in 1:counting]# .* im
	=#
	#for (idx,alpha) in enumerate(strens)
	#for (idx,num_particles) in enumerate(parts)
		#alpha = 0.0
		filename_dict = Dict([("layers",layer_count),("lr",longrange_dist),("particles",num_particles),("alpha",round(alpha,digits=4)),("if_periodic_phys",if_per_phys),("onsite_strength",onsite_strength),("hopping_anisotropy",anis)])
		if make_smaller_lattice != [sqrt(2^layer_count),sqrt(2^layer_count)]
			filename_dict["make_smaller_lattice"] = phys_edge_length
		end
		twist_angle != 0.0 ? filename_dict["twist_angle"] = twist_angle : nothing
		#if length(keys(params_dict)) == 0
		#	datafile_name = "layers-$layer_count-particles-$num_particles-mdim-$mdim-mag-$(!mag_off)-lr-$longrange_dist"
		#else
			datafile_name = make_parameters_filename(filename_dict)
		#end
		model_paras = (hopping_anisotropy=anis,
						syms=syms,
						cutoff=cutoff,
						twist_angle=twist_angle,
						if_continuous_saving=if_continuous_saving,
						nrgtol=nrgtol,
						if_densmat=if_densmat,
						restricted_size=make_smaller_lattice,
						centralflux_strength=centralflux_strength,
						if_pinning_pot=if_pinning,
						if_periodic_phys=if_per_phys,
						if_periodic_virt=if_per_virt,
						if_nn_int=if_NN,
						nn_int_strength=lr_scaling[2],
						chem_strength=mu,
						no_magF=mag_off,
						scaling=sc_type,
						scaling_dist=longrange_dist,
						onsite_strength=onsite_strength,
						cliff=if_cliff,
						trunc=trunc,
						if_change=if_change,
						change=change,
						if_gpu=if_gpu,
						noise=herenoise,
						if_save_data=save_data,
						if_save_fig=save_plot,
						if_sweep=evolve,
						sweep_type=sweep_type,
						expander=expan,
						max_occ=max_occ,
						mdim=mdim,
						num_sweeps=nswps,
						phi=alpha,
						output_level=0,
						name="ttn-"*datafile_name,
						location=loc)
		
		metadata_dict = merge(named_tuple_to_dict(model_paras),filename_dict)

		#
		println(datafile_name)
		if_exists,found_data = check_data_exists(filename_dict,"ttn"; location=loc,output_level=false)

		if if_exists
			println("Found Data")
			wavefunc = found_data[1]["ttn"]
			dens = found_data[1]["densmat"]

			rezobs = found_data[2]["observer"]
			ham = found_data[2]["ham"]
			#append!(wavefuncs,[wavefunc])
		else
			#title_string = "Np = $num_particles, LR = $longrange_dist at $limit"
			println("Starting Script using $num_particles particles on $tot_sites sites with $(!mag_off) Mag Field, Bond Dim = $mdim, and Long Range Dist = $longrange_dist")
			#if true
			#densmat = nothing
			starting = time()
			net = build_HH_net(layer_count; syms=syms, max_occ=max_occ)
			ham = long_range_HH_ham(net,ts,alpha; model_paras...)
			#display(ham)
			og_ttn, hamilt, dm_sp, rezobs, runtime, dens = find_ground_state(layer_count,num_particles; ttn_net=net,ham_op=ham,model_paras...,metadata=merge(metadata_dict,Dict([("ham",ham),("net",net),("t_strength",ts)])))
			total_time = time() - starting
			println("Running time = $total_time")
			wavefunc = dm_sp.ttn
			#append!(wavefuncs,[dm_sp.ttn])
		end

		#imshow(real.(dens))
		#colorbar()

		#scatter(stren,rezobs.nrg[end],c="b")
		#xlabel("Interaction Strength")
		#ylabel("Energy")

		#Profile.print()

		#=
		scatter([anis],[sum(dens) / (tot_sites * num_particles)],c="b")
		xlabel("Hopping Anisotropy")
		ylabel("Zero Momentum Occupation")
		xscale("log")
		=#

		#dcorrs = distance_correlation(dens,wavefunc,make_smaller_lattice[1],make_smaller_lattice[2],"y")
		#display(dcorrs)
		occs = get_occupancy(wavefunc; densmat=dens, plot_title="TTN")
		#rydberg_2pcorr(wavefunc)
		#=plot(collect(1:Int(sqrt(2^layer_count))),occs[4,:],label="$(round(num_particles/(alpha*tot_sites),digits=4))")
		legend()
		xlabel("Sites")
		ylabel("Occupancy")=#
		#=
		fig = figure()
		physical_distance_correlation(wavefunc; densmat=dens,if_plot=true,if_periodic_phys=if_per_phys,if_periodic_virt=if_per_virt)
		title("Filling = $(round(num_particles/(alpha*tot_sites),digits=4))")
		get_occupancy(wavefunc; plot_title = "Filling = $(round(num_particles/(alpha*tot_sites),digits=4))",densmat=dens)
		=#

		#occs = get_occupancy(wavefunc; plot_title="$stren")
		#densities[idx] = sum(occs) / tot_sites
		#scatter([mu],[densities[idx]],c="b")
		#xlabel("Chemical Potential")
		#ylabel("Density")

		#=sforderparams[idx] = abs(2*sum(dens)) / (2^layer_count)

		if idx > 1
			plot(num_particles ./ ((2^layer_count) .* [strens[idx-1],strens[idx]]),[sforderparams[idx-1],sforderparams[idx]],"-p",c="b")
		end=#

		#physical_distance_correlation(wavefunc; densmat=dens,if_plot=true)

		#get_occupancy(wavefunc)

		#append!(currents,[[ttn_current_site(dm_sp.ttn,i; centralflux_strength=centralflux_strength) for i in 1:edge_sites]])
		#append!(nrgs,[dm_sp.current_energy])

		#momentum_occupation(wavefuncs[idx],50,1.0; if_plot=true)

		
		#centermoms[idx] = minimum(abs.(rez[2]))

		#=
		if false
		allmoms = momentum_occupation(wavefunc,1,0.0; densmat=densmat)
		centermoms[idx] = allmoms[2][1]
		if idx > 1
			plot([num_particles/(strens[idx-1]*tot_sites),num_particles/(alpha*tot_sites)],[centermoms[idx-1],centermoms[idx]],"-p",c="b")
			#plot([(num_particles-1)/(strens[idx-1]*tot_sites),num_particles/(alpha*tot_sites)],[centermoms[idx-1],centermoms[idx]],"-p",c="b")
		else
			scatter([num_particles/(alpha*tot_sites)],[centermoms[idx]],c="b")
			#scatter([num_particles/(strens[idx]*tot_sites)],[centermoms[idx]],c="b")
		end
		end
		=#
		#get_occupancy(dm_sp.ttn; plot_title = "Alpha = $(round(alpha,digits=4))")
		#get_greenfunc(dm_sp.ttn,"phys")
		#get_greenfunc(dm_sp.ttn,"virt")
		#=
		specs = entanglement_spectrum(dm_sp.ttn)
		fig = figure()
		scatter(collect(1:mdim),-log.(specs))
		=#
	#end
end

#
#plot(strens,real.(centermoms),"-p")

#=
fig = figure()
for j in 1:edge_sites
plot(strens,[real(currents[i][j]) for i in 1:counting],"-p",label="$j")
legend()
end

fig2 = figure()
plot(strens,nrgs,"-p")
xlabel("Central Flux Strength")
ylabel("Energy")
=#
#=
include("time_evolution.jl")

psi_gs = wavefuncs[1]
if_states = true

alloccs = [get_occupancy(psi_gs; if_plot=false)]
whatmeasuring = Dict([("occs",(get_occupancy,(if_plot=false,))),("currents",(ttn_current,(centralflux_strength=centralflux_strength,)))])

if false
	time_end = 1.0
	time_change = 0.1
	tilt_stren = 0.00

	rez,swphndler = evolve_in_time(psi_gs,time_end,time_change,ham; if_states=if_states,obs_measures=whatmeasuring)
	append!(wavefuncs,rez.measurement_results["states"])
	append!(alloccs,rez.measurement_results["occs"])
	times = [[0.0]; rez.times]
end

#
fig = figure()
plot(times[2:end],rez.energy,"-p")
xlabel("Time")
ylabel("Energy")

avgpos = average_position(alloccs)
fig2 = figure()
plot(times,[avgpos[i][1] for i in 1:length(alloccs)],"-p",label="Real")
plot(times,[avgpos[i][2] for i in 1:length(alloccs)],"-p",label="Synth")
xlabel("Time")
legend()
ylabel("Average Position")

allcurrents = [[ttn_current(psi_gs; centralflux_strength=centralflux_strength)]; rez.measurement_results["currents"]]
alltotals = [allcurrents[i][1] for i in 1:length(allcurrents)]
allmats = [allcurrents[i][2] for i in 1:length(allcurrents)]
fig3 = figure()
plot(times,real.(alltotals),"-p")
xlabel("Time")
ylabel("Current")
=#




















"fin"

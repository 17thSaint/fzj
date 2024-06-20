#using Pkg
#Pkg.activate(".")
include("../review-practice-codes/ttn.jl")

function long_range_scaling_old(x_final,virt_edge_length,initial_strength; kwargs...)
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
	end
	
	if if_hard_cutoff
		strengths[x_final + 2:end] .= 0.0
	elseif if_rounding
		final_index = findfirst(x -> x .<= trunc,strengths)
		if !isnothing(final_index)
			strengths[final_index:end] .= 0.0
		end
	end
	
	if_plot || if_save_fig ? plot_long_range_scaling_old(strengths,virt_edge_length; kwargs...) : nothing
	#if_save_data ? save_long_range_scaling_old(strengths,virt_edge_length; kwargs...) : nothing

	return strengths
end

function long_range_scaling_ed(x_final::Int64,virt_edge_length::Int64,initial_strength::Float64; kwargs...)

	trunc = get(kwargs, :trunc_digits, 5)
	scaling_func = get(kwargs, :scaling, "flat")
	
	strengths = zeros(virt_edge_length)

    if x_final == 0.0 || initial_strength == 0.0
        strengths[1] = initial_strength != 0.0 ? initial_strength : 1.0
        return strengths
    end
	
	if scaling_func == "flat"
		strengths[1:x_final+1] .= initial_strength
	elseif scaling_func == "exp"
        corr_length = get(kwargs, :corr_length, virt_edge_length)
		strengths = map(1:virt_edge_length) do x
			initial_strength * exp(-(x-1)/corr_length)	
		end
	elseif scaling_func == "rydberg"
		blockade_radius = get(kwargs, :blockade_radius, 1.0)
		strengths = map(0:virt_edge_length-1) do x
			initial_strength * (blockade_radius^6) / (blockade_radius^6 + x^6)
		end
	end
	
	strengths = round.(strengths,digits=trunc)

	return strengths
end

function save_long_range_scaling_old(strengths,virt_edge_length; kwargs...)
	filename = get(kwargs, :name, "scaling-strength")
	location = get(kwargs, :location, pwd())
	xcoord = [i for i in 0:virt_edge_length-1]
	scaling_data = Dict([("strengths",strengths),("xcoord",xcoord)])
	write_data_jld2(filename,scaling_data,location=location; kwargs...)
	return
end

function plot_long_range_scaling_old(strengths,virt_edge_length; kwargs...)
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

function build_HH_net_old(num_layers; kwargs...)
	conserve_qns = get(kwargs, :syms, true)
	if_fermion = get(kwargs, :if_fermion, false)
	particle_type = if_fermion ? "Fermion" : "Boson"
	max_occ = get(kwargs,:max_occ,1)
	
	net = if_fermion ? TTNKit.BinaryRectangularNetwork(num_layers, TTNKit.ITensorNode, particle_type;conserve_nf=conserve_qns,conserve_nfparity=false) : TTNKit.BinaryRectangularNetwork(num_layers, TTNKit.ITensorNode, particle_type;conserve_qns=conserve_qns,dim=max_occ+1)
	
	return net
end

function get_interaction_coords_old(given_site,inter_dist,lat,if_periodic_virt) # written by ChatGPT 12.06.2023 then vastly edited 13.06.2023 by me
	virtual, physical = given_site
	#physical, virtual = given_site
	coordinates = []
    
	virt_edge_length, phys_edge_length = size(lat)
	#if typeof(given_site) == Int64
	#	given_site = TTNKit.coordinate(lat,given_site)
	#end

	for shift in [-inter_dist % virt_edge_length,inter_dist % virt_edge_length]
		new_virtual = virtual + shift
				
		# Apply periodic boundary conditions along the virtual-axis
		#=if if_periodic_virt
			if new_virtual < 1
				new_virtual += virt_edge_length
			elseif new_virtual > virt_edge_length
				new_virtual -= virt_edge_length
			end
		end=#

        # Check if new coordinates are within lattice dimensions
		if 1 <= new_virtual <= virt_edge_length && new_virtual != virtual
 			append!(coordinates, [[new_virtual,physical]])
			#append!(coordinates, [[physical,new_virtual]])
 		#else
 			#println("Still Outside Lattice")
		end

	end
	return unique(coordinates)
end

function get_inter_coeff_old(s1,s2,t_strength,phi,edge_length_x,edge_length_y; kwargs...) 
	hopping_anisotropy = get(kwargs, :hopping_anisotropy, 1.0)
	if_periodic_phys = get(kwargs, :if_periodic_phys, false)
	if_periodic_virt = get(kwargs, :if_periodic_virt, false)
	#t_strength_phys = t_strength * hopping_anisotropy
	if hopping_anisotropy < 1.0
		t_strength_synth = t_strength / hopping_anisotropy
		t_strength_phys = t_strength
	else
		t_strength_phys = t_strength * hopping_anisotropy
		t_strength_synth = t_strength
	end
	if get(kwargs, :no_magF, false)
		phi = 0.0
	end
	if s1[1] == s2[1] # Synthetic Dimension Hopping
		thetay = get(kwargs, :thetay, thetay_2)
		#=if ==(edge_length,s1[2])
			println("Using ThetaY")
		end
		=#
		stren = -t_strength_synth
		if !if_periodic_virt && if_periodic_phys
			stren *= exp(im*2*pi*phi*s1[1])
		end
		return round(stren,digits=10)
	elseif s1[2] == s2[2] # Physical Dimension Hopping
		thetax = get(kwargs, :thetax, thetax_2)
		#=if ==(edge_length,s1[1])
			println("Using ThetaX")
		end
		=#
		stren = -t_strength_phys * exp(im*2*pi*phi*s1[2])
		if !if_periodic_virt && if_periodic_phys
			stren = -t_strength_phys
		end
		return round(stren,digits=10)
	else
		return 0.0
	end

end

function long_range_HH_ham_old(net,t_strength,phi; kwargs...)
	resulting_ham = []
	phys_edge_length,virt_edge_length = get_lattice_dims(net)
	println("Phys = ",phys_edge_length,", Virt = ",virt_edge_length)
	
	scaling_distance = get(kwargs, :lr, 0)
	
	restricted_size = get(kwargs, :restricted_size, [phys_edge_length,virt_edge_length])
	if_periodic_virt = get(kwargs, :if_periodic_synth, false)
	if_periodic_phys = get(kwargs, :if_periodic_phys, false)
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
	
	long_range_strengths = long_range_scaling_ed(scaling_distance,virt_edge_length,onsite_strength; kwargs...)
	display(long_range_strengths)
	if_interaction = !all(long_range_strengths.==0)
	
	lat = TTNKit.physical_lattice(net)
	
	if if_hopping
		if_periodic_phys ? nothing : centralflux_strength = 0.0
		hopping = TTNKit.OpSum()
		#
		for (s1,s2) in TTNKit.nearest_neighbours(lat,collect(1:TTNKit.number_of_sites(lat)))
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)
						
			coeff = get_inter_coeff_old(s1_coord,s2_coord,t_strength,phi,phys_edge_length,virt_edge_length; kwargs...)
			
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
				s1_coord = (i,1)
				s2_coord = (i,restricted_size[2])
				coeff = get_inter_coeff_old(s1_coord,s2_coord,t_strength,phi,phys_edge_length,virt_edge_length; kwargs...)
				hopping += (coeff,"Adag",s1_coord,"A",s2_coord)
				hopping += (conj(coeff),"Adag",s2_coord,"A",s1_coord)
			end
		end

		if if_periodic_phys
			for i in 1:restricted_size[2]
				s1_coord = (1,i)
				s2_coord = (restricted_size[1],i)
				coeff = get_inter_coeff_old(s1_coord,s2_coord,t_strength,phi,phys_edge_length,virt_edge_length; kwargs...)
				coeff *= exp(im*2*pi*centralflux_strength/size(lat)[1])
				coeff *= exp(im*twist_angle*2*pi)
				hopping += (coeff,"Adag",s1_coord,"A",s2_coord)
				hopping += (conj(coeff),"Adag",s2_coord,"A",s1_coord)
			end
		end

		append!(resulting_ham,[hopping])
	end
	
	if if_interaction
		interaction = TTNKit.OpSum()
		for (idx,stren) in enumerate(long_range_strengths)
			if stren == 0.0
				continue
			else
				if idx == 1
					#=for j in TTNKit.eachindex(lat)
						interaction += (stren,"N * N",TTNKit.coordinate(lat,j))
					end=#
					continue
				else
					for j in TTNKit.eachindex(lat)
						s_coord = TTNKit.coordinate(lat,j)
						if s_coord[1] > phys_edge_length || s_coord[2] > virt_edge_length
							continue
						end
						interaction_sites = get_interaction_coords_old(s_coord,idx-1,lat,if_periodic_virt)
						
						for k in interaction_sites
							if k[1] > restricted_size[1] || k[2] > restricted_size[2]
								continue
							end
							interaction += (stren/2,"Adag * A",s_coord,"Adag * A",Tuple(k))
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














"fin"
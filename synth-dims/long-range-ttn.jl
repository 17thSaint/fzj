using PyPlot
include("../review-practice-codes/ttn.jl")

function long_range_scaling(x_final,virt_edge_length,initial_strength; kwargs...)
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

	return strengths
end

function build_HH_net(num_layers; kwargs...)
	conserve_qns = get(kwargs, :syms, true)
	max_occ = 1
	
	net = TTNKit.BinaryRectangularNetwork(num_layers, TTNKit.ITensorNode, "Boson";conserve_qns=conserve_qns,dim=max_occ+1)
	
	return net
end

function get_interaction_coords(given_site,inter_dist,lat; kwargs...) # written by ChatGPT 12.06.2023
	y, x = given_site
	coordinates = []
    
	phys_edge_length,virt_edge_length = get(kwargs, :edges, sort(size(lat)))
	if typeof(given_site) == Int64
		given_site = TTNKit.coordinate(lat,given_site)
	end

	for i = -inter_dist:inter_dist
		for j = -inter_dist:inter_dist
			if i^2 + j^2 == inter_dist^2
				new_x = x + i
				new_y = y + j
				
				# Apply periodic boundary conditions along the x-axis
				if new_x < 1
					new_x += virt_edge_length
				elseif new_x > virt_edge_length
					new_x -= virt_edge_length
				end

                		# Check if new coordinates are within lattice dimensions
				if 1 <= new_x <= virt_edge_length && 1 <= new_y <= phys_edge_length && new_x != x
                    			append!(coordinates, [[new_y,new_x]])
                		end
            		end
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
						interaction_sites = get_interaction_coords(TTNKit.coordinate(lat,j),i-1,lat; edges=(virt_edge_length,phys_edge_length))
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

#
if_per = false
mag_off = true
evolve = true
chemical = false
mu = 0.5
#max_occupation = 3
expan = TTNKit.DefaultExpander(0.5)
ts = 0.01
nu = 1/2
layers = 6
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
mdim = 100
nswps = 3
println("Using $num_particles particles on $tot_sites sites")

if_cliff = true

net = build_HH_net(layers; syms=true)
#=
all_ttns = []
for i in 0:2
	longrange_dist = i
	title_string = "LR = $longrange_dist"

	ham = long_range_HH_ham(net,ts,alpha; scaling="flat",scaling_dist=longrange_dist,cliff=if_cliff,if_periodic=if_per,if_chem=chemical,no_magF=mag_off)

	og_ttn, hamilt, dm_sp = build_full_harperhofstadter(layers,num_particles,ts,nu; ttn_net=net,ham_op=ham,max_dim=mdim, num_sweeps=nswps,phi=alpha, if_periodic=if_per,max_occ=1,if_sweep=evolve,sweep_type="dmrg",expander=expan,if_chem=chemical,chem_strength=mu,no_magF=mag_off,output_level=0)
	append!(all_ttns,[dm_sp.ttn])
	#
	rez1 = get_occupancy(dm_sp.ttn; plot_title=title_string)
	#rez2 = get_current_yfunc(dm_sp.ttn)
	#rez3_fqh = get_ydir_greenfunc(dm_sp.ttn)
	rez3_sf = get_ydir_greenfunc(dm_sp.ttn; plot_title=title_string)
	rez4 = get_xdir_greenfunc(dm_sp.ttn; plot_title=title_string)
	#
end
=#




























"fin"

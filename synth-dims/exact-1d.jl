using LinearAlgebra,PyPlot,Random

function log_prod(all_values)
	#=
	log_version = []
	all_neg_values = []
	for i in 1:length(all_values)
		if all_values[i] > 0.0
			append!(log_version,log(all_values[i]))
		else
			append!(all_neg_values,all_values[i])
		end
	end
	if length(all_neg_values) == 0
		return sum(log_version)
	elseif iseven(length(all_neg_values))
		append!(log_version,[log(prod(all_neg_values))])
		return sum(log_version)
	else 
		#println("Negative Wavefunction")
		return nothing
	end
	=#
	log_version = log.(all_values)
	return sum(log_version)
end

function log_add(a,b)
	if real(a) > real(b)
		ordered::Vector{typeof(a)} = [b,a]
	else
		ordered = [a,b]
	end
	result::ComplexF64 = ordered[2] + log(Complex(1 + exp(ordered[1] - ordered[2])))
	return Complex(result)
end

function log_sum(all_values)
	consecutive = [all_values[1],all_values[2]]
	for i in 3:length(all_values) + 1
		added_value = log_add(consecutive[1],consecutive[2])
		consecutive[1] = added_value
		consecutive[2] = i <= length(all_values) ? all_values[i] : 0.0
	end
	return consecutive[1]
end

function physical_part(particle_dictionary::Dict,L::Int; kwargs...)
	if_log = get(kwargs, :if_log, true)
	if_periodic = get(kwargs, :if_periodic, false)
	
	bb_part = []
	ff_part = []
	bf_part = []
	
	num_bosons,num_fermions = length(particle_dictionary["B"]),length(particle_dictionary["F"])
	
	for (species,locations) in particle_dictionary
		if species == "B"
			for i in 1:length(locations)
				for j in 1:i-1
					interaction = if_periodic ? sin(pi*(locations[j] - locations[i])/L) : abs(locations[j] - locations[i])
					append!(bb_part,[interaction])
				end
			end
		else
			for i in 1:length(locations)
				for j in 1:i-1
					interaction = if_periodic ? sin(pi*(locations[j] - locations[i])/L) : locations[j] - locations[i]
					append!(ff_part,[interaction])
				end
			end
		end
	end
	
	if num_bosons != 0 && num_fermions != 0
		for i in 1:length(particle_dictionary["B"])
			for j in 1:length(particle_dictionary["F"])
				interaction = if_periodic ? sin(pi*(particle_dictionary["B"][i] - particle_dictionary["F"][j])/L) : abs(particle_dictionary["B"][i] - particle_dictionary["F"][j])
				append!(bf_part,[interaction])
			end
		end
	end
	
	included_parts = []
	num_bosons != 0 ? append!(included_parts,bb_part) : nothing
	num_fermions != 0 ? append!(included_parts,ff_part) : nothing
	num_bosons != 0 && num_fermions != 0 ? append!(included_parts,bf_part) : nothing
	
	full_physical_wavefunc = if_log ? log_prod(Complex.(included_parts)) : prod(included_parts)
	
	return full_physical_wavefunc,bb_part,ff_part,bf_part
end

function orbital_part(particle_dictionary::Dict,L::Int; kwargs...)
	if_log = get(kwargs, :if_log, true)
	freq = get(kwargs, :freq, 1.0)
	mass = get(kwargs, :mass, 1.0)
	if_periodic = get(kwargs, :if_periodic, false)
	if if_periodic
		if if_log
			return log(1.0)
		else
			return 1.0
		end
	end
	
	num_bosons,num_fermions = length(particle_dictionary["B"]),length(particle_dictionary["F"])
	all_values = []
	xosc = sqrt(1/(mass*freq))
	hermite_part = 1.0
	for (species,locations) in particle_dictionary
		for loc in locations
			coeff_part = exp(-0.5*((loc-L/2)/xosc)^2) / sqrt((2^0)*factorial(0)*sqrt(pi)*xosc)
			part = coeff_part * hermite_part
			append!(all_values,[part])
		end
	end
	
	return if_log ? log_prod(Complex.(all_values)) : prod(all_values)	
end

function exchange_particles(species::String,particle_dictionary::Dict)
	if species == "M"
		part_locsB = particle_dictionary["B"]
		part_locsF = particle_dictionary["F"]
		num1, num2 = rand(1:length(part_locsB)),rand(1:length(part_locsF))

	    	num1_ogloc = particle_dictionary["B"][num1]
	    	num2_ogloc = particle_dictionary["F"][num2]
	    	exchanged_dict = particle_dictionary
	    	exchanged_dict["B"][num1] = num2_ogloc
    		exchanged_dict["F"][num2] = num1_ogloc
	else
		part_locs = particle_dictionary[species]
		num1, num2 = rand(1:length(part_locs), 2)
	    	while num1 == num2
			num2 = rand(1:length(part_locs))
	    	end
	    	num1_ogloc = particle_dictionary[species][num1]
	    	num2_ogloc = particle_dictionary[species][num2]
	    	exchanged_dict = particle_dictionary
	    	exchanged_dict[species][num1] = num2_ogloc
    		exchanged_dict[species][num2] = num1_ogloc
    	end

    	
    	return exchanged_dict
end

function normalize_wavefunc(wavefunc; kwargs...)
	if_log = get(kwargs, :if_log, true)
	return if_log ? wavefunc / 2 - conj(wavefunc) / 2 : nothing
end

function get_config_wavefunc(particle_dictionary::Dict,L::Int; kwargs...)
	if_log = get(kwargs, :if_log, true)
	freq = get(kwargs, :freq, 1.0)
	mass = get(kwargs, :mass, 1.0)
	if_periodic = get(kwargs, :if_periodic, false)
	
	phys = physical_part(particle_dictionary,L; kwargs...)[1]
	orb = orbital_part(particle_dictionary,L; kwargs...)
	if isnothing(phys) | isnothing(orb)
		return nothing
	else
		return if_log ? phys + orb : phys * orb
	end
end

function config_vector_to_string(config)
	return join(config,",")	
end

function config_string_to_vector(config)
	return parse.(Int64,split(config,","))
end

function make_config_dictionary(configB::Vector,configF::Vector)
	return Dict([("B",configB),("F",configF)])
end

function normalize_wavefunc_dict(wavefunc_dict::Dict)
	all_vals = collect(values(wavefunc_dict))
	log_psisquared = log_sum(all_vals .+ conj.(all_vals))
	norm_factor = -log_psisquared/2
	normed_wavefunc_dict = Dict(k => v + norm_factor for (k,v) in wavefunc_dict)
	return normed_wavefunc_dict
end

function get_wavefunc(configurations,L::Int; kwargs...)
	config_wavefuncs = Dict()
	for con in configurations
		config_dict = make_config_dictionary(con,[])
		config_wavefuncs[config_vector_to_string(con)] = get_config_wavefunc(config_dict,L; kwargs...)
	end
	return normalize_wavefunc_dict(config_wavefuncs)
end

function assign_locations(n_F::Int, n_B::Int, L::Int) # written by ChatGPT 24.07.2023
    if n_F + n_B > L
        throw(ArgumentError("Total number of particles exceeds the available locations (L)."))
    end

    # Create a shuffled array of all possible locations
    all_locations = randperm(L)

    # Assign locations for species F and B
    locations_dict = Dict("F" => all_locations[1:n_F], "B" => all_locations[n_F+1:n_F+n_B])

    return locations_dict
end

function configurations(N::Int, L::Int)
    # Initialize an array to store the current configuration
    current_config = zeros(Int, N)

    # Initialize an array to store all configurations
    all_configs = []

    # Recursive function to generate configurations
    function generate_configs(start_pos::Int, remaining_particles::Int)
        # Base case: all particles have been placed
        if remaining_particles == 0
            push!(all_configs, copy(current_config))
            return
        end

        # Recursive step: try all possible positions for the next particle
        for pos in start_pos:L
            # Check if the position is already occupied
            if pos in current_config[N - remaining_particles + 1:end]
                continue  # Skip this position if it's already occupied
            end

            current_config[N - remaining_particles + 1] = pos
            generate_configs(pos + 1, remaining_particles - 1)
        end
    end

    # Start the recursive function
    generate_configs(1, N)

    return all_configs
end

function cr_anh_pair_configs(all_configurations::Vector,cr_site::Int,anh_site::Int)
	kept_configs = []
	for c in all_configurations
		if cr_site in c && !(anh_site in c)
			append!(kept_configs,[c])
		end
	end
	return kept_configs
end

function momentum_dist_1d(position_dist,p_count,p_end,p_start=0.0)
	num_sites = length(position_dist)
	mom_occ = [0.0*im for i in 1:p_count]
	momenta = [p_start + (i-1)*(p_end - p_start)/(p_count-1) for i in 1:p_count]
	for i in 1:p_count
		momentum = momenta[i]
		exp_vec = [exp(im*momentum*(j - num_sites/2) + position_dist[j]) for j in 1:num_sites]
		local_value = sum(exp_vec)
		mom_occ[i] = local_value
	end
	return momenta,mom_occ
end

function normalize_log_occ_vector(full_vector,particle_count)
	log_of_norm_value = log(particle_count) - log_sum(full_vector)
	normed_log_vector = full_vector .+ log_of_norm_value
	return normed_log_vector
end

#=
pos_occs = []
mom_occs = []

#for L in [10,20,30,40,50]
L = 10
#println(L)
n_tot = 5
n_F = 0
n_B = n_tot - n_F
if_per = false

if true
locB_probs = [0.0 for i in 1:L]
#locF_probs = [0.0 for i in 1:L]
sites = [i-L/2 for i in 1:L]
nothing_counts = 0
all_configs = configurations(n_tot,L)
for i in 1:length(all_configs)
	local_config = all_configs[i]
	pd = Dict([("B",local_config),("F",[])])
	#pd = assign_locations(n_F,n_B,L)
	wavefunc = get_wavefunc(pd,L; if_periodic=if_per)
	loc_prob = wavefunc
	if i != 1
		#=
		for j in 1:n_F
			locF_probs[pd["F"][j]] = log_add(locF_probs[pd["F"][j]],loc_prob + conj(loc_prob))
		end
		=#
		for j in 1:n_B
			locB_probs[pd["B"][j]] = log_add(locB_probs[pd["B"][j]],loc_prob + conj(loc_prob))
		end
	else
		#=
		for j in 1:n_F
			locF_probs[pd["F"][j]] = loc_prob + conj(loc_prob)
		end
		=#
		for j in 1:n_B
			locB_probs[pd["B"][j]] = loc_prob + conj(loc_prob)
		end
	end
end
norm_locB = normalize_log_occ_vector(locB_probs,n_B)
plot(sites,exp.(norm_locB),"-p",label="B")
append!(pos_occs,[norm_locB])
#norm_locF = normalize_log_occ_vector(locF_probs,n_F)
#=
fig = figure()

plot(sites,exp.(norm_locF),"-p",label="F")
title("Bosons = $n_B,Fermions = $n_F")
legend()
=#
end
#=
mcount = 200
mfinal = 10
momB = momentum_dist_1d(norm_locB,mcount,mfinal)
append!(mom_occs,[momB[2]])
#momF = momentum_dist_1d(norm_locF,mcount,mfinal,-2)
#fig2 = figure()
plot(momB[1]./(pi*1),abs.(momB[2])./n_B,label="$L")
#plot(momF[1]./(pi*1),abs.(momF[2])./n_F,"-p",label="F")
end
legend()
=#


=#

































"fin"

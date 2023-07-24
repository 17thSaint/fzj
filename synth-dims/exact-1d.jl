using LinearAlgebra,PyPlot,Random

function log_prod(all_values)
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
	
	full_physical_wavefunc = if_log ? log_prod(included_parts) : prod(included_parts)
	
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
	
	return if_log ? log_prod(all_values) : prod(all_values)	
end

function get_wavefunc(particle_dictionary::Dict,L::Int; kwargs...)
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

L = 20
n_tot = 10
n_F = 5
n_B = n_tot - n_F
samples = 100000
locB_probs = [0.0 for i in 1:L]
locF_probs = [0.0 for i in 1:L]
sites = [i-L/2 for i in 1:L]

nothing_counts = 0
for i in 1:samples
	
	pd = assign_locations(n_F,n_B,L)
	wavefunc = get_wavefunc(pd,L; if_periodic=true)
	loc_prob = wavefunc
	if isnothing(loc_prob)
		global nothing_counts += 1
	else
		for j in 1:n_F
			locF_probs[pd["F"][j]] += 1/loc_prob
		end
		for j in 1:n_B
			locB_probs[pd["B"][j]] += 1/loc_prob
		end
	end
end
println("Nothing Counts = $nothing_counts")
plot(sites,locB_probs./locB_probs[1],"-p",label="B")
plot(sites,locF_probs./locF_probs[1],"-p",label="F")
legend()










































"fin"

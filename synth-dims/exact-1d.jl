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
	if length(all_values) < 2
		return all_values[1]
	else
		consecutive = [all_values[1],all_values[2]]
		for i in 3:length(all_values) + 1
			added_value = log_add(consecutive[1],consecutive[2])
			consecutive[1] = added_value
			consecutive[2] = i <= length(all_values) ? all_values[i] : 0.0
		end
		return consecutive[1]
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
			coeff_part = exp(-0.5*((loc-(L+1)/2)/xosc)^2) / sqrt((2^0)*factorial(0)*sqrt(pi)*xosc)
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
	if length(configurations) == 0
		return nothing
	else
		config_wavefuncs = Dict()
		for con in configurations
			config_dict = make_config_dictionary(con,[])
			config_wavefuncs[config_vector_to_string(con)] = get_config_wavefunc(config_dict,L; kwargs...)
		end
		return normalize_wavefunc_dict(config_wavefuncs)
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

function combine_two_vectors(vec1::Vector,vec2::Vector; kwargs...)
	if_log = get(kwargs, :if_log, true)
	
	halfway_vec = [0.0*im for j in 1:length(vec1)]
	for i in 1:length(vec1)
		halfway_vec[i] = if_log ? log_sum(vec2 .+ vec1[i]) : sum(vec2 .* vec1[i])
	end
	result = if_log ? log_sum(halfway_vec) : sum(halfway_vec)
	return result
end

function normalize_densmat(dens_mat::Matrix,part_count::Int; kwargs...)
	if_log = get(kwargs, :if_log, true)
	current_trace = if_log ? log_sum(diag(dens_mat)) : tr(dens_mat)
	println(current_trace)
	shift_mat = Diagonal([log(part_count) - current_trace for i in 1:L])
	norm_densmat = dens_mat + shift_mat
	return norm_densmat
end

function density_matrix(wavefunc,L::Int,part_count::Int; kwargs...)
	if_log = get(kwargs, :if_log, true)
	
	dens_mat = zeros(L,L) .* im
	for x in 1:L
		psi0s = []
		for (k,v) in wavefunc
			if string(x) in collect(split(k,","))
				append!(psi0s,[v])
			end
		end
		for xp in 1:L
			psi0ps = []
			for (k,v) in wavefunc
				if string(xp) in collect(split(k,","))
					append!(psi0ps,[conj(v)])
				end
			end
			dens_mat[x,xp] = combine_two_vectors(psi0s,psi0ps; kwargs...)
		end
	end
	return normalize_densmat(dens_mat,part_count; kwargs...)
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

function overlap_two_wavefuncs(wavefuncL,wavefuncR; kwargs...)
	if_log = get(kwargs, :if_log, true)
	if isnothing(wavefuncL) | isnothing(wavefuncR)
		return -Inf
	else
	
		overlaps = []
		for (k,v) in wavefuncL
			if k in keys(wavefuncR)
				config_overlap = if_log ? conj(wavefuncL[k]) + wavefuncR[k] : conj(wavefuncL[k]) + wavefuncR[k]
				append!(overlaps,[config_overlap])
			end
		end
		if if_log
			base_val = log_sum(overlaps)
			final_overlap = base_val + conj(base_val)
		else
			final_overlap = abs2.(sum(overlaps))
		end
		
		return final_overlap
	end
end

function position_occupancy(wavefunc::Dict,L::Int; kwargs...)
	if_log = get(kwargs, :if_log, true)
	if_plot = get(kwargs, :if_plot, true)
	all_sites = [i for i in 1:L]
	site_occs = [423.0 for i in 1:L]
	for (k,v) in wavefunc
		result = if_log ? v + conj(v) : v*conj(v)
		config = config_string_to_vector(k)
		for s in config
			if site_occs[s] == 423.0
				site_occs[s] = result
			else
				site_occs[s] = if_log ? log_add(site_occs[s],result) : site_occs[s] + result
			end
		end
	end
	
	if_plot ? plot_position_occupancy(all_sites,if_log ? exp.(site_occs) : site_occs; kwargs...) : nothing
	
	return all_sites,site_occs
end

function plot_position_occupancy(sites,occs; kwargs...)
	plot_label = get(kwargs, :plot_label, "")
	isempty(plot_label) ? fig = figure() : nothing
	title_string = "Position Occupancy, " * get(kwargs, :plot_title, "")
	plot(sites,occs,"-p",label=plot_label)
	title(title_string)
end

function get_CrAnh_correlation(wavefunc,all_configs::Vector,L::Int; kwargs...)
	all_hops = zeros(L,L)
	for i in 1:L
		for j in 1:L
			#println(i/L,", ",j/L)
			new_configs = cr_anh_pair_configs(all_configs,i,j)
			hop_wavefunc = get_wavefunc(new_configs,L;kwargs...)
			overlap_val = overlap_two_wavefuncs(wavefunc,hop_wavefunc; kwargs...)
			all_hops[i,j] = exp.(overlap_val)
		end
	end
	return all_hops
end

function momentum_dist_1d(wavefunc::Dict,part_count::Int,L::Int,p_count::Int,p_end::Float64,all_configs::Vector,p_start=0.0;kwargs...)
	if_plot = get(kwargs, :if_plot, true)
	num_sites = L
	mom_occ = [0.0*im for i in 1:p_count]
	momenta = [p_start + (i-1)*(p_end - p_start)/(p_count-1) for i in 1:p_count]
	pos_occ = get_CrAnh_correlation(wavefunc,all_configs,L; kwargs...)
	println("Made Position Occupation Matrix")
	for i in 1:p_count
		#println(round(100*i/p_count,digits=1),"%")
		momentum = momenta[i]
		exp_vect = zeros(num_sites,num_sites) .* im
		for j in 1:num_sites
			exp_vect[:,j] = [exp(im*momentum*(j-l)) for l in 1:num_sites]
		end
		
		mom_occ[i] = abs(sum(exp_vect .* pos_occ) / (num_sites))
	end
	
	if_plot ? plot_momentum(momenta,mom_occ,part_count; kwargs...) : nothing
	
	return momenta,mom_occ
end

function plot_momentum(momenta,mom_occ,part_count::Int; kwargs...)
	freq = get(kwargs, :freq, 1.0)
	mass = get(kwargs, :mass, 1.0)
	kosc = 2*pi/sqrt(1/(mass*freq))
	plot_label = get(kwargs, :plot_label, "")
	isempty(plot_label) ? fig = figure() : nothing
	title_string = "Momentum Distribution, " * get(kwargs, :plot_title, "")
	plot(momenta./pi,mom_occ./part_count,label=plot_label)
	if_logscale = get(kwargs, :if_logscale, true)
	if if_logscale
		yscale("log")
		xscale("log")
	end
	xlabel("p/pi")
	ylabel("Occupation / nparticles")
	title(title_string)
	
end

function normalize_log_occ_vector(full_vector,particle_count)
	log_of_norm_value = log(particle_count) - log_sum(full_vector)
	normed_log_vector = full_vector .+ log_of_norm_value
	return normed_log_vector
end

#
#for L in [10,20,30,40,50]
L = 10
n_tot = 5
n_F = 0
n_B = n_tot - n_F
if_per = false
pcount = 100
all_configs = configurations(n_tot,L)
println("Made Configs, total = ",length(all_configs))

#println(L)
moms = []
count = 10
fs = 0.01
fe = 2.0
omega = 1.0
#for omega in [fs + (i-1)*(fe-fs)/(count-1) for i in 1:count]
model_paras = (freq=omega,if_periodic=if_per,if_log=true,if_logscale=false)
gs_wavefunc = get_wavefunc(all_configs,L;model_paras...)
println("Made Wavefunc")
rho = density_matrix(gs_wavefunc,L,n_tot; model_paras...)
imshow(real.(exp.(rho)))
colorbar()
#position_occupancy(gs_wavefunc,L; model_paras...,plot_label="$(round(omega,digits=2))",plot_title="range HarmTrap Frequency, Nbosons=$n_tot")
#mrez = momentum_dist_1d(gs_wavefunc,n_tot,L,pcount,10.0,all_configs;model_paras...,plot_label="$(round(omega,digits=2))",plot_title=" range HarmTrap Frequency, Nbosons=$n_tot")
#append!(moms,[mrez[2]])
#end
#legend()



#

































"fin"

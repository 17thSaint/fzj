using LinearAlgebra,PyPlot,Random,NumericalIntegration

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

function remove_neg_infs(all_values)
	to_keep = []
	for i in all_values
		if !isinf(i)
			append!(to_keep,[i])
		end
	end
	return to_keep
end

function log_sum(all_values)
	all_values = remove_neg_infs(all_values)
	if 1 < length(all_values) < 2
		return all_values[1]
	elseif length(all_values) < 1
		return -Inf
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
	#hermite_part = 1.0
	for (species,locations) in particle_dictionary
		for loc in locations
			#coeff_part = exp(-0.5*((loc-(L+1)/2)/xosc)^2) / sqrt((2^0)*factorial(0)*sqrt(pi)*xosc)
			#part = hermite_part * coeff_part
			part = if_log ? -0.5 * ((loc-(L+1)/2)/xosc)^2 : exp(-0.5 * ((loc-(L+1)/2)/xosc)^2)
			append!(all_values,[part])
		end
	end
	
	return if_log ? sum(Complex.(all_values)) : prod(all_values)	
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

function stirling_formula(n::Int)
	if n == 0
		return 0.0
	end
	return n*log(n) - n
end

function get_coeff(num_parts::Int; kwargs...)
	if_log = get(kwargs, :if_log, true)
	freq = get(kwargs, :freq, 1.0)
	mass = get(kwargs, :mass, 1.0)
	
	xosc = sqrt(1/(mass*freq))
	if if_log
		p1 = 0.25*num_parts*(num_parts-1)*log(2)
		p2 = -0.5*num_parts*log(xosc)
		p3 = -0.5 * (0.5*num_parts*log(pi) + sum([log(factorial(i)) for i in 0:num_parts]))
		return p1+p2+p3
	else
		p1 = 2^(0.25*num_parts*(num_parts-1))
		p2 = xosc^(-0.5*num_parts)
		p3 = (1/sqrt(factorial(num_parts))) * (1/sqrt(prod([factorial(i)*sqrt(pi) for i in 0:num_parts-1])))
		return p1*p2*p3
	end
end

function get_config_wavefunc(particle_dictionary::Dict,L::Int; kwargs...)
	if_log = get(kwargs, :if_log, true)
	freq = get(kwargs, :freq, 1.0)
	mass = get(kwargs, :mass, 1.0)
	if_periodic = get(kwargs, :if_periodic, false)
	
	num_particles = length(particle_dictionary["B"]) + length(particle_dictionary["F"])
	phys = physical_part(particle_dictionary,L; kwargs...)[1]
	orb = orbital_part(particle_dictionary,L; kwargs...)
	coeff = get_coeff(num_particles; kwargs...)
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

function get_random_configurations(N::Int,L::Int,limit::Int)
	all_configs = []
	repeat_count = 0
	while length(all_configs) < limit && repeat_count < 50
		next_config = randperm(L)[1:N]
		if !(next_config in all_configs)
			append!(all_configs,[next_config])
		else
			repeat_count += 1
		end
		
	end
	return all_configs
end

function configurations(N::Int, L::Int; kwargs...)
    # Calculate number of configurations
    if L < 30
	    total_length = factorial(L)/factorial(L-N)
    else
    	    total_length = 100000.0
    end
    
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
    if_limit = get(kwargs, :if_limit, true)
    limit = get(kwargs, :limit, 5000)
    if total_length < 2*limit
    	generate_configs(1, N)
    	if length(all_configs) > limit
    		println("More than limit but accessible")
	    	which_keep = randperm(length(all_configs))[1:limit]
	    	kept_configs = [all_configs[i] for i in which_keep]
	    	return kept_configs
	else
		println("Less than limit")
		return all_configs
	end
    else
    	println("Beyond limit $total_length, random generation")
    	if if_limit
    		all_configs = get_random_configurations(N,L,limit)
    		return all_configs
    	else
    		generate_configs(1, N)
    		return all_configs
    	end
    end
    
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
	L = size(dens_mat)[1]
	current_trace = if_log ? log_sum(diag(dens_mat)) : tr(dens_mat)
	shift_mat = Diagonal([log(part_count) - current_trace for i in 1:L])
	norm_densmat = dens_mat + shift_mat
	return norm_densmat
end

function wrong_density_matrix(wavefunc,L::Int,part_count::Int; kwargs...)
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
	if cr_site == anh_site
		same_site = true
	else
		same_site = false
	end
	kept_configs = []
	for c in all_configurations
		if cr_site in c && (same_site ? true : !(anh_site in c))
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

function plot_position_occupancy(sites::Vector,occs::Vector; kwargs...)
	plot_label = get(kwargs, :plot_label, "")
	isempty(plot_label) ? fig = figure() : nothing
	title_string = "Position Occupancy, " * get(kwargs, :plot_title, "")
	plot(sites,occs,"-p",label=plot_label)
	title(title_string)
end

function plot_position_occupancy(dens_mat::Matrix; kwargs...)
	if_log = get(kwargs, :if_log, true)
	plot_label = get(kwargs, :plot_label, "")
	isempty(plot_label) ? fig = figure() : nothing
	title_string = "Position Occupancy, " * get(kwargs, :plot_title, "")
	sites = [i for i in 1:size(dens_mat)[1]]
	occs = if_log ? diag(exp.(dens_mat)) : diag(dens_mat) 
	plot(sites,occs,"-p",label=plot_label)
	title(title_string)
end

function pair_dist(x1,x2,L::Int; kwargs...)
	if_log = get(kwargs, :if_log, true)
	sp_dens_1 = orbital_part(Dict([("B",[x1]),("F",[])]),L; kwargs...)
	sp_dens_2 = orbital_part(Dict([("B",[x2]),("F",[])]),L; kwargs...)
end

function density_matrix(wavefunc,all_configs::Vector,L::Int; kwargs...)
	if_log = get(kwargs, :if_log, true)
	all_hops = zeros(L,L)
	for i in 1:L
		for j in 1:L
			#println(i,", ",j)
			new_configs = cr_anh_pair_configs(all_configs,i,j)
			hop_wavefunc = get_wavefunc(new_configs,L;kwargs...)
			overlap_val = overlap_two_wavefuncs(wavefunc,hop_wavefunc; kwargs...)
			all_hops[i,j] = overlap_val
		end
	end
	full_trace = if_log ? tr(exp.(all_hops)) : tr(all_hops)
	println("Trace of DensMat is $(round(full_trace,digits=0))")
	return all_hops
end

function get_dm_coeff(part_count::Int; kwargs...)
	freq = get(kwargs, :freq, 1.0)
	mass = get(kwargs, :mass, 1.0)
	
	xosc = 1/sqrt(mass*freq)
	p1 = log(part_count) + 0.5*part_count*(part_count-1)*log(2) - log(xosc) - 0.5*part_count*log(pi)
	p2 = -sum([log(factorial(i)) for i in 0:part_count])
	return p1+p2
end

function restrict_configs(all_configs::Vector,sites_not_allowed::Vector)
	remaining_configs = []
	for c in all_configs
		boolean_check = [!(k in c) for k in sites_not_allowed]
		if all(boolean_check)
			append!(remaining_configs,[c])
		end
	end
	return remaining_configs
end

function other_particles_part(config::Vector)
	all_vals = []
	for j in 1:length(config)
		for k in 1:j
			append!(all_vals,[log((config[j] - config[k])^2)])
		end
	end
	return sum(all_vals)
end

function matrix_elem_part(config::Vector,x::Float64,xp::Float64)
	all_vals = []
	for i in 1:length(config)
		local_val = -config[i] + log(abs(config[i] - x)) + log(abs(config[i] - xp))
		append!(all_vals,[local_val])
	end
	return sum(all_vals)
end

# this function will be done entirely in log form
function direct_density_matrix(all_configs::Vector,L::Int; kwargs...)
	freq = get(kwargs, :freq, 1.0)
	mass = get(kwargs, :mass, 1.0)
	
	xosc = 1/sqrt(mass*freq)
	coeff = get_dm_coeff(length(all_configs[1]); kwargs...)
	front_exps = zeros(L,L)
	integral_part = zeros(L,L)
	for x in 1:L
		for xp in 1:x
			println(x,", ",xp)
			all_vals = []
			front_exps[x,xp] = -((x/xosc)^2 + (xp/xosc)^2)/2
			front_exps[xp,x] = front_exps[x,xp]
			restricted_configs = all_configs ./ xosc#restrict_configs(all_configs,[x,xp]) ./ xosc
			for c in restricted_configs
				println(c)
				other_parts_part = other_particles_part(c)
				mat_elem_part = matrix_elem_part(c,x/xosc,xp/xosc)
				full_part = other_parts_part + mat_elem_part
				append!(all_vals,[full_part])
			end
			integral_part[x,xp] = log_sum(all_vals)
			integral_part[xp,x] = integral_part[x,xp]
		end
	end
	return coeff .+ (front_exps + integral_part)
end

function momentum_dist_1d(dens_mat::Matrix,part_count::Int,p_count::Int,p_end::Float64,p_start=0.0;kwargs...)
	if_plot = get(kwargs, :if_plot, true)
	if_norm_mom = get(kwargs, :if_norm_mom, true)
	num_sites = size(dens_mat)[1]
	mom_occ = [0.0 for i in 1:p_count]
	momenta = [p_start + (i-1)*(p_end - p_start)/(p_count-1) for i in 1:p_count]
	for i in 1:p_count
		#println(round(100*i/p_count,digits=1),"%")
		momentum = momenta[i]
		exp_vect = zeros(num_sites,num_sites) .* im
		for j in 1:num_sites
			exp_vect[:,j] = [exp(im*momentum*(l-j)) for l in 1:num_sites]
		end
		
		mom_occ[i] = abs(sum(exp_vect .* exp.(dens_mat)) / num_sites)
	end
	#=
	if if_norm_mom && momenta[end] > pi
		fin_elem = findfirst(x -> momenta[x] < pi && momenta[x+1] > pi,[i for i in 1:length(momenta)])
		norm_factor = 2*integrate(momenta[1:fin_elem+1],mom_occ[1:fin_elem+1])
		mom_occ ./= norm_factor
	end
	=#
	
	if_plot ? plot_momentum(momenta,mom_occ,part_count; kwargs...) : nothing
	
	return momenta,mom_occ
end

function momentum_dist_1d(wavefunc::Dict,part_count::Int,L::Int,p_count::Int,p_end::Float64,all_configs::Vector,p_start=0.0;kwargs...)
	if_log = get(kwargs, :if_log, true)
	if_plot = get(kwargs, :if_plot, true)
	num_sites = L
	mom_occ = [0.0 for i in 1:p_count]
	momenta = [p_start + (i-1)*(p_end - p_start)/(p_count-1) for i in 1:p_count]
	pos_occ = density_matrix(wavefunc,all_configs,L; kwargs...)
	if_log ? pos_occ = exp.(pos_occ) : nothing
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
	plot(momenta./pi,real.(mom_occ)./part_count,label=plot_label)
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

#=
if true
	include("../other-funcs/data-storage-funcs.jl")
	#effective_datadict = read_data_jld2("allmomdist-alpha-0.0-if_periodic-false-if_nn_int-false.jld2","../cluster-data/")
end

if false
params = []
for k in keys(effective_datadict[1])
	L,nflavors,nbosons = parse.(Int,split(k,","))
	append!(params,[[L,nbosons]])
end
params = unique(params)
#=dontkeep = []
for i in 1:length(params)
	if params[i][1] > 20
		append!(dontkeep,[i])
	end
end
deleteat!(params,dontkeep)
=#
display(params)
end

if true
suff_momdiff = 2*pi/200
#L = 20
if_per = false
pfinal = 10.0
pinit = 0.0
pcount = Int(ceil((pfinal-pinit)/suff_momdiff))
println("Running $pcount Mom Vals")
#=

#
count = 10
fs = 1/1000
fe = 1/100
omegas = [fs + (i-1)*(fe-fs)/(count-1) for i in 1:count]
=#
L = 30
#wavefuncs_dict = Dict()
#moms_dict = Dict()
#rhos_dict = Dict()
rhos = []
ns = [2,10,16,20]
omega = 1/10000.0#omegas[i]
model_paras = (freq=omega,if_periodic=if_per,if_log=true,if_logscale=false)
#densities = [1/i for i in 1:5]
#for i in 1:length(omegas)
for n_tot in ns
#keystring = join([string(L),string(n_tot)],",")
#println(keystring)
#n_tot = Int(round(dens*L,digits=0))
n_F = 0
n_B = n_tot - n_F
all_configs = configurations(n_tot,L)
println("Made Configs, total = ",length(all_configs))
gs_wavefunc = get_wavefunc(all_configs,L;model_paras...)
#wavefuncs_dict[keystring] = gs_wavefunc
println("Made Wavefunc")
rho = density_matrix(gs_wavefunc,all_configs,L; model_paras...)
append!(rhos,[rho])
#rhos_dict[keystring] = rho
#fig = figure()
fig = figure()
imshow(real.(rho))
colorbar()
#plottitle = "Freq = $(round(omega,digits=2))"
plottitle = "Dens Mat for Phys Dim = $L and Nbosons = $n_tot"
title(plottitle)
#plot_position_occupancy(rho; model_paras...,plot_label=plottitle)
#append!(rhos,[rho])
#position_occupancy(gs_wavefunc,L; model_paras...,plot_label="$n_tot")
#legend()
#mrez = momentum_dist_1d(gs_wavefunc,n_tot,L,pcount,pfinal,all_configs,pinit;model_paras...,plot_label="$(round(omega,digits=2))",plot_title=" Nbosons=$n_tot")
#mrez = momentum_dist_1d(rho,n_tot,pcount,pfinal,pinit;model_paras...,if_plot=true)
#moms_dict[keystring] = mrez
#end
#append!(moms,[mrez[2]])
#end
end
end

if false
for ke in keys(rhos_dict)
	phys_len,part_count = parse.(Int,split(ke,","))
	fig = figure()
	imshow(real.(rhos_dict[ke]))
	colorbar()
	title("Momentum Occupation for Phys Dim = $phys_len and Nbosons = $part_count")
	#=
	plot(moms_dict[ke][1]./pi,real.(moms_dict[ke][2])./part_count,label="Exact")
	exact_first_peak = real.(moms_dict[ke][2])[1]
	for km in keys(effective_datadict[1])
		if ke == join([split(km,",")[1],split(km,",")[3]],",")
			nflavs = parse(Int,split(km,",")[2])
			mps_first_peak = real.(effective_datadict[1][km][2])[1]
			plot(effective_datadict[1][km][1]./pi,real.(effective_datadict[1][km][2]).*(exact_first_peak/(mps_first_peak*part_count)),label="$nflavs")
		end
	end
	legend()
	xlabel("p/pi")
	ylabel("Occupation / nparticles")
	title("Momentum Occupation for Phys Dim = $phys_len and Nbosons = $part_count")
	=#
end
end

#

#legend()

#=
for i in 1:length(rhos)
	rho = rhos[i]
	omega = omegas[i]
	#plot_position_occupancy(rho; model_paras...,plot_label="$n_tot",plot_title="L=$L range Density")
	mrez = momentum_dist_1d(rho,n_tot,pcount,pfinal,pinit;if_periodic=if_per,if_log=true,plot_label="$(round(omega,digits=2))",plot_title="L=$L range HarmTrap Mass",if_logscale=false)
end
legend()

=#



=#





























"fin"

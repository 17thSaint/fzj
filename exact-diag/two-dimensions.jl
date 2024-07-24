using Pkg
Pkg.activate(".")
using LinearAlgebra,KrylovKit,Combinatorics,SparseArrays

function find_center()
	all_folders = split(pwd(),"/")
	if "fzj" in all_folders
		return "fzj"
	elseif "local" in all_folders
		return all_folders[findfirst(x -> all_folders[x] == "local",1:length(all_folders))+1]
	elseif "Local" in all_folders
		return all_folders[findfirst(x -> all_folders[x] == "Local",1:length(all_folders))+1]
	else
		println("Not sure where the center is: $(pwd())")
	end
end

function include_other_files(all_files,output_level=0)
	center = find_center()
	get_to_fzj = split(pwd(),center)[1]
	if typeof(all_files) == String
		all_files = [all_files]
	end
	for file in all_files
		occursin("main-git",pwd()) ? include(get_to_fzj * center * "/main-git/" * file) : include(get_to_fzj * center * "/" * file)
		output_level > 0 ? println("Included $file") : nothing
	end
end

include_other_files(["other-funcs/data-storage-funcs.jl"])

######## hopping is wrong because has amplitude even if no particle at starting location #########
#=
function right_x_hopping_matrix(site::Tuple{Int64,Int64},Lx::Int64,Ly::Int64)
    
    # the matrix is mostly close to identity
    totalSites = Lx * Ly
    m::Matrix{Int64} = Matrix(I, totalSites, totalSites)

    # hop to the right and use periodic boundary if necessary, then get snake linear indices
    next_site::Tuple{Int64,Int64} = (mod1(site[1]+1,Lx),site[2])
    site_linear = linear_index(site,Lx,Ly)
    next_site_linear = linear_index(next_site,Lx,Ly)

    # update the matrix
    m[next_site_linear,next_site_linear] = -1
    m[site_linear,site_linear] = 0
    m[next_site_linear,site_linear] = 1

    return m
end

function left_x_hopping_matrix(site::Tuple{Int64,Int64},Lx::Int64,Ly::Int64)
    
    # the matrix is mostly close to identity
    totalSites = Lx * Ly
    m::Matrix{Int64} = Matrix(I, totalSites, totalSites)

    # hop to the left and use periodic boundary if necessary, then get snake linear indices
    next_site::Tuple{Int64,Int64} = (mod1(site[1]-1,Lx),site[2])
    site_linear = linear_index(site,Lx,Ly)
    next_site_linear = linear_index(next_site,Lx,Ly)

    # update the matrix
    m[next_site_linear,next_site_linear] = -1
    m[site_linear,site_linear] = 0
    m[next_site_linear,site_linear] = 1

    return m
end

# hopping matrix for x direction from given site
function x_hopping_matrix(site::Tuple{Int64,Int64},direction::Int64,Lx::Int64,Ly::Int64)
    if direction == 1
        return right_x_hopping_matrix(site,Lx,Ly)
    elseif direction == -1
        return left_x_hopping_matrix(site,Lx,Ly)
    else
        error("direction should be 1 or -1")
    end
end

function up_y_hopping_matrix(site::Tuple{Int64,Int64},Lx::Int64,Ly::Int64)
    
    # the matrix is mostly close to identity
    totalSites = Lx * Ly
    m::Matrix{Int64} = Matrix(I, totalSites, totalSites)

    # hop to the right and use periodic boundary if necessary, then get snake linear indices
    next_site::Tuple{Int64,Int64} = (site[1],mod1(site[2]+1,Ly))
    site_linear = linear_index(site,Lx,Ly)
    next_site_linear = linear_index(next_site,Lx,Ly)

    # update the matrix
    m[next_site_linear,next_site_linear] = -1
    m[site_linear,site_linear] = 0
    m[next_site_linear,site_linear] = 1

    return m
end

function down_y_hopping_matrix(site::Tuple{Int64,Int64},Lx::Int64,Ly::Int64)
    
    # the matrix is mostly close to identity
    totalSites = Lx * Ly
    m::Matrix{Int64} = Matrix(I, totalSites, totalSites)

    # hop to the right and use periodic boundary if necessary, then get snake linear indices
    next_site::Tuple{Int64,Int64} = (site[1],mod1(site[2]-1,Ly))
    site_linear = linear_index(site,Lx,Ly)
    next_site_linear = linear_index(next_site,Lx,Ly)

    # update the matrix
    m[next_site_linear,next_site_linear] = -1
    m[site_linear,site_linear] = 0
    m[next_site_linear,site_linear] = 1

    return m
end

# hopping matrix for y direction from given site
function y_hopping_matrix(site::Tuple{Int64,Int64},direction::Int64,Lx::Int64,Ly::Int64)
    if direction == 1
        return up_y_hopping_matrix(site,Lx,Ly)
    elseif direction == -1
        return down_y_hopping_matrix(site,Lx,Ly)
    else
        error("direction should be 1 or -1")
    end
end

# number operator matrix for given site
function number_operator_matrix(site::Tuple{Int64,Int64},Lx::Int64,Ly::Int64)
    totalSites = Lx * Ly
    m::Matrix{Int64} = zeros(Int64, totalSites, totalSites)

    site_linear = linear_index(site,Lx,Ly)
    m[site_linear,site_linear] = 1

    return m
end

function full_tensor_matrix(ham_comp,which_particle::Int64,N::Int64)
    # different for single particle
    if N == 1
        return ham_comp
    end

    totalSites = size(ham_comp)[1]
    H_part = which_particle == 1 ? ham_comp : Matrix(I,totalSites,totalSites)
    for j in 2:N
        if j == which_particle
            H_part = kron(H_part,ham_comp)
        else
            H_part = kron(H_part,Matrix(I,totalSites,totalSites))
        end
    end
    return H_part
end

function naive_build_ham(lattice_params::Dict,hamilt_params::Dict; kwargs...)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]
    if_periodic_x = lattice_params["if_periodic_x"]
    if_periodic_y = lattice_params["if_periodic_y"]
    totalSites = Lx * Ly

    alpha = get(hamilt_params,"alpha",0.0)
    tx = get(hamilt_params,"tx",1.0)
    ty = get(hamilt_params,"ty",1.0)
    U = get(hamilt_params,"U",zeros(Ly))

    H = sparse(zeros(Float64,totalSites,totalSites))

    for j in 1:Lx
        for s in 1:Ly
            site = (j,s)
            # x-direction hopping (physical)
            coeff = -tx * exp(im*alpha*(s - Ly/2)*2*pi)
            if abs(coeff) > 1e-5
                for dir in [1,-1]
                    dir == -1 ? coeff = conj(coeff) : nothing
                    ham_component = coeff * sparse(x_hopping_matrix(site,dir,Lx,Ly))
                    
                    # skip term if at boundary and no periodic boundary
                    if j == Lx && !if_periodic_x && dir == 1
                        continue
                    elseif j == 1 && !if_periodic_x && dir == -1
                        continue
                    end

                    H += ham_component
                end
            end

            # y-direction hopping (synthetic)
            coeff = -ty
            if abs(coeff) > 1e-5
                for dir in [1,-1]
                    dir == -1 ? coeff = conj(coeff) : nothing
                    ham_component = coeff * sparse(y_hopping_matrix(site,dir,Lx,Ly))

                    # skip term if at boundary and no periodic boundary
                    if s == Ly && !if_periodic_y && dir == 1
                        continue
                    elseif s == 1 && !if_periodic_y && dir == -1
                        continue
                    end

                    H += ham_component
                end
            end

            # interaction
            for ss in 1:Ly
                site2 = (j,ss)
                if s != ss
                    coeff = U[abs(s-ss) + 1]
                    if abs(coeff) > 1e-4
                        ham_component = coeff * sparse(number_operator_matrix(site2,Lx,Ly))

                        H += ham_component
                    end
                end
            end

        end
    end

    return H
end
=#

# builds up the full basis of N hard-core particles in a Lx by Ly lattice (no symmetries yet)
function generate_basis_naive(Lx::Int64,Ly::Int64,N::Int64; kwargs...)

    output_level = get(kwargs,:output_level,1)
    totalSites = Lx * Ly
    n_states = binomial(totalSites,N)
    output_level > 0 ? println("Numbers of basis states: ",n_states) : nothing

    #bit_strings = Array{String,2}(undef,totalSites,n_states)
    bit_states    = Array{Int,2}(undef,   totalSites,n_states)

    basis_dict  = Dict{String,Int}()

    subsets = combinations([1:1:totalSites;],N)

    jj = 1
    for subs in subsets
        #bit_strings[:,jj]  = fill("0",totalSites)
        bit_states[:,jj] = fill(0,  totalSites)
        for en in subs
            #bit_strings[en,jj]  = "1"
            bit_states[en,jj] =  1
        end
        basis_dict[prod(string.(bit_states[:,jj]))] = jj
        jj = jj + 1
        if output_level > 0 && isapprox(jj/length(subsets)*100 % 10,0.0,atol=1e-4)
            println(round(jj/length(subsets)*100,digits=2),"% done.")
        end
    end

    return bit_states,basis_dict
end

generate_basis_naive(L::Int64,N::Int64) = generate_basis_naive(L,L,N)

function get_occupancy_naive(x::Vector{ComplexF64},lattice_params::Dict{String,Any}; kwargs...)
    basis_dict = lattice_params["basis_dict"]
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]

    occupancy = zeros(Float64,Ly,Lx)
    wavefunc = zeros(ComplexF64,Lx*Ly)
    for (key,val) in basis_dict
        wavefunc += x[val] .* [parse(Int64,b) for b in key]
    end
    normalize!(wavefunc)
    wavefunc .*= sqrt(N)

    for j in 1:Lx
        for s in 1:Ly
            site = (j,s)
            occupancy[s,j] = abs2(wavefunc[linear_index(site,Lx,Ly)])
        end
    end

    if_plot = get(kwargs,:if_plot,true)
    if_plot ? plot_occupancy(occupancy; kwargs...) : nothing

    return occupancy
end

function applyHam_naive(which_basis::Int64,lattice_params::Dict,hamilt_params::Dict)
    
    output_states = Array{Int64,1}(undef,0)
    output_weights = Array{ComplexF64,1}(undef,0)
    # get the basis state
    basis_state = lattice_params["full_basis"][:,which_basis]
    basis_dict = lattice_params["basis_dict"]
    
    if_periodic_x = lattice_params["if_periodic_x"]
    if_periodic_y = lattice_params["if_periodic_y"]
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

    tx = hamilt_params["tx"]
    ty = hamilt_params["ty"]
    alpha = hamilt_params["alpha"]
    U = hamilt_params["U"]
    interaction_cutoff = hamilt_params["interaction_cutoff"]

    particle_locations_linear = findall(x->x==1,basis_state)
    particle_locations_coordinate = coordinate.(particle_locations_linear,Lx,Ly)

    # get the hopping weights
    for (idx,starting_site) in enumerate(particle_locations_coordinate)
        
        # x-direction hopping (physical)
        for dir in [1,-1]

            # skip term if at boundary and no periodic boundary
            if starting_site[1] == Lx && !if_periodic_x && dir == 1
                continue
            elseif starting_site[1] == 1 && !if_periodic_x && dir == -1
                continue
            end

            next_site = (mod1(starting_site[1]+dir,Lx),starting_site[2])

            # enforce hard-core constraint
            if next_site in particle_locations_coordinate
                continue
            end

            coeff = -tx * exp(im*alpha*starting_site[2]*2*pi)
            dir == -1 ? coeff = conj(coeff) : nothing
            push!(output_weights,coeff)

            
            output_basis_state = zeros(Int64,length(basis_state)) + basis_state
            output_basis_state[linear_index(starting_site,Lx,Ly)] = 0
            output_basis_state[linear_index(next_site,Lx,Ly)] = 1
            output_basis_state_index = basis_dict[prod(string.(output_basis_state))]
            push!(output_states,output_basis_state_index)
        end

        # y-direction hopping (synthetic)
        for dir in [1,-1]

            # skip term if at boundary and no periodic boundary
            if starting_site[2] == Ly && !if_periodic_y && dir == 1
                continue
            elseif starting_site[2] == 1 && !if_periodic_y && dir == -1
                continue
            end

            next_site = (starting_site[1],mod1(starting_site[2]+dir,Ly))

            # enforce hard-core constraint
            if next_site in particle_locations_coordinate
                continue
            end

            coeff = -ty #* exp(im*alpha*starting_site[1]*2*pi)
            dir == -1 ? coeff = conj(coeff) : nothing
            push!(output_weights,coeff)

            output_basis_state = zeros(Int64,length(basis_state)) + basis_state
            output_basis_state[linear_index(starting_site,Lx,Ly)] = 0
            output_basis_state[linear_index(next_site,Lx,Ly)] = 1
            output_basis_state_index = basis_dict[prod(string.(output_basis_state))]
            push!(output_states,output_basis_state_index)
        end

    end
    #

    # interaction
    lr_dist = sum(U .> interaction_cutoff) - 1
    if length(particle_locations_linear) > 1 && lr_dist > 0
        #println("Doing Interactions")
        for phys_loc in 1:Lx

            # find interacting particles at given physical site
            interacting_particles = findall(x->x[1]==phys_loc,particle_locations_coordinate)
            
            if length(interacting_particles) > 1 # need more than 1 particle to interact
                for i in 1:length(interacting_particles) # loop over all pairs of interacting particles
                    for j in i+1:length(interacting_particles)
                        dist = abs(particle_locations_coordinate[interacting_particles[i]][2] - particle_locations_coordinate[interacting_particles[j]][2])
                        if_periodic_y ? dist = min(dist,Ly-dist) : nothing
                        if dist <= lr_dist && U[dist+1] > interaction_cutoff
                            push!(output_weights,U[dist+1])
                            push!(output_states,which_basis)
                        end
                    end
                end
            end

        end
    end

        
    return output_states,output_weights
end

function buildHam_naive(lattice_params::Dict,hamilt_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    full_basis = lattice_params["full_basis"]

    ham = spzeros(ComplexF64,size(full_basis)[2],size(full_basis)[2])

    for j in 1:size(full_basis)[2]
        output_states,output_weights = applyHam_naive(j,lattice_params,hamilt_params)
        for (idx,state) in enumerate(output_states)
            ham[j,state] += output_weights[idx]
        end
        if output_level > 0 && isapprox(j/size(full_basis)[2]*100 % 10,0.0,atol=1e-4)
            println(round(j/size(full_basis)[2]*100,digits=2),"% done.")
        end
    end

    return ham
end

function get_old_basis_version(basis::Array)
    previous_iteration = sparse(basis[end] == 0 ? [1 0] : [0 1])
    for i in 1:length(basis)-1
        if basis[length(basis)-i] == 0
            previous_iteration = kron(previous_iteration,sparse([1 0]))
        else
            previous_iteration = kron(previous_iteration,sparse([0 1]))
        end
    end
    return previous_iteration
end

function get_old_basis_version(operator::Matrix,which_site::Int64,total_sites::Int64)
    previous_iteration = sparse(which_site == total_sites ? operator : Matrix(I,2,2))
    for i in 1:total_sites-1
        if total_sites - i == which_site
            previous_iteration = kron(previous_iteration,sparse(operator))
        else
            previous_iteration = kron(previous_iteration,sparse(Matrix(I,2,2)))
        end
    end
    return previous_iteration
end

function get_old_basis_version(operator1::Matrix,operator2::Matrix,site1::Int64,site2::Int64,total_sites::Int64)
    
    if site1 == site2
        return get_old_basis_version(operator1*operator2,site1,total_sites)
    end
    
    previous_iteration = sparse(Matrix(I,2,2))
    if site1 == total_sites
        previous_iteration = sparse(operator1)
    elseif site2 == total_sites
        previous_iteration = sparse(operator2)
    end
    for i in 1:total_sites-1
        if total_sites-i == site1
            previous_iteration = kron(previous_iteration,sparse(operator1))
        elseif total_sites-i == site2
            previous_iteration = kron(previous_iteration,sparse(operator2))
        else
            previous_iteration = kron(previous_iteration,sparse(Matrix(I,2,2)))
        end
    end
    return previous_iteration
end

get_old_basis_version(operator1::Matrix,operator2::Matrix,site1::Tuple{Int64,Int64},site2::Tuple{Int64,Int64},Lx::Int64,Ly::Int64) = get_old_basis_version(operator1,operator2,linear_index(site1,Lx,Ly),linear_index(site2,Lx,Ly),Lx*Ly)
get_old_basis_version(operator::Matrix,which_site::Tuple{Int64,Int64},Lx::Int64,Ly::Int64) = get_old_basis_version(operator,linear_index(which_site,Lx,Ly),Lx*Ly)

function convert_new_basis_to_old(basis::Array,lattice_params::Dict)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = length(basis)

    old_basis_version = zeros(Int64,Lx*Ly)
    for i in 1:N
        old_basis_version[basis[i]] = 1
    end

    return old_basis_version
end

function convert_full_basis_new_to_old(lattice_params::Dict)
    full_basis = lattice_params["new_basis"]
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]

    old_basis = Array{Int64,2}(undef,Lx*Ly,binomial(Lx*Ly,N))
    for i in 1:size(full_basis)[2]
        old_basis[:,i] = convert_new_basis_to_old(full_basis[:,i],lattice_params)
    end

    return old_basis
end

function change_of_basis_matrix(lattice_params::Dict{String,Any})
    full_basis = lattice_params["full_basis"]
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]

    total_sites = Lx * Ly
    change_of_basis = spzeros(ComplexF64,2^(Lx*Ly),size(full_basis)[2])

    for j in 1:size(full_basis)[2]
        change_of_basis[:,j] = transpose(get_old_basis_version(full_basis[:,j]))
    end

    return change_of_basis#,sparse(pinv(Matrix(change_of_basis)))
end

function hopping_matrix_naive(site1::Tuple{Int64,Int64},site2::Tuple{Int64,Int64},lattice_params::Dict{String,Any})
    full_basis = lattice_params["full_basis"]
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    change_of_basis = lattice_params["change_of_basis"]
    change_of_basis_inv = lattice_params["change_of_basis_inv"]

    creation = [0 0; 1 0]
    annihilation = [0 1; 0 0]

    hopping_old_basis = get_old_basis_version(creation,annihilation,site1,site2,Lx,Ly)

    new_basis_hopping = change_of_basis_inv * hopping_old_basis * change_of_basis
    return new_basis_hopping
end

function density_matrix_naive(x::Vector{ComplexF64},lattice_params::Dict{String,Any})
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

    rho = Array{ComplexF64,2}(undef,Lx*Ly,Lx*Ly)

    for i in 1:Lx*Ly
        for j in 1:i
            rho[i,j] = conj(transpose(x)) * hopping_matrix_naive(coordinate(i,Lx,Ly),coordinate(j,Lx,Ly),lattice_params) * x
            rho[j,i] = conj(rho[i,j])
        end
    end

    return rho
end

function buildHopping_old(lattice_params::Dict,site1::Int64,site2::Int64; kwargs...)
    output_level = get(kwargs,:output_level,1)
    full_basis = lattice_params["full_basis"]
    basis_dict = lattice_params["basis_dict"]

    hop = spzeros(ComplexF64,size(full_basis)[2],size(full_basis)[2])

    if site1 == site2
        for j in 1:size(full_basis)[2]
            this_basis_state = full_basis[:,j]
            if this_basis_state[site1] == 1
                output_state = basis_dict[prod(string.(this_basis_state))]
                hop[j,output_state] = 1.0+0.0*im
            end
        end
    else
        for j in 1:size(full_basis)[2]
            this_basis_state = full_basis[:,j]
            if this_basis_state[site1] == 1 && this_basis_state[site2] == 0 # check if particle at starting location and check hard-core constraint
                
                output_state = this_basis_state .+ 0
                output_state[site1] = 0
                output_state[site2] = 1

                output_state = basis_dict[prod(string.(output_state))]
                hop[j,output_state] = 1.0+0.0*im
            end
        end
    end

    return hop
end

function hopping_probability_old(x::Vector{ComplexF64},site1::Tuple{Int64,Int64},site2::Tuple{Int64,Int64},lattice_params::Dict{String,Any}; kwargs...)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

    s1_linear = linear_index(site1,Lx,Ly)
    s2_linear = linear_index(site2,Lx,Ly)

    hopping_operator = buildHopping_old(lattice_params,s1_linear,s2_linear; kwargs...)
    hopping_prob = conj(transpose(x)) * hopping_operator * x

    return hopping_prob
end

function density_matrix_old(x::Vector{ComplexF64},lattice_params::Dict{String,Any}; kwargs...)
    output_level = get(kwargs,:output_level,1)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

    rho = Array{ComplexF64,2}(undef,Lx*Ly,Lx*Ly)

    start_time = time()
    for i in 1:Lx*Ly
        for j in 1:i
            rho[i,j] = hopping_probability_old(x,coordinate(i,Lx,Ly),coordinate(j,Lx,Ly),lattice_params; kwargs...)
            rho[j,i] = conj(rho[i,j])
        end
        output_level > 0 ? println(round(i/(Lx*Ly)*100,digits=2),"% done.") : nothing
    end

    output_level > 0 ? println("Density Matrix: Elapsed time: ",time()-start_time) : nothing

    return rho
end

function density_matrix_slow(x::Vector{ComplexF64},lattice_params::Dict{String,Any}; kwargs...)
    output_level = get(kwargs,:output_level,1)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    if_periodic_x = lattice_params["if_periodic_x"]
    if_periodic_y = lattice_params["if_periodic_y"]
    if_periodic = (if_periodic_x,if_periodic_y)
    edges = (Lx,Ly)
    N = lattice_params["N"]
    full_basis = lattice_params["full_basis"]
    dimHilb = size(full_basis)[2]
    output_level > 0 ? println("Hilbert Space Dimension: ",dimHilb) : nothing

    rho = Array{ComplexF64,2}(undef,Lx*Ly,Lx*Ly)
    all_hopping_operators = spzeros(Int64,dimHilb*Lx*Ly,dimHilb*Lx*Ly)

    start_time = time()
    for i in 1:dimHilb
        this_basis = full_basis[:,i]
        for n in 1:N
            for j in 1:Lx*Ly

                # enforce hard-core constraint but keep diagonal elements
                if j in this_basis && this_basis[n] != j
                    continue
                end

                # find basis of next configuration
                next_basis = this_basis .+ 0
                next_basis[n] = j
                sort!(next_basis,rev=true)
                next_basis_index = find_basis_index(next_basis)

                # add element to hopping operator
                local_indices = (i,next_basis_index)
                shift = ((this_basis[n]-1)*dimHilb,(j-1)*dimHilb)
                index_result = local_indices .+ shift
                all_hopping_operators[index_result[1],index_result[2]] = 1
            end
        end
    end

    for i in 1:Lx*Ly
        for j in 1:i
            rho[i,j] = conj(transpose(x)) * all_hopping_operators[(i-1)*dimHilb+1:i*dimHilb,(j-1)*dimHilb+1:j*dimHilb] * x
            rho[j,i] = conj(rho[i,j])
        end
        output_level > 0 ? println(round(i/(Lx*Ly)*100,digits=2),"% done.") : nothing
    end

    output_level > 0 ? println("Density Matrix: Elapsed time: ",time()-start_time) : nothing

    return rho

end

######## This is the better working way to do things ########

# finds the linear index assuming jump snake mapping with site 1 at bottom left corner
function linear_index(site::Tuple{Int64,Int64},Lx::Int64,Ly::Int64)
    return (site[2] - 1)*Lx + site[1]
end

# finds the site assuming jump snake mapping with site 1 at bottom left corner
function coordinate(site::Int64,Lx::Int64,Ly::Int64)
    x = mod1(site,Lx)
    y = Int((site - x) / Lx + 1)
    return (x,y)
end

function find_basis_index(basis::Vector{Int64})
    p = basis[1]
    lower_limit = 1 + binomial(p-1,length(basis))
    for i in 1:length(basis)-1
        q = basis[i+1]
        lower_limit += binomial(q-1,length(basis)-i)
    end
    
    return lower_limit
end

function single_particle_basis(start_point::Int64,end_point::Int64)
    return [[i] for i in start_point:end_point]
end

function two_particle_basis(Lx::Int64,Ly::Int64)
    sp_basis = single_particle_basis(1,Lx*Ly)
    two_particle_basis = Array{Vector{Int64},1}(undef,0)
    for i in 1:length(sp_basis)
        p = sp_basis[i][1]
        for q in 1:p-1
            push!(two_particle_basis,[p,q])
        end
    end
    return two_particle_basis
end

function save_basis(full_basis::Array{Int64,2},N::Int64,Lx::Int64,Ly::Int64,dataloc::String)
    data_dict = Dict([("full_basis",full_basis)])
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N)])
    filename = "basis-N-"*string(N)*"-Lx-"*string(Lx)*"-Ly-"*string(Ly)*".jld2"
    full_loc = join([dataloc,filename],"/")
    write_data_jld2(full_loc,data_dict,metadata_dict)
end

function n_particle_basis(N::Int64,Lx::Int64,Ly::Int64; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if_save_data = get(kwargs,:if_save_data,true)
    dataloc = get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag"))
    if_find_existing = get(kwargs,:if_find_existing,true)

    if if_find_existing
        metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N)])
        if_exists,data = check_data_exists(metadata_dict,"basis"; location=dataloc,output_level=output_level)
        if if_exists
            output_level > 0 ? println("Found existing file with basis data") : nothing
            output_level > 0 ? println("Basis has ",size(data[1]["full_basis"],2)," states") : nothing
            return data[1]["full_basis"]
        end
    end

    n2_basis = two_particle_basis(Lx,Ly)
    for i in 3:N
        n_particle_basis = Array{Vector{Int64},1}(undef,0)
        for j in 1:length(n2_basis)
            q = n2_basis[j][end]
            for r in 1:q-1
                push!(n_particle_basis,[n2_basis[j];[r]])
            end
        end
        n2_basis = n_particle_basis
    end

    output_level > 0 ? println("Basis has ",length(n2_basis)," states") : nothing

    full_basis = Array{Int64,2}(undef,N,length(n2_basis))
    for i in 1:length(n2_basis)
        full_basis[:,i] = zeros(Int64,N) + n2_basis[i]
    end

    if_save_data ? save_basis(full_basis,N,Lx,Ly,dataloc) : nothing

    return full_basis
end

function n_particle_basis(lattice_params::Dict; kwargs...)
    N = lattice_params["N"]
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    return n_particle_basis(N,Lx,Ly; kwargs...)
end

function long_range_scaling(x_final::Int64,virt_edge_length::Int64,initial_strength::Float64; kwargs...)

	trunc = get(kwargs, :trunc_digits, 5)
	scaling_func = get(kwargs, :scaling, "flat")
	
	strengths = zeros(virt_edge_length)

    if x_final == 0.0 || initial_strength == 0.0
        strengths[1] = initial_strength != 0.0 ? initial_strength : 1.0
        return strengths
    end
	
	if scaling_func == "flat"
		strengths[1:x_final+1] .= initial_strength
    elseif scaling_func == "gaussian"
        sigma = get(kwargs, :sigma, 1.0)
        strengths = map(0:virt_edge_length-1) do x
            #(initial_strength/(sqrt(2*pi)*sigma)) * exp(-x^2/(2*sigma^2))
            initial_strength * exp(-x^2/(2*sigma^2))
        end
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

function applyHam(which_basis::Int64,lattice_params::Dict,hamilt_params::Dict)
    
    output_states = Array{Int64,1}(undef,0)
    output_weights = Array{ComplexF64,1}(undef,0)
    # get the basis state
    basis_state = lattice_params["full_basis"][:,which_basis]
    
    if_periodic_x = lattice_params["if_periodic_x"]
    if_periodic_y = lattice_params["if_periodic_y"]
    if_per = (if_periodic_x,if_periodic_y)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    Ls = (Lx,Ly)
    twist_angle = lattice_params["twist_angle"]

    tx = hamilt_params["tx"]
    ty = hamilt_params["ty"]
    alpha = hamilt_params["alpha"]
    flux_direction = hamilt_params["flux_direction"]
    U = hamilt_params["U"]
    interaction_cutoff = hamilt_params["interaction_cutoff"]
    which_dir = hamilt_params["which_dir"]

    particle_locations_linear = basis_state
    particle_locations_coordinate = coordinate.(particle_locations_linear,Lx,Ly)

    # get the hopping weights
    for (idx,starting_site) in enumerate(particle_locations_coordinate)
        
        # x-direction hopping (physical)
        for dir in [1,-1]

            # skip term if at boundary and no periodic boundary
            if starting_site[1] == Lx && dir == 1
                twist = true
                if !if_periodic_x
                    continue
                end
            elseif starting_site[1] == 1 && dir == -1
                twist = true
                if !if_periodic_x
                    continue
                end
            else
                twist = false
            end

            next_site = (mod1(starting_site[1]+dir,Lx),starting_site[2])

            # enforce hard-core constraint
            if next_site in particle_locations_coordinate
                continue
            end

            coeff = -tx
            flux_direction == "x" ? coeff *= exp(im*alpha*starting_site[2]*2*pi) : nothing

            # apply twist at boundary
            if twist
                #println("Applying twist at $starting_site moving to $next_site along X of $(twist_angle[1])")
                coeff *= exp(im*2*pi*twist_angle[1])
            end
            
            dir == -1 ? coeff = conj(coeff) : nothing
            push!(output_weights,round(real(coeff),digits=10) + im*round(imag(coeff),digits=10))

            
            output_basis_state = zeros(Int64,length(basis_state)) + basis_state
            output_basis_state[idx] = linear_index(next_site,Lx,Ly)
            sort!(output_basis_state,rev=true)
            output_basis_state_index = find_basis_index(output_basis_state)
            push!(output_states,output_basis_state_index)

        end

        # y-direction hopping (synthetic)
        for dir in [1,-1]

            # skip term if at boundary and no periodic boundary
            if starting_site[2] == Ly && dir == 1
                twist = true
                if !if_periodic_y
                    continue
                end
            elseif starting_site[2] == 1 && dir == -1
                twist = true
                if !if_periodic_y
                    continue
                end
            else
                twist = false
            end

            next_site = (starting_site[1],mod1(starting_site[2]+dir,Ly))

            # enforce hard-core constraint
            if next_site in particle_locations_coordinate
                continue
            end

            coeff = -ty 
            flux_direction == "y" ? coeff *= exp(im*alpha*starting_site[1]*2*pi) : nothing

            # apply twist at boundary
            if twist
                #println("Applying twist at $starting_site moving to $next_site along Y of $(twist_angle[2])")
                coeff *= exp(im*2*pi*twist_angle[2])
            end

            dir == -1 ? coeff = conj(coeff) : nothing
            push!(output_weights,round(real(coeff),digits=10) + im*round(imag(coeff),digits=10))

            output_basis_state = zeros(Int64,length(basis_state)) + basis_state
            output_basis_state[idx] = linear_index(next_site,Lx,Ly)
            sort!(output_basis_state,rev=true)
            output_basis_state_index = find_basis_index(output_basis_state)
            push!(output_states,output_basis_state_index)
        end

    end
    #

    # interaction
    lr_dist = sum(abs.(U) .> interaction_cutoff) - 1
    if length(particle_locations_linear) > 1 && lr_dist > 0
        #println("Doing Interactions")
        if which_dir == "virt"
            which_loc = 1
            other_loc = 2
        elseif which_dir == "phys"
            which_loc = 2
            other_loc = 1
        else
            error("which_dir must be either 'virt' or 'phys'")
        end
        for loc in 1:Ls[which_loc]

            # find interacting particles at given physical site
            interacting_particles = findall(x->x[which_loc]==loc,particle_locations_coordinate)

            if length(interacting_particles) > 1 # need more than 1 particle to interact
                for i in 1:length(interacting_particles) # loop over all pairs of interacting particles
                    for j in i+1:length(interacting_particles)
                        dist = abs(particle_locations_coordinate[interacting_particles[i]][other_loc] - particle_locations_coordinate[interacting_particles[j]][other_loc])
                        #if if_per[other_loc] && which_loc == 2
                        #    dist = min(dist,Ls[other_loc]-dist)
                        #end
                        if dist <= lr_dist && abs(U[dist+1]) > interaction_cutoff
                            #println("Interacting between ",particle_locations_coordinate[interacting_particles[i]]," and ",particle_locations_coordinate[interacting_particles[j]])
                            push!(output_weights,U[dist+1])
                            push!(output_states,which_basis)
                        end
                    end
                end
            end

        end
    end

    if hamilt_params["disorder_strength"] != 0.0
        local_disorder_strength = rand() * hamilt_params["disorder_strength"] * 2 - hamilt_params["disorder_strength"]
        push!(output_weights,local_disorder_strength)
        push!(output_states,which_basis)
    end

    return output_states,output_weights
end

function buildHam(lattice_params::Dict,hamilt_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    full_basis = lattice_params["full_basis"]

    ham = spzeros(ComplexF64,size(full_basis)[2],size(full_basis)[2])

    for j in 1:size(full_basis)[2]
        output_states,output_weights = applyHam(j,lattice_params,hamilt_params)
        for (idx,state) in enumerate(output_states)
            ham[j,state] += output_weights[idx]
        end
        if output_level > 0 && isapprox(j/size(full_basis)[2]*100 % 10,0.0,atol=1e-4)
            println(round(j/size(full_basis)[2]*100,digits=2),"% done.")
        end
    end

    return ham
end

function buildHopping(lattice_params::Dict,site1::Int64,site2::Int64; kwargs...)
    output_level = get(kwargs,:output_level,1)
    full_basis = lattice_params["full_basis"]

    hop = spzeros(Int64,size(full_basis)[2],size(full_basis)[2])

    if site1 == site2
        for j in 1:size(full_basis)[2]
            site1 in full_basis[:,j] ? hop[j,j] += 1.0+0.0*im : nothing
        end
    else
        for j in 1:size(full_basis)[2]
            this_basis_state = full_basis[:,j]
            if site1 in this_basis_state && !(site2 in this_basis_state) # check if particle at starting location and check hard-core constraint
                output_state = zeros(Int64,length(this_basis_state)) + this_basis_state
                output_state[findfirst(x->this_basis_state[x]==site1,1:length(this_basis_state))] = site2
                sort!(output_state,rev=true)
                hop[j,find_basis_index(output_state)] = 1
            end
        end
    end

    return hop
end

function hopping_probability(x::Vector{ComplexF64},site1::Tuple{Int64,Int64},site2::Tuple{Int64,Int64},lattice_params::Dict{String,Any}; kwargs...)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

    s1_linear = linear_index(site1,Lx,Ly)
    s2_linear = linear_index(site2,Lx,Ly)

    hopping_operator = buildHopping(lattice_params,s1_linear,s2_linear; kwargs...)
    hopping_prob = conj(transpose(x)) * hopping_operator * x

    return hopping_prob
end

function density_matrix(x::Vector{ComplexF64},lattice_params::Dict{String,Any}; kwargs...)
    output_level = get(kwargs,:output_level,1)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

    rho = Array{ComplexF64,2}(undef,Lx*Ly,Lx*Ly)

    start_time = time()
    for i in 1:Lx*Ly
        for j in 1:i
            rho[i,j] = hopping_probability(x,coordinate(i,Lx,Ly),coordinate(j,Lx,Ly),lattice_params; kwargs...)
            rho[j,i] = conj(rho[i,j])
        end
        #rho[i,:] ./= abs(hopping_probability(x,coordinate(i,Lx,Ly),coordinate(i,Lx,Ly),lattice_params; kwargs...))
        output_level > 0 ? println(round(i/(Lx*Ly)*100,digits=2),"% done.") : nothing
    end

    output_level > 0 ? println("Density Matrix: Elapsed time: ",time()-start_time) : nothing

    return rho
end

function site_occupation(x::Vector{ComplexF64},site::Tuple{Int64,Int64},lattice_params::Dict{String,Any})
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    site_linear = linear_index(site,Lx,Ly)
    return conj(transpose(x)) * buildHopping(lattice_params,site_linear,site_linear) * x
end

function find_eigenstates(nev::Int,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if_densmat = get(kwargs,:if_densmat,true)
    if_save_data = get(kwargs,:if_save_data,false)
    if_exact = get(kwargs,:if_exact,false)
    if_function = get(kwargs,:if_function,true)

    metadata_dict = merge(merge(lattice_params,hamilt_params),named_tuple_to_dict(kwargs))

    metadata_displaying::Dict = copy(metadata_dict)
    delete!(metadata_displaying,"full_basis")
    output_level > 0 ? display(metadata_displaying) : nothing

    dimHilb = size(lattice_params["full_basis"])[2]
    
    start_time = time()
    if if_function
        
        function ham_func(x)
            output_vector::Vector{ComplexF64} = zeros(ComplexF64,dimHilb)
            for i in 1:dimHilb
                output_states,output_weights = applyHam(i,lattice_params,hamilt_params)
                for (idx,state) in enumerate(output_states)
                    output_vector[state] += output_weights[idx]*x[i]
                end
            end
            return output_vector
        end

        kdim = get(kwargs,:kdim,nev+10)
        x0 = rand(Float64,size(lattice_params["full_basis"])[2])
        rez = eigsolve(ham_func,x0,nev,:SR; krylovdim=kdim)
    else
        H = buildHam(lattice_params,hamilt_params; output_level)
        #display(H)
        metadata_dict["H"] = H
        output_level > 0 ? println("Sparsity = ",SparseArrays.nnz(H)/size(H)[1]^2) : nothing

        if if_exact
            everything = eigen(Matrix(H))
            rez = (everything.values,everything.vectors)
        else
            x0 = rand(Float64,size(lattice_params["full_basis"])[2])
            rez = eigsolve(H,x0,nev,:SR,Lanczos())
        end
    end
    output_level > 0 ? println("Ground State: Elapsed time: ",time()-start_time) : nothing
    metadata_dict["runtime"] = time()-start_time

    sorted_indices = sortperm(real.(rez[1]))
    states = rez[2][sorted_indices][1:nev]
    nrgs = rez[1][sorted_indices][1:nev]

    rhos = []
    if if_densmat
        for state in states
            append!(rhos,[density_matrix(state,lattice_params)])
        end
    end

    if_save_data ? save_eigenstates(states,rhos,nrgs,metadata_dict) : nothing

    if if_function
        return states,real.(nrgs),rhos
    else
        return states,real.(nrgs),rhos,H
    end
end

function rerun_eigenstates(nev::Int,lattice_params::Dict,hamilt_params::Dict,metadata::Dict,data_dict::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if_densmat = get(kwargs,:if_densmat,true)
    if_save_data = get(kwargs,:if_save_data,false)
    if_exact = get(kwargs,:if_exact,false)
    if_function = get(kwargs,:if_function,true)


    metadata_displaying::Dict = copy(metadata)
    delete!(metadata_displaying,"full_basis")
    delete!(metadata_displaying,"H")
    output_level > 0 ? display(metadata_displaying) : nothing
    
    dimHilb = size(lattice_params["full_basis"])[2]

    start_time = time()
    if if_function

        function ham_func(x)
            output_vector::Vector{ComplexF64} = zeros(ComplexF64,dimHilb)
            for i in 1:dimHilb
                output_states,output_weights = applyHam(i,lattice_params,hamilt_params)
                for (idx,state) in enumerate(output_states)
                    output_vector[state] += output_weights[idx]*x[i]
                end
            end
            return output_vector
        end

        kdim = get(kwargs,:kdim,nev+10)
        x0 = rand(Float64,size(lattice_params["full_basis"])[2])
        rez = eigsolve(ham_func,x0,nev,:SR; krylovdim=kdim)
    else
        H = metadata["H"]
        output_level > 0 ? println("Sparsity = ",nnz(H)/size(H)[1]^2) : nothing
        if if_exact
            everything = eigen(Matrix(H))
            rez = (everything.values,everything.vectors)
        else
            x0 = rand(Float64,size(lattice_params["full_basis"])[2])
            rez = eigsolve(H,x0,nev,:SR,Lanczos())
        end
    end
    previous_nev = metadata["nev"]
    output_level > 0 ? println("Ground State: Elapsed time: ",time()-start_time) : nothing

    sorted_indices = sortperm(real.(rez[1]))
    states = rez[2][sorted_indices][1:nev]
    nrgs = rez[1][sorted_indices][1:nev]

    rhos = []
    if if_densmat
        for (idx,state) in enumerate(states)
            idx <= previous_nev ? append!(rhos,[data_dict["densmat"][idx]]) : append!(rhos,[density_matrix(state,lattice_params)])
        end
    end

    if if_save_data
        full_loc = join([metadata["dataloc"],metadata["filename"]],"/")
        modify_data_jld2("nev",nev,full_loc,"metadata")
        new_data_dict = Dict([("state",states),("nrg",nrgs),("densmat",rhos)])
        modify_data_jld2(new_data_dict,full_loc,"all_data")
    end

    return states,nrgs,rhos
end

find_ground_state(lattice_params::Dict,hamilt_params::Dict; kwargs...) = find_eigenstates(1,lattice_params,hamilt_params; kwargs...)

function save_eigenstates(states,densmats,nrgs,metadata::Dict)
    dataloc = get(metadata,"dataloc",get_folder_location("cluster-data/exact-diag"))
    data_dict = Dict([("state",states),("nrg",nrgs),("densmat",densmats)])
    lattice_params,hamilt_params = make_latticehamilt_params_from_metadata(metadata)
    filename_dict = make_filename_dict(lattice_params,hamilt_params)
    filename = join(["ed",make_parameters_filename(filename_dict)],"-")
    metadata["filename"] = filename
    full_loc = join([dataloc,filename],"/")
    println("Filename: ",full_loc)
    metadata["full_basis"] = nothing
    write_data_jld2(full_loc,data_dict,metadata)
end

function make_latticehamilt_params_from_metadata(metadata::Dict)
    hamilt_params = Dict([("tx",metadata["tx"]),("ty",metadata["ty"]),("hopping_anisotropy",metadata["hopping_anisotropy"]),("alpha",metadata["alpha"]),("U",metadata["U"]),("interaction_cutoff",metadata["interaction_cutoff"]),("which_dir",metadata["which_dir"])])
    lattice_params = get_lattice_params_from_metadata(metadata)
    flux_dir = get(metadata,"flux_direction","x")
    if lattice_params["if_periodic_y"] && !lattice_params["if_periodic_x"]
        flux_dir = "y"
    elseif !lattice_params["if_periodic_y"] && lattice_params["if_periodic_x"]
        flux_dir = "x"
    end
    hamilt_params["flux_direction"] = flux_dir
    hamilt_params["disorder_strength"] = "disorder_strength" in keys(metadata) ? metadata["disorder_strength"] : 0.0
    return lattice_params,hamilt_params
end

function physical_correlation(densmat::Array{ComplexF64,2},Lx::Int64,Ly::Int64; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    phys_corrs = Array{Float64,2}(undef,Lx,Ly)
    for j in 1:Lx
        for s in 1:Ly
            phys_corr = densmat[linear_index((j,s),Lx,Ly),linear_index((1,s),Lx,Ly)]
            phys_corr /= sqrt(densmat[linear_index((j,s),Lx,Ly),linear_index((j,s),Lx,Ly)] * densmat[linear_index((1,s),Lx,Ly),linear_index((1,s),Lx,Ly)])
            phys_corrs[j,s] = abs(phys_corr)
        end
    end

    if_plot ? plot_physical_correlation(phys_corrs; kwargs...) : nothing

    return phys_corrs
end

function synthetic_correlation(densmat::Array{ComplexF64,2},Lx::Int64,Ly::Int64; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    syn_corrs = Array{Float64,2}(undef,Lx,Ly)
    for j in 1:Lx
        for s in 1:Ly
            syn_corr = densmat[linear_index((j,s),Lx,Ly),linear_index((j,1),Lx,Ly)]
            syn_corr /= sqrt(densmat[linear_index((j,s),Lx,Ly),linear_index((j,s),Lx,Ly)] * densmat[linear_index((j,1),Lx,Ly),linear_index((j,1),Lx,Ly)])
            syn_corrs[j,s] = abs(syn_corr)
        end
    end

    if_plot ? plot_synthetic_correlation(syn_corrs; kwargs...) : nothing

    return syn_corrs
end

function plot_synthetic_correlation(syn_corrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(syn_corrs)[1]
        plot(0:size(syn_corrs)[2]-1,syn_corrs[i,:],"-p",label="$i")
    end
    xlabel("Synthetic Distance")
    ylabel("Correlation")
    title("Synthetic Correlation "*plot_title)
    legend()
    yscale("log")
    return nothing
end

function plot_physical_correlation(phys_corrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(phys_corrs)[2]
        plot(0:size(phys_corrs)[1]-1,phys_corrs[:,i],"-p",label="$i")
    end
    xlabel("Physical Distance")
    ylabel("Correlation")
    title("Physical Correlation "*plot_title)
    legend()
    yscale("log")
    return nothing
end

function physical_current(densmat::Array{ComplexF64,2},lattice_params::Dict{String,Any}; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    if_periodic_x = lattice_params["if_periodic_x"]

    xlen = if_periodic_x ? Lx : Lx-1

    currents = Array{Float64,2}(undef,xlen,Ly)
    for s in 1:Ly
        for j in 1:xlen
            site1 = linear_index((j,s),Lx,Ly)
            site2 = linear_index((mod1(j+1,Lx),s),Lx,Ly)
            current_val = imag(densmat[site1,site2] - densmat[site2,site1])
            current_normalization = densmat[site1,site1] + densmat[site2,site2]
            current_val /= current_normalization
            currents[j,s] = real(current_val)
        end
    end

    if_plot ? plot_physical_current(currents; kwargs...) : nothing

    return currents
end

function synthetic_current(densmat::Array{ComplexF64,2},lattice_params::Dict{String,Any}; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    if_periodic_y = lattice_params["if_periodic_y"]

    ylen = if_periodic_y ? Ly : Ly-1

    currents = Array{Float64,2}(undef,Lx,ylen)
    for j in 1:Lx
        for s in 1:ylen
            site1 = linear_index((j,s),Lx,Ly)
            site2 = linear_index((j,mod1(s+1,Ly)),Lx,Ly)
            current_val = imag(densmat[site1,site2] - densmat[site2,site1])
            current_normalization = densmat[site1,site1] + densmat[site2,site2]
            current_val /= current_normalization
            currents[j,s] = real(current_val)
        end
    end

    if_plot ? plot_synthetic_current(currents; kwargs...) : nothing

    return currents
end

function plot_synthetic_current(currents::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(currents)[2]
        plot(1:size(currents)[1],currents[:,i],"-p",label="$i")
    end
    xlabel("Physical Site")
    ylabel("Current")
    title("Synthetic Current "*plot_title)
    legend()
    return nothing
end

function plot_physical_current(currents::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(currents)[1]
        plot(1:size(currents)[2],currents[i,:],"-p",label="$i")
    end
    xlabel("Synthetic Site")
    ylabel("Current")
    title("Physical Current "*plot_title)
    legend()
    return nothing
end

function check_fluxes(alpha,Lx::Int64,Ly::Int64,if_periodic_x::Bool,if_periodic_y::Bool,flux_direction::String)
    if alpha == 0.0
        return nothing
    end
    if alpha > 0.4
        error("Alpha is too large: ",alpha)
    end
    x_shift,y_shift = !if_periodic_x, !if_periodic_y
    num_fluxes = round(alpha*(Lx - x_shift) * (Ly - y_shift),digits=5)
    println("Number of Fluxes = ",num_fluxes," for Lx = ",Lx," and Ly = ",Ly)
    if !isinteger(num_fluxes)
        error("Number of fluxes is not an integer")
    end

    println("Checking fluxes only along Gauge Direction")
    if flux_direction == "y"
        if if_periodic_x && !isinteger(num_fluxes/Ly)
            if if_periodic_y && isinteger(num_fluxes/Lx)
                flux_direction = "x"
                println("Fluxes don't fit, changing to X direction")
            else
                error("Number of fluxes is not an integer multiple of Lx")
            end
        end
    elseif flux_direction == "x"
        if if_periodic_y && !isinteger(num_fluxes/Lx)
            if if_periodic_x && isinteger(num_fluxes/Ly)
                flux_direction = "y"
                println("Fluxes don't fit, changing to Y direction")
            else
                error("Number of fluxes is not an integer multiple of Ly")
            end
        end
    else
        error("Flux direction is not valid")
    end

    #=println("Not properly checking fluxes")
    println("Checking fluxes only along Gauge Direction")
    if if_periodic_x && flux_direction == "y" && !isinteger(num_fluxes/Ly)
        error("Number of fluxes is not an integer multiple of Lx")
    end

    if if_periodic_y && flux_direction == "x" && !isinteger(num_fluxes/Lx)
        error("Number of fluxes is not an integer multiple of Ly")
    end=#

    return flux_direction
end

function get_occupancy(x::Vector{ComplexF64},lattice_params::Dict{String,Any}; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    full_basis = lattice_params["full_basis"]

    occs = zeros(Float64,Ly,Lx)
    
    for i in 1:length(x)
        basis_state = full_basis[:,i]
        for n in 1:length(basis_state)
            site = coordinate(basis_state[n],Lx,Ly)
            occs[site[2],site[1]] += abs(x[i])^2
        end
    end

    if_plot ? plot_occupancy(occs; kwargs...) : nothing

    return occs
end

function get_occupancy(rho::Array{ComplexF64,2},lattice_params::Dict{String,Any}; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

    occs = zeros(Float64,Ly,Lx)
    all_occs = diag(rho)

    for (idx,val) in enumerate(all_occs)
        site = coordinate(idx,Lx,Ly)
        occs[site[2],site[1]] = real(val)
    end

    if_plot ? plot_occupancy(occs; kwargs...) : nothing

    return occs
end

function plot_occupancy(exp_occ; kwargs...)
	fig = figure()
	imshow(exp_occ)
	colorbar()
	plot_title = get(kwargs, :plot_title, "")
	title_string = "Occupancy, " * plot_title
	title(title_string)
	ylabel("Synthetic")
	xlabel("Physical")

    return nothing
end

function make_filename_dict(lattice_params::Dict,hamilt_params::Dict)
    if hamilt_params["U"][2] == 0.0
        intstren = 0.0
    else
        intstren = hamilt_params["U"][1]
    end
    fdict = Dict([("Lx",lattice_params["Lx"]),("Ly",lattice_params["Ly"]),("N",lattice_params["N"]),("alpha",hamilt_params["alpha"]),("hopping_anisotropy",hamilt_params["tx"]/hamilt_params["ty"]),("interaction_strength",intstren),("if_periodic_x",lattice_params["if_periodic_x"]),("if_periodic_y",lattice_params["if_periodic_y"])])
    if lattice_params["twist_angle"] != [0.0,0.0]
        fdict["twist_angle1"] = lattice_params["twist_angle"][1]
        fdict["twist_angle2"] = lattice_params["twist_angle"][2]
    end
    if hamilt_params["disorder_strength"] != 0.0
        fdict["disorder_strength"] = hamilt_params["disorder_strength"]
    end
    return fdict
end

function get_lattice_params_from_metadata(metadata::Dict)
    lat_paras = Dict([("Lx",metadata["Lx"]),("Ly",metadata["Ly"]),("N",metadata["N"]),("if_periodic_x",metadata["if_periodic_x"]),("if_periodic_y",metadata["if_periodic_y"]),("full_basis",metadata["full_basis"]),("twist_angle",metadata["twist_angle"])])
    if isnothing(lat_paras["full_basis"])
        lat_paras["full_basis"] = n_particle_basis(lat_paras)
    end
    return lat_paras
end

function hh_gap_exact(hopping_anisotropy::Float64,alpha::Float64)
    a2 = 1/(4*(sin(pi*alpha)^2))
    a1 = (1+2*(1/hopping_anisotropy^2)*(a2^2))^-0.5
    val = 4*alpha*(1/hopping_anisotropy^2)*(a1^4)*(a2^2)
    if val > 0.5
        return 1/(1+2*(1/hopping_anisotropy^2)*(a2^2))
    else
        return val
    end
end

function hh_gap_fit(x,p)
    return p[1] .* hh_gap_exact.(x,0.25) .+ p[2]
end

function distance_correlation(rho::Matrix,lattice_params::Dict,direction::String="x")
    if_periodic_x = lattice_params["if_periodic_x"]
    if_periodic_y = lattice_params["if_periodic_y"]
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

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
            s1 = direction == "x" ? linear_index((x1,y1),Lx,Ly) : linear_index((y1,x1),Lx,Ly)
            for x2 in 1:len-1
                if x1 == x2
                    continue
                end
                s2 = direction == "x" ? linear_index((x2,y1),Lx,Ly) : linear_index((y1,x2),Lx,Ly)
                dist_corrs[Int(abs(x1-x2))] += abs(rho[s1,s2])
                corr_counts[Int(abs(x1-x2))] += 1
            end
        end
    end

    dist_corrs ./= corr_counts

    return dist_corrs
end

distance_correlation(x::Vector{ComplexF64},lattice_params::Dict,direction::String="x") = distance_correlation(density_matrix(x,lattice_params),lattice_params,direction)

function find_dist(p1::Tuple{Int,Int}, p2::Tuple{Int,Int}, size::Tuple{Int,Int}, periodic::Tuple{Bool,Bool}=(false, false))
    dx = abs(p1[1] - p2[1])
    dy = abs(p1[2] - p2[2])

    if periodic[1]
        dx = min(dx, size[1] - dx)
    end

    if periodic[2]
        dy = min(dy, size[2] - dy)
    end

    return sqrt(dx^2 + dy^2),(dx,dy),(p1[1]-p2[1],p1[2]-p2[2])
end

function cdw_structure_factor(rho::Array{ComplexF64,2},qvec::Vector,psi::Vector{ComplexF64},lattice_params::Dict; kwargs...)
    if_periodic_phys = lattice_params["if_periodic_x"]
    if_periodic_synth = lattice_params["if_periodic_y"]
    phys_len = lattice_params["Lx"]
    synth_len = lattice_params["Ly"]

	occs = get_occupancy(psi,lattice_params; if_plot=false,densmat=rho)

	struc_fact = 0.0
	for j in 1:phys_len
		for s in 1:synth_len
			p1 = (j,s)
			p1_linear = linear_index(p1,phys_len,synth_len)
			for jj in 1:phys_len
				for ss in 1:synth_len
					p2 = (jj,ss)
					p2_linear = linear_index(p2,phys_len,synth_len)
					dist = find_dist(p1, p2, (phys_len,synth_len), (if_periodic_phys,if_periodic_synth))[2]
					struc_fact += occs[p1[1],p1[2]] * occs[p2[1],p2[2]] * exp(im * dot(qvec,dist))
				end
			end
		end
	end
	return struc_fact / sum(occs)
end

function cdw_sf(rho::Array{ComplexF64,2},psi::Vector{ComplexF64},lattice_params::Dict,qend=(0.0,3.0),howmany=50; kwargs...)
	if_plot = get(kwargs, :if_plot, true)

    if qend[1] == 0.0
        direction = 2
        axis_dir = "Y"
        qmax = qend[2]
    elseif qend[2] == 0.0
        direction = 1
        axis_dir = "X"
        qmax = qend[1]
    else
        direction = 3
    end

    if direction < 3
        qvec = zeros(Float64,2)
        qs = range(0.0,stop=qmax,length=howmany)
        struct_factor = zeros(ComplexF64,howmany)
        for (i,qx) in enumerate(qs)
            qvec[direction] = qx
            struct_factor[i] = cdw_structure_factor(rho,qvec,psi,lattice_params; kwargs...)
        end

        if if_plot
            plot_title = get(kwargs, :plot_title, "")
            plot_label = get(kwargs, :plot_label, nothing)
            isnothing(plot_label) ? fig = figure() : nothing
            isnothing(plot_label) ? plot(qs,abs.(struct_factor)) : plot(qs,abs.(struct_factor),label=plot_label)
            xlabel("q"*axis_dir)
            ylabel("Structure Factor")
            title("CDW Structure Factor, " * plot_title)
            legend()
        end
    else
        struct_factor = zeros(ComplexF64,howmany,howmany)
        qxs = range(-qend[1],stop=qend[1],length=howmany)
        qys = range(-qend[2],stop=qend[2],length=howmany)
        for (i,qx) in enumerate(qxs)
            for (j,qy) in enumerate(qys)
                qvec = [qx,qy]
                struct_factor[i,j] = cdw_structure_factor(rho,qvec,psi,lattice_params; kwargs...)
            end
        end

        if if_plot
            plot_title = get(kwargs, :plot_title, "")
            fig = figure()
            imshow(abs.(struct_factor),extent=[qxs[1],qxs[end],qys[1],qys[end]])
            colorbar()
            xlabel("qx")
            ylabel("qy")
            title("CDW Structure Factor, " * plot_title)
        end
    end

	return struct_factor
end

function spin_stiffness(energies::Vector,thetas::Vector; kwargs...)
    if_plot = get(kwargs, :if_plot, true)

    stiffnesses = zeros(Float64,length(thetas)-2)
    for i in 2:length(thetas)-1
        stiffnesses[i-1] = (energies[i+1] - 2*energies[i] + energies[i-1]) / (thetas[i+1] - thetas[i])^2
    end

    if if_plot
        fig = figure()
        plot(thetas[2:end-1],stiffnesses,"-p")
        xlabel("Theta")
        ylabel("Spin Stiffness")
        title("Spin Stiffness")
    end

    return stiffnesses
end

function get_normal_model_params_ed(params_dict::Dict)

    # set lattice parameters
    Lx = get(params_dict, "Lx", 4)
    Ly = get(params_dict, "Ly", Lx)
    N = get(params_dict, "N", 2)
    if_periodic_x = get(params_dict, "if_periodic_x", false)
    if_periodic_y = get(params_dict, "if_periodic_y", false)
    twist_angle = [get(params_dict, "tw1", 0.0),get(params_dict, "tw2", 0.0)]
    #=if typeof(twist_angle) == String
        println("Twist angle string is $twist_angle")
        tw_str = split(twist_angle,"c")
        tws = tryparse.(Float64,tw_str)
        println("Twist angles found are $tws")
        for (idx,twa) in enumerate(tws)
            if isnothing(twa)
                parts = parse.(Float64,split(tw_str[idx],"p"))
                tws[idx] = sum(parts .* [10.0^(-(i-1)) for i in 1:length(parts)])
            end
        end
        println("Twist angles are $tws")
        twist_angle = tws
    end=#
    expected_dimHilb = binomial(Lx*Ly,N)


    # set running operation parameters
    nev = get(params_dict,"nev",1)
    if_save_data = get(params_dict, "if_save_data", true)
    if if_periodic_x && if_periodic_y
		dataloc = get_folder_location("cluster-data/exact-diag/torus")
	elseif if_periodic_x || if_periodic_y
		dataloc = get_folder_location("cluster-data/exact-diag")
	elseif !if_periodic_x && !if_periodic_y
		dataloc = get_folder_location("cluster-data/exact-diag/obc")
	end
    dataloc = get(params_dict, "dataloc", dataloc)
    if occursin("geraghty1",dataloc)
        basis_dataloc = "/p/project/netenesyquma/geraghty1/data/data-ed/basis-files"
    else
        basis_dataloc = get_folder_location("cluster-data/exact-diag")#dataloc
    end
    opl = get(params_dict, "output_level", 1)
    if_exact = get(params_dict, "if_exact", false)
    if_densmat = get(params_dict, "if_densmat", false)
    if_find_data = get(params_dict, "if_find_data", true)
    if_function = get(params_dict, "if_function", false)
    running_args = (nev=nev,
                    if_exact=if_exact,
                    if_function=if_function,
                    if_densmat=if_densmat,
                    if_find_data=if_find_data,
                    if_save_data=if_save_data,
                    dataloc=dataloc,
                    basis_dataloc=basis_dataloc,
                    output_level=opl)


    opl > 0 ? println("Using ",N," particles with density ",round(N/(Lx*Ly),digits=3)) : nothing

    # build lattice parameters dictionary
    lattice_params::Dict{String,Any} = Dict("Lx"=>Lx,
                        "Ly"=>Ly,
                        "N"=>N,
                        "if_periodic_x"=>if_periodic_x,
                        "if_periodic_y"=>if_periodic_y,
                        "twist_angle"=>twist_angle)

    #println("Finished Building Lattice Parameters")
    # build long range interaction parameters
    stren = get(params_dict,"interaction_strength", 0.0)
    lr_dist = get(params_dict,"lr", "all")
    lr_dist == "all" ? lr_dist = Ly-1 : nothing
    scaling_type = get(params_dict,"scaling_type","flat")
    other_params_dict::Dict{String,Any} = Dict([("scaling",scaling_type)])
    if scaling_type != "flat"
        corr_length = get(params_dict,"corr_length",Ly)
        sigma = get(params_dict, "sigma", 1.0)
        blockade_radius = get(params_dict, "blockade_radius", 1.0)
        other_params_dict["corr_length"] = corr_length
        other_params_dict["sigma"] = sigma
        other_params_dict["blockade_radius"] = blockade_radius
    end
    us = long_range_scaling(lr_dist,Ly,stren; dict_to_symbols(other_params_dict)...)

    # get hopping anisotropy values
    hopping_anisotropy = get(params_dict,"hopping_anisotropy",1.0)
    if hopping_anisotropy < 1.0
		ty = 1.0 / hopping_anisotropy
		tx = 1.0
	else
		tx = 1.0 * hopping_anisotropy
		ty = 1.0
	end
    #=println("Using Alberto's Hopping Anisotropy")
    tx = hopping_anisotropy
    ty = 1/hopping_anisotropy=#
    
    # build magnetic field parameters
    alpha = get(params_dict,"alpha",nothing)
    if_bc_shift = get(params_dict,"if_bc_shift",true)
    x_shift,y_shift = if_bc_shift ? (!if_periodic_x, !if_periodic_y) : (0.0,0.0)
    if isnothing(alpha)
        filling = get(params_dict,"filling",0.5)
        alpha = N / (filling * (Lx - x_shift) * (Ly - y_shift))
        filling == 0.0 ? alpha = 0.0 : nothing
    end

    # build hamiltonian parameters dictionary and check fluxes for periodicity
    int_cutoff = get(params_dict,"interaction_cutoff",1e-5)
    which_dir = get(params_dict,"which_dir","virt")
    flux_dir = get(params_dict,"flux_direction","x")
    if if_periodic_y && !if_periodic_x
        flux_dir = "y"
    elseif !if_periodic_y && if_periodic_x
        flux_dir = "x"
    end
    if_check_fluxes = get(params_dict,"if_check_fluxes",true)
    if_check_fluxes ? flux_dir = check_fluxes(alpha,Lx,Ly,if_periodic_x,if_periodic_y,flux_dir) : nothing

    disorder_strength = get(params_dict,"disorder_strength",0.0)
    hamilt_params = Dict("alpha"=>alpha,
                        "flux_direction"=>flux_dir,
                        "tx"=>tx,
                        "ty"=>ty,
                        "hopping_anisotropy"=>hopping_anisotropy,
                        "disorder_strength"=>disorder_strength,
                        "U"=>us,
                        "which_dir"=>which_dir,
                        "interaction_cutoff"=>int_cutoff)

    println("Finished Building Model")
    return lattice_params,hamilt_params,running_args

end

function bin_values(vector, num_bins)
    min_val = minimum(vector)
    max_val = maximum(vector)
    bin_width = (max_val - min_val) / num_bins

    bins = Dict()
    for val in vector
        bin_index = min(floor(Int, (val - min_val) / bin_width) + 1, num_bins)
        bins[val] = bin_index
    end

    return bins
end

function charge_polarization(psi::Vector{ComplexF64},lattice_params::Dict; kwargs...)
    occs = get_occupancy(psi,lattice_params; if_plot=false)
    Ly = lattice_params["Ly"]

    cppul = sum([sum(occs[m,:] .* m) for m in 1:Ly])/(Ly) # charge polarization per unit length
    return cppul
end

function bulk_density(psi::Vector{ComplexF64},lattice_params::Dict,bulk_width_phys=1,bulk_width_synth=1; kwargs...)
    occ_mat = get_occupancy(psi,lattice_params; if_plot=false)
    #size(occ_mat)[1] == size(occ_mat)[2] ? bulk_width_virt = bulk_width_phys : nothing
	bulk_occ_mat = occ_mat[1+bulk_width_phys:end-bulk_width_phys,1+bulk_width_virt:end-bulk_width_virt]
	bulk_density = sum(bulk_occ_mat)/prod(size(bulk_occ_mat))
	return bulk_density
end

function twopointcorrelator(densmat::Array{ComplexF64,2},lattice_params::Dict; kwargs...)
    !isodd(lattice_params["Lx"]) && !isodd(lattice_params["Ly"]) ? error("Must have odd sized lattice") : nothing
    lattice_params["Lx"] != lattice_params["Ly"] ? error("Not implemented for non-square lattice") : nothing

    if_plot = get(kwargs,:if_plot,true)

    center_site = (Int((lattice_params["Lx"]+1)/2),Int((lattice_params["Ly"]+1)/2))
    center_site_linear = linear_index(center_site,lattice_params["Lx"],lattice_params["Ly"])

    twopointcorrs = zeros(Float64,lattice_params["Lx"],lattice_params["Ly"])
    for j in 1:lattice_params["Lx"]
        for s in 1:lattice_params["Ly"]
            next_site = (j,s)
            next_site_linear = linear_index(next_site,lattice_params["Lx"],lattice_params["Ly"])
            twopointcorrs[next_site[1],next_site[2]] = abs2(densmat[center_site_linear,next_site_linear])
        end
    end
    twopointcorrs ./= twopointcorrs[center_site[1],center_site[2]]

    if_plot ? plot_twopointcorrelator(twopointcorrs; kwargs...) : nothing

    return twopointcorrs
end

function plot_twopointcorrelator(twopointcorrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    imshow(twopointcorrs)
    colorbar()
    title("Two Point Correlator"*plot_title)
    xlabel("Physical")
    ylabel("Synthetic")
    return nothing
end

function pointfourhopping(psi::Vector{ComplexF64},s1::Int64,s2::Int64,lattice_params::Dict)
    hopmat = buildHopping(lattice_params,s1,s2)

    if s1 != s2
        return conj(transpose(psi)) * (conj(transpose(hopmat)) * hopmat) * psi
    else
        return conj(transpose(psi)) * (conj(transpose(hopmat)) * hopmat - hopmat) * psi
    end
end

function fourpointcorrelator(psi::Vector{ComplexF64},lattice_params::Dict; kwargs...)
    !isodd(lattice_params["Lx"]) && !isodd(lattice_params["Ly"]) ? error("Must have odd sized lattice") : nothing
    lattice_params["Lx"] != lattice_params["Ly"] ? error("Not implemented for non-square lattice") : nothing

    if_plot = get(kwargs,:if_plot,true)

    center_site = (Int((lattice_params["Lx"]+1)/2),Int((lattice_params["Ly"]+1)/2))
    center_site_linear = linear_index(center_site,lattice_params["Lx"],lattice_params["Ly"])

    fourpointcorrs = zeros(Float64,lattice_params["Lx"],lattice_params["Ly"])
    for j in 1:lattice_params["Lx"]
        for s in 1:lattice_params["Ly"]
            next_site = (j,s)
            next_site_linear = linear_index(next_site,lattice_params["Lx"],lattice_params["Ly"])
            fourpointcorrs[next_site[1],next_site[2]] = abs2(pointfourhopping(psi,center_site_linear,next_site_linear,lattice_params))
        end
    end

    #fourpointcorrs ./= fourpointcorrs[center_site[1],center_site[2]]

    if_plot ? plot_fourpointcorrelator(fourpointcorrs; kwargs...) : nothing

    return fourpointcorrs
end

function plot_fourpointcorrelator(fourpointcorrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    imshow(fourpointcorrs)
    colorbar()
    title("Four Point Correlator" * plot_title)
    xlabel("Physical")
    ylabel("Synthetic")
    return nothing
end

if false
    lx = 4
    ly = 8
    N = 4

    cols = ["b","r","g","m","c"]
    if 10 > length(cols)
        cols = repeat(cols,ceil(Int,10/length(cols)))
    end

    nrgs1 = []
    nrgs2 = []
    nrgs3 = []
    nrgs4 = []
    xs = []
    ys = []
    
    params_dict = Dict([("Lx",lx),("Ly",ly),("N",N),("if_periodic_x",true),("if_periodic_y",true)])
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)

        for f in files
            data,metadata = read_data_jld2(f,dataloc; output_level=0)
            nrgs = data["nrg"]
            intstren = metadata["U"][end]
            anis = metadata["hopping_anisotropy"]

            if anis > 1.3
                continue
            end

            append!(xs,[anis])
            append!(ys,[intstren])
            #append!(nrgs1,[nrgs[2] - nrgs[1]])
            append!(nrgs2,[nrgs[3] - nrgs[2]])
            #append!(nrgs3,[nrgs[4] - nrgs[1]])
            #append!(nrgs4,[nrgs[5] - nrgs[1]])

        end
    
    bin_count = 100
    data_dict = bin_values(nrgs2,bin_count)
    bv = [data_dict[val] for val in nrgs2]
    min_nrgs2, max_nrgs2 = minimum(nrgs2), maximum(nrgs2)
    normalized_bv = [(val - minimum(bv)) / (maximum(bv) - minimum(bv)) * (max_nrgs2 - min_nrgs2) + min_nrgs2 for val in bv]

    #= now try to find linear fit to maximum of nrgs2 at each hopping anis row
    hanises = unique(xs)
    #max_intstrens = zeros(Float64,length(hanises))
    min_intstrens = zeros(Float64,length(hanises))
    for (j,hanis) in enumerate(hanises)
        indices = findall(x -> x == hanis, xs)
        if length(indices) == 1
            #println("Only one data point for Anis=$hanis")
            continue
        end
        #max_index = findfirst(x -> x == maximum(nrgs2[i] for i in indices), nrgs2)
        #max_intstrens[j] = ys[max_index]
        relevant_nrgs = [nrgs2[i] for i in indices]
        min_index = findfirst(i -> isapprox(relevant_nrgs[i],0.0,atol=1e-6),1:length(relevant_nrgs))
        println("Found Level Crossing for Anis=$hanis at index $min_index")
        println("It has value $(relevant_nrgs[min_index])")
        println("This is intstrength of $(ys[indices[min_index]]) \n")
        min_intstrens[j] = ys[indices[min_index]]
    end

    linfit(x,p) = p[1] .* x .+ p[2]
    fit = curve_fit(linfit,min_intstrens,hanises,[1.0,0.0])
    fit_xs = range(minimum(ys),2.0,length=10)
    plot(fit_xs,linfit(fit_xs,fit.param),c="r",label="m=$(round(fit.param[1],digits=3))")=#

    fig = figure()
    scatter(ys, xs, c=normalized_bv, cmap="viridis")
    colorbar()
    
    ylim(minimum(xs)-0.05,0.05+maximum(xs))

    ylabel("Hopping Anisotropy")
    xlabel("Interaction Strength")
    legend()
    title("Energy Gap btw 2nd and 3rd 4x8 N=4")

end

if false
    lx = 5
    N = 5

    cols = ["b","r","g","m","c"]
    if 10 > length(cols)
        cols = repeat(cols,ceil(Int,10/length(cols)))
    end
    
    anises = [1/0.6,1.2,1.0,0.8,0.6]

    nrgs1 = []
    nrgs2 = []
    nrgs3 = []
    nrgs4 = []
    xs = []
    ys = []
    #for anis in anises
        params_dict = Dict([("Lx",lx),("Ly",lx),("N",N),("if_periodic_x",true),("if_periodic_y",true)])#,("hopping_anisotropy",anis)])
        dataloc = get_folder_location("cluster-data/exact-diag/torus")
        files = find_data_file(params_dict,"ed",dataloc; output_level=0)

        #fig = figure()
        for f in files
            data,metadata = read_data_jld2(f,dataloc; output_level=0)
            nrgs = data["nrg"]
            intstren = metadata["U"][end]
            anis = metadata["hopping_anisotropy"]

            if intstren > 6.0
                continue
            elseif anis > 2.0 || anis <= 1.0
                continue
            end

            append!(xs,[anis])
            append!(ys,[intstren])
            #append!(nrgs1,[nrgs[2] - nrgs[1]])
            append!(nrgs2,[nrgs[3] - nrgs[2]])
            #append!(nrgs3,[nrgs[4] - nrgs[1]])
            #append!(nrgs4,[nrgs[5] - nrgs[1]])
            

            #=cols = ["b","r","g","m","c"]
            if length(nrgs) > length(cols)
                cols = repeat(cols,ceil(Int,length(nrgs)/length(cols)))
            end

            for i in 1:length(nrgs)
                scatter(intstren,nrgs[i] - nrgs[1],c=cols[i])
            end
            xlabel("Interaction Strength")
            ylabel("Energy - E1")
            title("Energy Spectrum for Hopping Anisotropy = $(round(anis,digits=4))")=#
        end
    #end
    #=plot3D(xs,ys,nrgs1,"p")
    plot3D(xs,ys,nrgs2,"p")
    plot3D(xs,ys,nrgs3,"p")
    plot3D(xs,ys,nrgs4,"p")
    xlabel("Hopping Anisotropy")
    ylabel("Interaction Strength")
    zlabel("Energy Difference")=#

    bin_count = 100
    data_dict = bin_values(nrgs2,bin_count)
    bv = [data_dict[val] for val in nrgs2]
    min_nrgs2, max_nrgs2 = minimum(nrgs2), maximum(nrgs2)
    normalized_bv = [(val - minimum(bv)) / (maximum(bv) - minimum(bv)) * (max_nrgs2 - min_nrgs2) + min_nrgs2 for val in bv]

    # now try to find linear fit to maximum of nrgs2 at each hopping anis row
    hanises = unique(xs)
    max_intstrens = zeros(Float64,length(hanises))
    for (j,hanis) in enumerate(hanises)
        indices = findall(x -> x == hanis, xs)
        max_index = findfirst(x -> x == maximum(nrgs2[i] for i in indices), nrgs2)
        max_intstrens[j] = ys[max_index]
    end

    #=linfit(x,p) = p[1] .* x .+ p[2]
    fit = curve_fit(linfit,max_intstrens,hanises,[1.0,0.0])
    fit_xs = range(minimum(ys),2.0,length=10)
    plot(fit_xs,linfit(fit_xs,fit.param),c="r",label="m=$(round(fit.param[1],digits=3))")=#

    fig = figure()
    scatter(ys, xs, c=normalized_bv, cmap="viridis")
    colorbar()
    
    ylim(minimum(xs)-0.05,0.05+maximum(xs))

    ylabel("Hopping Anisotropy")
    xlabel("Interaction Strength")
    legend()
    title("Energy Gap with linear fit at maximum")

end

if false
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    params_dict = Dict([("if_periodic_x",true)])
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    for f in files
        data,metadata = read_data_jld2(f,dataloc)
        if_func = metadata["if_function"]
        lattice_params = get_lattice_params_from_metadata(metadata)
        dimHilb = binomial(lattice_params["Lx"]*lattice_params["Ly"],lattice_params["N"])

        col = if_func ? "r" : "b"
        scatter(dimHilb,metadata["runtime"],c=col)
    end
    xlabel("Hilbert Space Dimension")
    ylabel("Runtime")
    xscale("log")
    yscale("log")
end

if false
    params_dict = Dict([("Lx",5),("Ly",8),("N",4),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)


    cols = ["b","g","r","m","c"]
    if 10 > length(cols)
        cols = repeat(cols,ceil(Int,10/length(cols)))
    end

    for f in files
        d,m = read_data_jld2(f,dataloc)
        intstren = m["U"][end]

        for i in 1:length(d["nrg"])
            scatter(intstren,d["nrg"][i] - d["nrg"][1],c=cols[i])
        end
    end
    xlabel("Interaction Strength")
    ylabel("Energy")

end
        
    

if true
    #fig = figure()
    #xlabel("Hopping Anisotropy")
    #ylabel("Gap")
    #lx = 6
    #n = 3
    #for (idx,n) in enumerate([2,3,4,5])
    #intstrens = range(0.0,3.0,length=20)
    #other_intstrens = range(2.0,10.0,length=37)
    #intstrens = sort([intstrens; other_intstrens])
    #change = 0.001
    #real_alphas = [range(0.1,0.21,length=10); range(0.22,0.28,length=10); range(0.29,0.35,length=5)]
    #howmany = length(real_alphas)
    #alphas = [real_alphas; real_alphas .+ change]
    #all_bds = zeros(Float64,length(alphas))
    #thetas = range(0.01,0.5,length=50)
    #all_nrgs = zeros(Float64,length(thetas))
    #anises = range(1.0,5.0,length=20)
    #nus = range(0.4,0.6,length=100)
    #alphas = range(0.16,2*3/(6*5),length=30)
    #for (idx,alpha) in enumerate(alphas)
    #for (idx,ly) in enumerate([5,6,7,8])
    #for (idx,nu) in enumerate(nus)
    #for (idx,anis) in enumerate(anises)
    #sigmas = vcat([1/i for i in 1:5],[i for i in 2:5])#vcat(range(1.0,5.0,length=5),[100.0])
    #for (idx,intstren) in enumerate(intstrens)
    #for (idx2,sigma) in enumerate(sigmas)
    #for lrd in [0,1]
    #tws = range(0.0,1.0,length=20)
    #cps = zeros(Float64,length(tws))
    #for (idx,tw1) in enumerate(tws)
    #for (idx2,tw2) in enumerate(tws)
    #for tw1 in tws
    #for ii in 1:1
        #change_nrgs = zeros(Float64,3)
        #for (ii,change) in enumerate([0,0.0001,0.0002])
        change = 0.0
        params_dict = Dict([("Lx",4),("Ly",4),("N",2),("tw1",0.0),("tw2",0.0),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",10000.0),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
        #params_dict = make_args_dict(ARGS)

        # set number of open cores
        open_cores = get(params_dict, "open_cores", 5)
        if typeof(open_cores) != String
            BLAS.set_num_threads(open_cores)
            display(BLAS.get_config())
        end

        lattice_params,hamilt_params,running_args = get_normal_model_params_ed(params_dict)    
        basis_dataloc = running_args.basis_dataloc

        # build filename dictionary
        filename_dict = make_filename_dict(lattice_params,hamilt_params)
        #display(filename_dict)
        if_exists,found_data = running_args.if_find_data ? check_data_exists(filename_dict,"ed"; location=running_args.dataloc,output_level=false) : (false,nothing)

        # some old data has bad naming with int_stren = 1.0 even though rest of Us is zeros
        if params_dict["interaction_strength"] == 1.0 && if_exists
            if found_data[2]["U"] != hamilt_params["U"]
                println("Found Data has wrong Interaction Potential")
                if_exists = false
            end
        end

        # check if data exists and rerun if need more eigenstates
        if if_exists
            start_time = time()
            println("Found existing data: ",found_data[2]["filename"])
            if running_args.nev > found_data[2]["nev"]
                running_args.output_level > 0 ? println("Asking for more eigenstates than in file, rerunning") : nothing
                start_time = time()
                full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
                running_args.output_level > 0 ? println("Made basis in ",time()-start_time) : nothing
                lattice_params["full_basis"] = full_basis 
                states,nrgs,rhos = rerun_eigenstates(running_args.nev,lattice_params,hamilt_params,found_data[2],found_data[1]; running_args...)
            elseif running_args.nev < found_data[2]["nev"]
                running_args.output_level > 0 ? println("Asking for fewer eigenstates than in file, using existing data") : nothing
                states = found_data[1]["state"][1:running_args.nev]
                nrgs = found_data[1]["nrg"][1:running_args.nev]
                #rhos = found_data[1]["densmat"][1:running_args.nev]
                full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
                running_args.output_level > 0 ? println("Made basis in ",time()-start_time) : nothing
                lattice_params["full_basis"] = full_basis 
            else
                states = found_data[1]["state"]
                nrgs = found_data[1]["nrg"]
                rhos = found_data[1]["densmat"]
                full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
                running_args.output_level > 0 ? println("Made basis in ",time()-start_time) : nothing
                lattice_params["full_basis"] = full_basis 
            end
        else

            # make basis only if data doesn't exist
            start_time = time()
            full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
            running_args.output_level > 0 ? println("Made basis in ",time()-start_time) : nothing
            lattice_params["full_basis"] = full_basis 

            # run exact diagonalization for find eigenstates
            if running_args.nev == 1
                states,nrgs,rhos = find_ground_state(lattice_params,hamilt_params; running_args...)
            else
                if running_args.if_function
                    states,nrgs,rhos = find_eigenstates(running_args.nev,lattice_params,hamilt_params; running_args...)
                else
                    states,nrgs,rhos,hh = find_eigenstates(running_args.nev,lattice_params,hamilt_params; running_args...)
                end
            end
        end
        #change_nrgs[ii] = nrgs[2] - nrgs[1]
    #end

        #if idx == 1
        #    get_occupancy(states[1],lattice_params; plot_title=" LR=$(params_dict["interaction_strength"])")
        #    fig = figure()
        #end

        #=rho1 = density_matrix(states[1],lattice_params)
        rho2 = density_matrix(states[2],lattice_params)
        rez21 = twopointcorrelator(rho1,lattice_params; if_plot=false)
        rez22 = twopointcorrelator(rho2,lattice_params; if_plot=false)
        plot_twopointcorrelator((rez21+rez22) ./ 2; plot_title=" LR=$(params_dict["interaction_strength"]), Anis=$(params_dict["hopping_anisotropy"])")
        rez41 = fourpointcorrelator(states[1],lattice_params; if_plot=false)
        rez42 = fourpointcorrelator(states[2],lattice_params; if_plot=false)
        plot_fourpointcorrelator((rez41+rez42) ./ 2; plot_title=" LR=$(params_dict["interaction_strength"]), Anis=$(params_dict["hopping_anisotropy"])")=#

        #rho1 = density_matrix(states[1],lattice_params; output_level=0)
        #rez21 = twopointcorrelator(rho1,lattice_params; if_plot=true)
        #isotropicness = sum(abs.(rez21 .- transpose(rez21))) / sum(rez21)
        #title("2PT LR=$(params_dict["interaction_strength"]) Sigma=$(params_dict["sigma"]), Anis=$(params_dict["hopping_anisotropy"]), Iso=$(round(isotropicness,digits=3))")

        #rez41 = fourpointcorrelator(states[1],lattice_params; plot_title=" LR=$(params_dict["interaction_strength"]) $(params_dict["scaling_type"]), Anis=$(params_dict["hopping_anisotropy"])")
        
        cols = ["b","g","r","m","c"]
        if running_args.nev > length(cols)
            cols = repeat(cols,ceil(Int,running_args.nev/length(cols)))
        end#

        #=xxs = intstrens
        if idx == 1
            #if idx2 == length(sigmas)
            #    scatter(xxs[idx],isotropicness,c=cols[idx2],label="Flat") 
            #else
                scatter(xxs[idx],isotropicness,c=cols[idx2],label="$(round(sigma,digits=2))")
            #end
        else
            scatter(xxs[idx],isotropicness,c=cols[idx2])    
        end
        xlabel("Interaction Strength")
        ylabel("Anisotropicness")
        legend()=#
        #=scatter(sigma,isotropicness,c="b")
        xlabel("Sigma")
        ylabel("Anisotropicness")
        yscale("log")
        xscale("log")
        title("2pt Anisotropicness vs Sigma")=#

        #=nrgs_derivative = (change_nrgs[3] - change_nrgs[1]) / 0.0002
        if idx == 1
            scatter(intstren,nrgs_derivative,c="b",label="Deriv")
        else
            scatter(intstren,nrgs_derivative,c="b")
        end

        nrgs_2nd_derivative = (change_nrgs[3] - 2*change_nrgs[2] + change_nrgs[1]) / 0.0001^2
        if idx == 1
            scatter(intstren,nrgs_2nd_derivative,c="r",label="2nd Deriv")
        else
            scatter(intstren,nrgs_2nd_derivative,c="r")
        end=#

        #=xxs = intstrens
        for i in 1:running_args.nev
            change = abs(xxs[1] - xxs[2])
            xval = xxs[idx]
            shift = (i - running_args.nev/2) * ((0.1*change)/(running_args.nev/2))
            scatter(xval + shift,nrgs[i]-nrgs[1],c=cols[i])
        end
        #legend()
        #xlabel("System Size")
        xlabel("Interaction Strength")
        #xlabel("Flux")
        #xlabel("Theta_x / 2pi")
        #ylabel("Theta_y")
        ylabel("NRG")=#
        #xlabel("Hopping Anisotropy tx/ty")
        #title("4x4 N=2, Anis=$(hamilt_params["hopping_anisotropy"])")
        #title("Topological Degeneracy Closing in Thermodynamic Limit")#
        #title("Spectrum Twist BC $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"]) Anis=$(params_dict["hopping_anisotropy"])")

        #=for i in 1:running_args.nev
            scatter3D(tw1,tw2,nrgs[i],c=cols[i])
        end
        xlabel("Theta1 / 2pi")
        ylabel("Theta2 / 2pi")
        zlabel("NRG")=#

        #=cps[idx] = charge_polarization(states[1],lattice_params)

        scatter(tw1,cps[idx],c="b")
        xlabel("Theta")
        ylabel("Charges Pumped")
        title("Synth Length = $(params_dict["Ly"])")=#

        #=xxs = tws
        if idx == length(xxs)
            get_occupancy(states[1],lattice_params; plot_title="$(xxs[end])")
        end=#


    #end
    #end

    #=plot3D(xs,ys,nrgs1 .- nrgs1,label="E1")
    plot3D(xs,ys,nrgs2 .- nrgs1,label="E2")
    plot3D(xs,ys,nrgs3 .- nrgs1,label="E3")
    xlabel("Theta1 / 2pi")
    ylabel("Theta2 / 2pi")
    zlabel("NRG")=#


    #bdderivs = (all_bds[howmany+1:end] .- all_bds[1:howmany]) ./ change
    #fillings = n ./ (alphas[1:howmany] .* ((lx-1)*(lx-1)))
    #fig = figure()
    #scatter(fillings,bdderivs)
    #xlabel("Filling")
    #title("Derivative of Bulk Density")

    #th_alphas = range(minimum(alphas),maximum(alphas),length=100)
    #plot(th_alphas,4*(sin.(pi .* th_alphas)).^2,label="Theory",c="b")
    #legend()

    #plot(alphas,4 .* (sin.(pi .* alphas)).^2,c="r",label="Theory")
    #legend()


    #=maybefit = curve_fit(hh_gap_fit,anises[end-20:end],gapvals[end-20:end],[1.0,0.0])
    fig = figure()
    plot(anises,hh_gap_fit(anises,maybefit.param) ./ maybefit.param[1],c="r",label="Fit")
    legend()
    title("$lx x $lx N=$n HH Gap Exact: U_eff = $(round(maybefit.param[1],digits=2))")
    #ylim(0.0,1.5)

    scatter(anises,gapvals ./ maybefit.param[1],c="b")

    xlabel("Hopping Anisotropy tx/ty")
    ylabel("E2 - E0")=#

end





































"fin"
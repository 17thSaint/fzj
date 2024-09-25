#####################################################
#=

This file is for building the ED Hamiltonian and then using it to calculate the eigenstates

Depends on:
    other-funcs/data-storage-funcs.jl

=#
######################################################



#using Pkg
#Pkg.activate(".")
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

# pretty plot of N=4 4x8 and 8x4 on torus for poster
if false
    params_dict = Dict([("Lx",4),("Ly",8),("N",4),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    cols = ["b","r","k","k","k","k","k","k","k","k"]

    nrgs::Dict{String,Vector} = Dict([("intstrens",[])])
    etop::Int = params_dict["Lx"] == 8 ? 10 : 5
    for i in 1:etop
        nrgs["E$(i-1)"] = []
    end
    
    for f in files
        d,m = read_data_jld2(f,dataloc)
        intstren = m["U"][end]
        if intstren > 2.0
            continue
        end
        append!(nrgs["intstrens"],intstren)
        for i in 1:etop
            append!(nrgs["E$(i-1)"],d["nrg"][i] - d["nrg"][1])
        end
    end
    
    for i in 1:etop
        if i > 3
            scatter(nrgs["intstrens"],nrgs["E$(i-1)"],c=cols[i])
        elseif i == 3
            scatter(nrgs["intstrens"],nrgs["E$(i-1)"],c=cols[i],label="E2-$(etop-1)")
        else
            scatter(nrgs["intstrens"],nrgs["E$(i-1)"],c=cols[i],label="E$(i-1)")
        end
    end
    legend()
    ylim([-0.01,0.4])
    xlabel("Long Range Interaction Strength V/t")
    ylabel("E - E0")
    title("Energy Spectrum N=4 $(params_dict["Lx"])x$(params_dict["Ly"])")
end

# quick plot energies of N=4 4x8 and 8x4 on torus
if false
    params_dict = Dict([("Lx",4),("Ly",8),("N",4),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    files = find_data_file(params_dict,"ed",dataloc; output_level=0)


    cols = ["b","g","r","m","c"]
    if 10 > length(cols)
        cols = repeat(cols,ceil(Int,10/length(cols)))
    end

    
    for f in files
        d,m = read_data_jld2(f,dataloc)
        intstren = m["U"][end]
        if intstren > 2.0
            continue
        end
        for i in 1:length(d["nrg"])
            scatter(intstren,d["nrg"][i] - d["nrg"][1],c=cols[i])
        end
    end
    xlabel("Interaction Strength")
    ylabel("Energy")

end
       
#= finding input parameters to do finite size scaling in constant alpha
all_sizes = collect(4:20)
good_alphas = []
for a in 1:50
    for b in 1:50
        if a/b < 0.4
            if_works = isinteger.((a/b) .* all_sizes)
            working_vals = all_sizes[if_works]
            particles = (0.5*a/b) .* (working_vals .^2)
            dilute_enough = (particles ./ working_vals) .!= 1.0
            mostly_good = working_vals[dilute_enough]
            if length(mostly_good) > 3
                if !(a/b in good_alphas)
                    append!(good_alphas,[a/b])
                    println("Found alpha = $a / $b with working values $mostly_good and particles $(particles[dilute_enough])")
                end
            end
        end
    end
end=#
    

# run data collection with for loops
if true
    #fig = figure()
    #xlabel("Hopping Anisotropy")
    #ylabel("Gap")
    #lx = 6
    #n = 3
    #for (idx,n) in enumerate([2,3,4,5])
    #intstrens = range(0.0,2.0,length=10)
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
    #for (idx,ly) in enumerate(lys)
    #for (idx,nu) in enumerate(nus)
    #for (idx,anis) in enumerate(anises)
    #sigmas = vcat([1/i for i in 1:5],[i for i in 2:5])#vcat(range(1.0,5.0,length=5),[100.0])
    #for (idx,intstren) in enumerate(intstrens)
    #for (idx2,sigma) in enumerate(sigmas)
    #for lrd in [0,1]

    #args_dict = Dict([("which_twist_angle",1)])
    lx,ly,n = 4,4,2

    tws = range(0.0,1.0,length=11)
    #omegas = zeros(ComplexF64,length(tws),length(tws))
    #gammas1::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws))
    #gammas2::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws))
    ref_multiplets,rm1_name,rm2_name = get_reference_multiplets(lx,ly,n)
    #cps = zeros(Float64,length(tws))    tws[args_dict["which_twist_angle"]]
    for (idx,tw1) in enumerate(tws)
    for (idx2,tw2) in enumerate(tws)
    #for tw1 in tws
    #for ii in 1:1
        params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("tw1",tw1),("tw2",tw2),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",0.0),("lr",0),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",false)])
        #params_dict = make_args_dict(ARGS)

        # set number of open cores
        open_cores = get(params_dict, "open_cores", 5)
        if typeof(open_cores) != String
            BLAS.set_num_threads(open_cores)
            display(BLAS.get_config())
        end

        states,nrgs,rhos,filename = run_normal_ed(params_dict)

        cols = ["b","g","r","m","c"]
        if params_dict["nev"] > length(cols)
            cols = repeat(cols,ceil(Int,params_dict["nev"]/length(cols)))
        end#

        #=omega = get_omega(states[1],states[2],ref_multiplets)
        gamma1 = get_gamma(states[1],states[2],ref_multiplets[1:2])
        gamma2 = get_gamma(states[1],states[2],ref_multiplets[3:4])
        hatsugai_datadict = Dict([("omega",omega),("gamma1",gamma1),("gamma2",gamma2),("rm1_name",rm1_name),("rm2_name",rm2_name)])
        modify_data_jld2(hatsugai_datadict,running_args[:dataloc]*"/"*filename,"metadata")=#

        xxs = tws
        for i in 1:params_dict["nev"]
            change = abs(xxs[1] - xxs[2])
            xval = xxs[idx]
            shift = (i - params_dict["nev"]/2) * ((0.1*change)/(params_dict["nev"]/2))
            scatter(xval + shift,nrgs[i],c=cols[i])
        end
        #legend()
        #xlabel("System Size")
        #xlabel("Interaction Strength")
        #xlabel("Flux")
        xlabel("Theta_x / 2pi")
        #ylabel("Theta_y")
        ylabel("NRG")
        #yscale("log")
        #xlabel("Hopping Anisotropy tx/ty")
        #title("4x4 N=2, Anis=$(hamilt_params["hopping_anisotropy"])")
        #title("Topological Degeneracy Closing in Thermodynamic Limit")#
        #title("Spectrum Twist BC $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"]) Anis=$(params_dict["hopping_anisotropy"])")


    end
    end

end

#=fig = figure()
imshow(angle.(omegas))
colorbar()
title("Omega Phase")

fig = figure()
imshow(abs2.(omegas))
colorbar()
title("Omega Magnitude")

fig = figure()
imshow(abs2.(gammas1))
colorbar()
title("Gamma1 Magnitude")

fig = figure()
imshow(abs2.(gammas2))
colorbar()
title("Gamma2 Magnitude")=#





































"fin"
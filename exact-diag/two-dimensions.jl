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

            # apply magnetic flux but switch for correct cylinder flux counting
            coeff = -tx * exp(im*alpha*starting_site[2]*2*pi)
            if if_periodic_y && !if_periodic_x
                coeff = -tx
            end

            # apply twist at boundary
            if twist
                coeff *= exp(im*2*pi*twist_angle)
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
            # apply magnetic flux for correct cylinder flux counting
            if if_periodic_y && !if_periodic_x
                coeff *= exp(im*alpha*starting_site[1]*2*pi)
            end

            # apply twist at boundary
            if twist
                coeff *= exp(im*2*pi*twist_angle)
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
            site1 in full_basis[:,j] ? hop[j,j] = 1.0+0.0*im : nothing
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
        output_level > 0 ? println(round(i/(Lx*Ly)*100,digits=2),"% done.") : nothing
    end

    output_level > 0 ? println("Density Matrix: Elapsed time: ",time()-start_time) : nothing

    return rho
end

function find_eigenstates(nev::Int,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if_densmat = get(kwargs,:if_densmat,true)
    if_save_data = get(kwargs,:if_save_data,false)
    if_exact = get(kwargs,:if_exact,false)

    metadata_dict = merge(merge(lattice_params,hamilt_params),named_tuple_to_dict(kwargs))

    metadata_displaying::Dict = copy(metadata_dict)
    delete!(metadata_displaying,"full_basis")
    output_level > 0 ? display(metadata_displaying) : nothing
    
    start_time = time()
    H = buildHam(lattice_params,hamilt_params; output_level)
    dimHilb = size(H)[1]
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
    output_level > 0 ? println("Ground State: Elapsed time: ",time()-start_time) : nothing

    sorted_indices = sortperm(rez[1])
    states = rez[2][sorted_indices][1:nev]
    nrgs = rez[1][sorted_indices][1:nev]

    rhos = []
    if if_densmat
        for state in states
            append!(rhos,[density_matrix(state,lattice_params)])
        end
    end

    if_save_data ? save_eigenstates(states,rhos,nrgs,metadata_dict) : nothing

    return states,nrgs,rhos,H
end

function rerun_eigenstates(nev::Int,lattice_params::Dict,hamilt_params::Dict,metadata::Dict,data_dict::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if_densmat = get(kwargs,:if_densmat,true)
    if_save_data = get(kwargs,:if_save_data,false)

    metadata_displaying::Dict = copy(metadata)
    delete!(metadata_displaying,"full_basis")
    delete!(metadata_displaying,"H")
    output_level > 0 ? display(metadata_displaying) : nothing
    
    H = metadata["H"]
    output_level > 0 ? println("Sparsity = ",nnz(H)/size(H)[1]^2) : nothing

    start_time = time()
    previous_nev = metadata["nev"]
    rez = eigsolve(H,nev+3)
    output_level > 0 ? println("Ground State: Elapsed time: ",time()-start_time) : nothing

    sorted_indices = sortperm(rez[1])
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
    filename_dict = Dict([("Lx",metadata["Lx"]),("Ly",metadata["Ly"]),("N",metadata["N"]),("alpha",metadata["alpha"]),("hopping_anisotropy",metadata["hopping_anisotropy"]),("interaction_strength",metadata["U"][1]),("if_periodic_x",metadata["if_periodic_x"]),("if_periodic_y",metadata["if_periodic_y"])])
    filename = join(["ed",make_parameters_filename(filename_dict)],"-")
    metadata["filename"] = filename
    full_loc = join([dataloc,filename],"/")
    println("Filename: ",full_loc)
    write_data_jld2(full_loc,data_dict,metadata)
end

function physical_correlation(densmat::Array{ComplexF64,2},Lx::Int64,Ly::Int64; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    phys_corrs = Array{Float64,2}(undef,Lx-1,Ly)
    for j in 2:Lx
        for s in 1:Ly
            phys_corr = densmat[linear_index((j,s),Lx,Ly),linear_index((1,s),Lx,Ly)]
            phys_corr /= sqrt(densmat[linear_index((j,s),Lx,Ly),linear_index((j,s),Lx,Ly)] * densmat[linear_index((1,s),Lx,Ly),linear_index((1,s),Lx,Ly)])
            phys_corrs[j-1,s] = abs(phys_corr)
        end
    end

    if_plot ? plot_physical_correlation(phys_corrs; kwargs...) : nothing

    return phys_corrs
end

function synthetic_correlation(densmat::Array{ComplexF64,2},Lx::Int64,Ly::Int64; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    syn_corrs = Array{Float64,2}(undef,Lx,Ly-1)
    for j in 1:Lx
        for s in 2:Ly
            syn_corr = densmat[linear_index((j,s),Lx,Ly),linear_index((j,1),Lx,Ly)]
            syn_corr /= sqrt(densmat[linear_index((j,s),Lx,Ly),linear_index((j,s),Lx,Ly)] * densmat[linear_index((j,1),Lx,Ly),linear_index((j,1),Lx,Ly)])
            syn_corrs[j,s-1] = abs(syn_corr)
        end
    end

    if_plot ? plot_synthetic_correlation(syn_corrs; kwargs...) : nothing

    return syn_corrs
end

function plot_synthetic_correlation(syn_corrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(syn_corrs)[1]
        plot(1:size(syn_corrs)[2],syn_corrs[i,:],"-p",label="$i")
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
        plot(1:size(phys_corrs)[1],phys_corrs[:,i],"-p",label="$i")
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
    return Dict([("Lx",lattice_params["Lx"]),("Ly",lattice_params["Ly"]),("N",lattice_params["N"]),("alpha",hamilt_params["alpha"]),("hopping_anisotropy",hamilt_params["tx"]/hamilt_params["ty"]),("interaction_strength",hamilt_params["U"][1]),("if_periodic_x",lattice_params["if_periodic_x"]),("if_periodic_y",lattice_params["if_periodic_y"])])
end

function get_lattice_params_from_metadata(metadata::Dict)
    return Dict([("Lx",metadata["Lx"]),("Ly",metadata["Ly"]),("N",metadata["N"]),("if_periodic_x",metadata["if_periodic_x"]),("if_periodic_y",metadata["if_periodic_y"]),("full_basis",metadata["full_basis"])])
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

    return sqrt(dx^2 + dy^2),(dx,dy)
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

#=Lx = 4
Ly = 4
N = 2
lattice_params::Dict{String,Any} = Dict("Lx"=>Lx,
                        "Ly"=>Ly,
                        "N"=>N,
                        "if_periodic_x"=>true,
                        "if_periodic_y"=>true)
lattice_params["full_basis"] = n_particle_basis(N,Lx,Ly)
hamilt_params::Dict{String,Any} = Dict("tx"=>1.0,
                        "ty"=>0.0,
                        "alpha"=>0.0,
                        "U"=>[0.0],
                        "interaction_cutoff"=>0.0)
ttham = buildHam(lattice_params,hamilt_params)

for i in 1:Ly
    for j in i+1:Ly
        thisstate = zeros(ComplexF64,size(ttham,1))

        firstcolinds = linear_index((1,i),Lx,Ly) .+ [k for k in 0:Lx-1]
        for k in 1:Lx
            for k2 in 1:Lx
                if k == k2
                    continue
                end
                thisstate[find_basis_index(sort([firstcolinds[k],firstcolinds[k2]],rev=true))] += 1.0
            end
        end

        secondcolinds = linear_index((1,j),Lx,Ly) .+ [k for k in 0:Lx-1]
        for k in 1:Lx
            for k2 in 1:Lx
                if k == k2
                    continue
                end
                thisstate[find_basis_index(sort([secondcolinds[k],secondcolinds[k2]],rev=true))] += 1.0
            end
        end

        #=for k in 1:Lx
            for k2 in 1:Lx
                thisstate[find_basis_index(sort([firstcolinds[k],secondcolinds[k2]],rev=true))] += 1.0
            end
        end=#

        thisstate ./= sqrt(norm(thisstate))

        nrg = round(real((transpose(thisstate) * ttham * thisstate)),digits=3)
        get_occupancy(thisstate,lattice_params; plot_title="Energy = $nrg")
    end
end=#


#anises = range(1.0,100.0,length=30)
#gapvals = zeros(Float64,length(anises))
#lx = 6
#n = 4
#=ac = n / (0.5 * (lx-1)*(lx-1))
as = n / (0.4 * (lx-1)*(lx-1))
af = n / (0.6 * (lx-1)*(lx-1))
hm = 20=#
#alphas = range(0.25,1.25,length=40)#[1/(i*lx/10) for i in 1:10]#[1/(1*lx),1/(0.5*lx),1/(2*lx),1/(0.1*lx),1/(0.2*lx)]#[1/(i*lx) for i in 1:5]

#=for thisalpha in alphas
fig = figure()
ks = range(-pi,pi,length=1000)
realks = [2*pi*i/lx - pi for i in 0:lx]
for s in 1:4
    plot(ks,-2 .* cos.(2*pi*s*thisalpha .- ks),label="$s")
    scatter(realks,-2 .* cos.(2*pi*s*thisalpha .- realks),label="Discrete")
end
xlabel("Momentum")
ylabel("Energy")
title("Alpha = $(thisalpha)")
end=#

if false
lx = 6
n = 3
#for (idx,n) in enumerate([2,3,4,5])
#intstrens = [100.0]
#change = 0.001
#real_alphas = [range(0.1,0.21,length=10); range(0.22,0.28,length=10); range(0.29,0.35,length=5)]
#howmany = length(real_alphas)
#alphas = [real_alphas; real_alphas .+ change]
#all_bds = zeros(Float64,length(alphas))
#thetas = range(0.01,0.5,length=50)
#all_nrgs = zeros(Float64,length(thetas))
#anises = [1.0,10.0,100.0]
#for (idx,alph) in enumerate(alphas)
#for (idx,lx) in enumerate(4:1:30)
#for (idx,theta) in enumerate(thetas)
#for (idx,anis) in enumerate(anises)
#for (idx,intstren) in enumerate(intstrens)
#for lrd in [0,1]
    #for change in [0,0.0001]true
    params_dict = Dict([("Lx",lx),("N",n),("if_periodic_x",false),("if_periodic_y",false),("hopping_anisotropy",1.0),("interaction_strength",10.0),("lr","all"),("filling",0.5),("nev",10),("if_save_data",false)])
    #params_dict = make_args_dict(ARGS)

    # set number of open cores
    open_cores = get(params_dict, "open_cores", 5)
	if typeof(open_cores) != String
		BLAS.set_num_threads(open_cores)
		display(BLAS.get_config())
	end

    # set lattice parameters
    Lx = get(params_dict, "Lx", 4)
    Ly = get(params_dict, "Ly", Lx)
    N = get(params_dict, "N", 2)
    if_periodic_x = get(params_dict, "if_periodic_x", false)
    if_periodic_y = get(params_dict, "if_periodic_y", false)
    twist_angle = get(params_dict, "twist_angle", 0.0)


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
        basis_dataloc = dataloc
    end
    opl = get(params_dict, "output_level", 1)
    running_args = (nev=nev,
                    if_exact=false,
                    if_densmat=false,
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

    # build long range interaction parameters
    stren = get(params_dict,"interaction_strength", 0.0)
    lr_dist = get(params_dict,"lr", "all")
    lr_dist == "all" ? lr_dist = Ly : nothing
    us = [i < lr_dist+2 ? stren : 0.0 for i in 1:Ly]
    if lr_dist == 0.0
        us[1] = stren == 0.0 ? 1.0 : stren
    end
    #us[1] = 0.0

    # get hopping anisotropy values
    hopping_anisotropy = get(params_dict,"hopping_anisotropy",1.0)
    if hopping_anisotropy < 1.0
		ty = 1.0 / hopping_anisotropy
		tx = 1.0
	else
		tx = 1.0 * hopping_anisotropy
		ty = 1.0
	end
    #tx = 0.0
    #ty = 0.0
    
    # build magnetic field parameters and check fluxes for periodicity
    alpha = get(params_dict,"alpha",nothing)
    x_shift,y_shift = !if_periodic_x, !if_periodic_y
    if isnothing(alpha)
        filling = get(params_dict,"filling",0.5)
        alpha = N / (filling * (Lx - x_shift) * (Ly - y_shift))
    end
    check_fluxes(alpha,Lx,Ly,if_periodic_x,if_periodic_y)

    # build hamiltonian parameters dictionary
    hamilt_params = Dict("alpha"=>alpha,
                        "tx"=>tx,
                        "ty"=>ty,
                        "hopping_anisotropy"=>hopping_anisotropy,
                        "U"=>us,
                        "interaction_cutoff"=>1e-5)

    # build filename dictionary
    filename_dict = make_filename_dict(lattice_params,hamilt_params)
    if_exists,found_data = check_data_exists(filename_dict,"ed"; location=dataloc,output_level=false)

    # check if data exists and rerun if need more eigenstates
    if if_exists
        println("Found existing data: ",found_data[2]["filename"])
        if nev > found_data[2]["nev"]
            opl > 0 ? println("Asking for more eigenstates than in file, rerunning") : nothing
            start_time = time()
            full_basis = n_particle_basis(N,Lx,Ly; output_level=opl,dataloc=basis_dataloc)
            opl > 0 ? println("Made basis in ",time()-start_time) : nothing
            lattice_params["full_basis"] = full_basis 
            states,nrgs,rhos = rerun_eigenstates(nev,lattice_params,hamilt_params,found_data[2],found_data[1]; running_args...)
        elseif nev < found_data[2]["nev"]
            opl > 0 ? println("Asking for fewer eigenstates than in file, using existing data") : nothing
            states = found_data[1]["state"][1:nev]
            nrgs = found_data[1]["nrg"][1:nev]
            rhos = found_data[1]["densmat"][1:nev]
            full_basis = n_particle_basis(N,Lx,Ly; output_level=opl,dataloc=basis_dataloc)
            opl > 0 ? println("Made basis in ",time()-start_time) : nothing
            lattice_params["full_basis"] = full_basis 
        else
            states = found_data[1]["state"]
            nrgs = found_data[1]["nrg"]
            rhos = found_data[1]["densmat"]
            full_basis = n_particle_basis(N,Lx,Ly; output_level=opl,dataloc=basis_dataloc)
            opl > 0 ? println("Made basis in ",time()-start_time) : nothing
            lattice_params["full_basis"] = full_basis 
        end
    else

        # make basis only if data doesn't exist
        start_time = time()
        full_basis = n_particle_basis(N,Lx,Ly; output_level=opl,dataloc=basis_dataloc)
        opl > 0 ? println("Made basis in ",time()-start_time) : nothing
        lattice_params["full_basis"] = full_basis 

        # run exact diagonalization for find eigenstates
        if nev == 1
            states,nrgs,rhos = find_ground_state(lattice_params,hamilt_params; running_args...)
        else
            states,nrgs,rhos,hh = find_eigenstates(nev,lattice_params,hamilt_params; running_args...)
        end
    end

    #dcorrs = distance_correlation(states[1],lattice_params,"x")
    #display(dcorrs)
    #dcorrsy = distance_correlation(states[1],lattice_params,"y")
    #display(dcorrsy)
    #bulkdens = sum(occs[2,2:5]) + sum(occs[5,2:5]) + sum(occs[3:4,2]) + sum(occs[3:4,5])
    #bulkdens = iseven(lx) ? sum(occs[Int(lx/2):Int(lx/2)+1,Int(lx/2):Int(lx/2)+1]) : sum(occs[Int(floor(lx/2)):Int(floor(lx/2))+2,Int(floor(lx/2)):Int(floor(lx/2))+2])
    #append!(all_bulkdens,[bulkdens/n])

    #display(nrgs)

    #for i in 1:nev
    #    get_occupancy(rhos[i],lattice_params; if_plot=true,plot_title="Alpha=$alpha, E=$(round(nrgs[i]/Ly,digits=5))")
    #end

    #occs1 = get_occupancy(rhos[1],lattice_params; if_plot=true,plot_title="tx=$tx ty=$ty")
    #append!(all_bulkdens,[sum(occs1[3:4,3:4])])

    #for i in 1:1#nev
    #    occs = get_occupancy(rhos[i],lattice_params; if_plot=true,plot_title="$i E=$(round(nrgs[i],digits=5))")
    #end
    #coeff = (maximum(nrgs) .- nrgs[1]) / hh_gap_exact(anis,alpha)
    #append!(coeffs,[coeff])
    cols = ["b","r","g","m","c"]
    if nev > length(cols)
        cols = repeat(cols,ceil(Int,nev/length(cols)))
    end

    #xx = N / (alpha * (Lx - x_shift) * (Ly - y_shift))

    #gapvals[idx] = maximum(nrgs) - minimum(nrgs)
    
    #=nrggaps = nrgs .- nrgs[1]
    theory = 4*(sin(pi*alpha)^2)
    which_is_close = argmin(abs.(nrggaps .- theory))

    for i in 1:nev
        if i == which_is_close
            scatter(alpha,nrggaps[i],c="r")
        else
            scatter(alpha,nrggaps[i],c="k")
        end
    end=#

    #=if idx == 1
        for i in 1:nev
            scatter(intstren,nrgs[i] - nrgs[1],c=cols[i],label="E$i - E1")
        end
        #scatter(id2 == 1 ? anis : -anis,nrgs[2] - nrgs[1],c=cols[id2*2-1],label="E2 - E1")
        #scatter(id2 == 1 ? anis : -anis,nrgs[3] - nrgs[1],c=cols[2*id2],label="E3 - E1")
        #scatter(id2 == 1 ? anis : -anis,nrgs[1],c="b",label="E0")
        #scatter(id2 == 1 ? anis : -anis,nrgs[2],c="r",label="E1")
        #scatter(id2 == 1 ? anis : -anis,nrgs[3],c="g",label="E2")
        #scatter(intstren,nrgs[4],c="k",label="E3")
    else
        for i in 1:nev
            scatter(intstren,nrgs[i] - nrgs[1],c=cols[i])
        end
        #scatter(id2 == 1 ? anis : -anis,nrgs[2] - nrgs[1],c=c=cols[id2*2-1])
        #plot(anises[idx-1:idx],[hh_gap_exact(anises[idx-1],alpha),hh_gap_exact(anises[idx],alpha)],c="r")
        #scatter(id2 == 1 ? anis : -anis,nrgs[3] - nrgs[1],c=c=cols[id2*2])
        #scatter(id2 == 1 ? anis : -anis,nrgs[1],c="b")
        #scatter(id2 == 1 ? anis : -anis,nrgs[2],c="r")
        #scatter(id2 == 1 ? anis : -anis,nrgs[3],c="g")
        #scatter(intstren,nrgs[4],c="k")
    end
    legend()
    #xlabel("System Size")
    xlabel("Interaction Strength")
    #xlabel("Flux")
    #xlabel("Hopping Anisotropy tx/ty")
    ylabel("Energy - E1")=#

    #cdw = cdw_sf(rhos[1],states[1],lattice_params,(3.0,0.0); if_plot=true,plot_label="$anis")
    #occs = get_occupancy(states[1],lattice_params; if_plot=false,plot_title="ED")
    #all_bds[idx] = sum(occs[3:4,3:4]) / 4
    #scurr = synthetic_current(rhos[1],lattice_params; if_plot=true,plot_title="Int Stren=$intstren")
    #pcurr = physical_current(rhos[1],lattice_params; if_plot=true,plot_title="Int Stren=$intstren")


    #title("Energy Gap for 2 x $Lx decoupled wires")
    #scatter(intstren,nrgs[2] .- nrgs[1],c="b")
    #scatter(intstren,nrgs[3] .- nrgs[2],c="r")
    #if idx == 1
    #    plot(anises,ones(length(anises)) .* 4*(sin(pi*alpha)^2),c="g",label="TT-Gap")
    #end
    #scatter(anis,nrgs[2] - nrgs[1],c="b")
    #scatter(anis,nrgs[3] - nrgs[2],c="r")
    #scatter(anis,nrgs[3],c="g")
    #scatter(anis,hh_gap_exact(anis,alpha),c="g")
    #xlabel("Hopping Anisotropy")
    #xlabel("Interaction Strength")
    #ylabel("Energy")
    #ylim(-9.5,-8.5)
    #=corrs = physical_correlation(rhos,Lx,Ly; if_plot=true)
    currents = physical_current(rhos,lattice_params; if_plot=true)
    corrs_syn = synthetic_correlation(rhos,Lx,Ly; if_plot=true)
    currents_syn = synthetic_current(rhos,lattice_params; if_plot=true,plot_title="Int Stren=$stren")=#
#end

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
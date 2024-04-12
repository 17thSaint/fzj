using Pkg
Pkg.activate(".")
using LinearAlgebra,KrylovKit,Combinatorics,SparseArrays

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

# builds up the full basis of N hard-core particles in a Lx by Ly lattice (no symmetries yet)
function generate_basis(Lx::Int64,Ly::Int64,N::Int64; kwargs...)

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

generate_basis(L::Int64,N::Int64) = generate_basis(L,L,N)

function get_occupancy(x::Vector{ComplexF64},lattice_params::Dict{String,Any}; kwargs...)
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

function applyHam(which_basis::Int64,lattice_params::Dict,hamilt_params::Dict)
    
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

    return change_of_basis,sparse(pinv(Matrix(change_of_basis)))
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

    if_plot ? plot_physical_correlation(phys_corrs) : nothing

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

    if_plot ? plot_synthetic_correlation(syn_corrs) : nothing

    return syn_corrs
end

function plot_synthetic_correlation(syn_corrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(syn_corrs)[1]
        plot(1:size(syn_corrs)[2],syn_corrs[i,:],"-p",label="$i")
    end
    xlabel("Synthetic Site")
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

function buildHopping(lattice_params::Dict,site1::Int64,site2::Int64; kwargs...)
    output_level = get(kwargs,:output_level,1)
    full_basis = lattice_params["full_basis"]

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
    for i in 1:size(currents)[1]
        plot(1:size(currents)[2],currents[i,:],"-p",label="$i")
    end
    xlabel("Synthetic Site")
    ylabel("Current")
    title("Synthetic Current "*plot_title)
    legend()
    return nothing
end

function plot_physical_current(currents::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(currents)[2]
        plot(1:size(currents)[1],currents[:,i],"-p",label="$i")
    end
    xlabel("Physical Site")
    ylabel("Current")
    title("Physical Current "*plot_title)
    legend()
    return nothing
end

function check_fluxes(alpha::Float64,Lx::Int64,Ly::Int64,if_periodic_x::Bool,if_periodic_y::Bool)
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


if true

#density = 1/4
Lx,Ly = 8,8
N = 4#Int(floor(density*Lx*Ly))
println("Using ",N," particles with density ",round(N/(Lx*Ly),digits=3))
if_periodic_x,if_periodic_y = true,true
start_time = time()
full_basis,basis_dict = generate_basis(Lx,Ly,N; output_level=1)
println("Made basis in ",time()-start_time)
lattice_params::Dict{String,Any} = Dict("Lx"=>Lx,
                      "Ly"=>Ly,
                      "N"=>N,
                      "if_periodic_x"=>if_periodic_x,
                      "if_periodic_y"=>if_periodic_y,
                      "full_basis"=>full_basis,
                      "basis_dict"=>basis_dict)

stren = 0.0
lr_dist = "all"
lr_dist == "all" ? lr_dist = Ly : nothing
us = [i < lr_dist+1 ? stren : 0.0 for i in 1:Ly]    
filling = 0.5
x_shift,y_shift = !if_periodic_x, !if_periodic_y
alpha = N / (filling * (Lx - x_shift) * (Ly - y_shift))
check_fluxes(alpha,Lx,Ly,if_periodic_x,if_periodic_y)
hamilt_params = Dict("alpha"=>alpha,
                     "tx"=>1.0,
                     "ty"=>1.0,
                     "U"=>us,
                     "interaction_cutoff"=>1e-5)

#
H = buildHam(lattice_params,hamilt_params)
println("Sparsity = ",nnz(H)/size(H)[1]^2)

nev = 20
rez = eigsolve(H,nev)
println("Ground State: Elapsed time: ",time()-start_time)
gs = rez[2][findfirst(x->x==minimum(rez[1]),rez[1])]#


rho = density_matrix(gs,lattice_params)
occs = get_occupancy(rho,lattice_params; if_plot=true)
corrs = physical_correlation(rho,Lx,Ly; if_plot=true)
currents = physical_current(rho,lattice_params; if_plot=true)
corrs_syn = synthetic_correlation(rho,Lx,Ly; if_plot=true)
currents_syn = synthetic_current(rho,lattice_params; if_plot=true)


end

#=nrgs_krylov = filter(x->x<0.0,rez[1])
fig = figure()
scatter(1:length(nrgs_krylov),nrgs_krylov,c="b",label="Krylov")
xlabel("Eigenvalue Index")
ylabel("Energy")
legend()

fig2 = figure()
scatter(1:10,rr[end-9:end],c="b")
xlabel("Eigenvalue Index")
ylabel("Eigenvalue")
title("Density Matrix Eigenvalues")
yscale("log")=#







































"fin"
#####################################################
#=

This file is for building more efficient ED Hamiltonians by pulling from saved undressed hopping and interaction matrices

Depends on:
    other-funcs/data-storage-funcs.jl

=#
######################################################


function buildHopping(lattice_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if_save = get(kwargs,:if_save,true)
    if output_level > 1
        println("Building hopping matrix")
    end
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    if_periodic_x = lattice_params["if_periodic_x"]
    if_periodic_y = lattice_params["if_periodic_y"]
    

    # initialize the hopping Matrix
    hopping = spzeros(ComplexF64,size(lattice_params["full_basis"])[2],size(lattice_params["full_basis"])[2])

    for j in 1:Lx
        for s in 1:Ly
            starting_site::Tuple{Int64,Int64} = (j,s)
            for dir in [(1,0),(-1,0),(0,1),(0,-1)]
            
                # skip if at boundary and no periodic boundary
                # x-direction
                if !if_periodic_x && ((starting_site .+ dir)[1] < 1 || (starting_site .+ dir)[1] > Lx)
                    continue
                end

                # y-direction
                if !if_periodic_y && ((starting_site .+ dir)[2] < 1 || (starting_site .+ dir)[2] > Ly)
                    continue
                end

                # find the next site modulo the system size
                next_site = (mod1(starting_site[1]+dir[1],Lx),mod1(starting_site[2]+dir[2],Ly))

                # add local hopping sparse matrix to full one
                hopping .+= buildHopping(lattice_params,linear_index(starting_site,Lx,Ly),linear_index(next_site,Lx,Ly))

            end
        end
    end

    if_save && saveHopping(hopping,get_folder_location("cluster-data/exact-diag"),lattice_params; kwargs...)

    return hopping
end

function dressHopping(hamilt_params::Dict,lattice_params::Dict,hopping::SparseMatrixCSC{ComplexF64}; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if output_level > 1
        println("Dressing hopping matrix")
    end
    
    rows,cols,vals = findnz(hopping)

    twist_angle = lattice_params["twist_angle"]
    tx = hamilt_params["tx"]
    ty = hamilt_params["ty"]
    alpha = hamilt_params["alpha"]
    flux_direction = hamilt_params["flux_direction"]

    for idx in 1:length(rows)

        # find the starting and ending basis
        starting_basis::Vector{Int64} = lattice_params["full_basis"][:,rows[idx]]
        ending_basis::Vector{Int64} = lattice_params["full_basis"][:,cols[idx]]

        # find the hopping vector
        starting_linear_position = setdiff(starting_basis,ending_basis)[1]
        starting_site = coordinate(starting_linear_position,lattice_params["Lx"],lattice_params["Ly"])
        ending_linear_position = setdiff(ending_basis,starting_basis)[1]
        dir = coordinate(ending_linear_position,lattice_params["Lx"],lattice_params["Ly"]) .- coordinate(starting_linear_position,lattice_params["Lx"],lattice_params["Ly"])
        
        # correct hopping vector for periodic boundary conditions
        (dir[1] > 1) && (dir = (-1,0))
        (dir[1] < -1) && (dir = (1,0))
        (dir[2] > 1) && (dir = (0,-1))
        (dir[2] < -1) && (dir = (0,1))

        # hopping amplitude from tx/ty
        coeff = abs(dir[1]) == 0 ? -ty : -tx

        # flux attachment
        coeff *= exp(im*dot(alpha,dir)*dot(starting_site,abs.(reverse(dir)))*2*pi)

        # boundary condition twisting
        coeff *= exp(im*2*pi*dot(twist_angle ./ (lattice_params["Lx"],lattice_params["Ly"]),dir))

        # add the hopping term to the Hamiltonian
        hopping[rows[idx],cols[idx]] *= coeff

    end

    return hopping
end

function saveHopping(hopping::SparseMatrixCSC{ComplexF64},dataloc::String,lattice_params::Dict; kwargs...)
    data_dict = Dict([("hopping_matrix",hopping)])
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N)])
    filename = "hopping-N-"*string(N)*"-Lx-"*string(Lx)*"-Ly-"*string(Ly)*".jld2"
    full_loc = join([dataloc,filename],"/")
    write_data_jld2(full_loc,data_dict,metadata_dict; kwargs...)
end

function findHopping(lattice_params::Dict; kwargs...)
    dataloc = get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag"))
    output_level = get(kwargs,:output_level,1)

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N)])
    return check_data_exists(metadata_dict,"hopping"; location=dataloc,output_level=output_level)
end

function getHopping(lattice_params::Dict,hamilt_params::Dict; kwargs...)
    if_found,hopping_data = findHopping(lattice_params; kwargs...)

    !if_found && (hopping = buildHopping(lattice_params; kwargs...))
    if_found && (hopping = hopping_data[1]["hopping_matrix"])
    
    hopping = dressHopping(hamilt_params,lattice_params,hopping; kwargs...)
    
    return hopping
end

function buildInteraction(lattice_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if_save::Bool = get(kwargs,:if_save,true)
    if output_level > 1
        println("Building interaction matrix")
    end

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

    # initialize the interaction Matrix
    interaction = spzeros(ComplexF64,size(lattice_params["full_basis"])[2],size(lattice_params["full_basis"])[2])

    Threads.@threads for idx in 1:size(interaction,1)
        basis = lattice_params["full_basis"][:,idx]
        for s1 in basis
            for s2 in basis
                # skip if the same particle
                s1 == s2 && continue

                # find the distance between the two particles
                dist = coordinate(s2,Lx,Ly) .- coordinate(s1,Lx,Ly)

                # skip if the particles aren't on the same physical index
                dist[1] != 0 && continue

                # add the interaction term to the Hamiltonian
                interaction[idx,idx] += 0.5

            end
        end
    end
    
    if_save && saveInteraction(interaction,get_folder_location("cluster-data/exact-diag"),lattice_params; kwargs...)

    return interaction
end

function dressInteraction(hamilt_params::Dict,lattice_params::Dict,interaction::SparseMatrixCSC{ComplexF64}; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if output_level > 1
        println("Dressing interaction matrix")
    end

    U::Vector{Float64} = hamilt_params["U"]
    Lx::Int64 = lattice_params["Lx"]
    Ly::Int64 = lattice_params["Ly"]

    if length(unique(U[2:end])) == 1
        interaction = dressInteraction(U[2],lattice_params,interaction; kwargs...)
        
        # onsite pinning potential
        (haskey(hamilt_params,"if_pinning") && hamilt_params["if_pinning"]) && (addPinning(interaction,lattice_params,hamilt_params))

        return interaction
    end

    rows,cols,vals = findnz(interaction)

    for idx in 1:length(rows)
        basis = lattice_params["full_basis"][:,rows[idx]]
        for s1 in basis
            for s2 in basis
                # skip if the same particle
                s1 == s2 && continue

                # find the distance between the two particles
                dist = coordinate(s2,Lx,Ly) .- coordinate(s1,Lx,Ly)

                # skip if the particles aren't on the same physical index
                dist[1] != 0 && continue

                # add the interaction term to the Hamiltonian
                interaction[idx,idx] *= U[abs(dist[2])]

            end
        end
    end

    # onsite pinning potential
    (haskey(hamilt_params,"if_pinning") && hamilt_params["if_pinning"]) && (addPinning(interaction,lattice_params,hamilt_params))

    return interaction
end

function dressInteraction(flat_intstren::Float64,lattice_params::Dict,interaction::SparseMatrixCSC{ComplexF64}; kwargs...)
    interaction .*= flat_intstren
    return interaction
end

function saveInteraction(interaction::SparseMatrixCSC{ComplexF64},dataloc::String,lattice_params::Dict; kwargs...)
    data_dict = Dict([("interaction_matrix",interaction)])
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N)])
    filename = "interaction-N-"*string(N)*"-Lx-"*string(Lx)*"-Ly-"*string(Ly)*".jld2"
    full_loc = join([dataloc,filename],"/")
    write_data_jld2(full_loc,data_dict,metadata_dict; kwargs...)
end

function findInteraction(lattice_params::Dict; kwargs...)
    dataloc = get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag"))
    output_level = get(kwargs,:output_level,1)

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N)])
    return check_data_exists(metadata_dict,"interaction"; location=dataloc,output_level=output_level)
end

function getInteraction(lattice_params::Dict,hamilt_params::Dict; kwargs...)
    if_found,interaction_data = findInteraction(lattice_params; kwargs...)

    !if_found && (interaction = buildInteraction(lattice_params; kwargs...))
    if_found && (interaction = interaction_data[1]["interaction_matrix"])
    
    interaction = dressInteraction(hamilt_params,lattice_params,interaction; kwargs...)
    
    return interaction
end

function getHamiltonian(lattice_params::Dict,hamilt_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if output_level > 1
        println("Building Hamiltonian by dressing")
    end

    if lattice_params["twist_angle"] != [0.0,0.0]
        error("Twistings not yet implemented")
    end

    hopping::SparseMatrixCSC{ComplexF64} = getHopping(lattice_params,hamilt_params; kwargs...)
    interaction::SparseMatrixCSC{ComplexF64} = getInteraction(lattice_params,hamilt_params; kwargs...)

    return hopping + interaction
end




































"fin"
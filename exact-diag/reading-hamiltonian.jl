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

    if_save && saveHopping(hopping,get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag")),lattice_params; kwargs...)

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

        output_level > 0 && idx % 1000 == 0 && println("$((idx/length(rows))*100)% complete")

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

        # flux attachment (0-indexed coordinates to match applyHam) and boundary condition twisting
        phase(d) = exp(im*dot(alpha,d)*dot(starting_site .- 1,abs.(reverse(d)))*2*pi) * exp(im*2*pi*dot(twist_angle ./ (lattice_params["Lx"],lattice_params["Ly"]),d))

        # a periodic direction of length 2 merges the +dir and -dir bonds into one entry, so average their phases
        if_double_bond = abs(dir[1]) == 1 ? (lattice_params["Lx"] == 2 && lattice_params["if_periodic_x"]) : (lattice_params["Ly"] == 2 && lattice_params["if_periodic_y"])
        coeff *= if_double_bond ? (phase(dir) + phase(dir .* -1))/2 : phase(dir)

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

    # the undressed hopping structure depends on the boundary conditions so they belong in the cache key
    if_periodic_x = lattice_params["if_periodic_x"]
    if_periodic_y = lattice_params["if_periodic_y"]
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",if_periodic_x),("if_periodic_y",if_periodic_y)])
    filename = "hopping-N-"*string(N)*"-Lx-"*string(Lx)*"-Ly-"*string(Ly)*"-if_periodic_x-"*string(if_periodic_x)*"-if_periodic_y-"*string(if_periodic_y)*".jld2"
    full_loc = join([dataloc,filename],"/")
    write_data_jld2(full_loc,data_dict,metadata_dict; kwargs...)
end

function findHopping(lattice_params::Dict; kwargs...)
    dataloc = get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag"))
    output_level = get(kwargs,:output_level,1)

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",lattice_params["if_periodic_x"]),("if_periodic_y",lattice_params["if_periodic_y"])])
    return check_data_exists(metadata_dict,"hopping"; location=dataloc,output_level=output_level,file_type="jld2")
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

    # count interacting pairs per basis state (plain vector so the threaded loop is race-free)
    dimHilb = size(lattice_params["full_basis"])[2]
    pair_counts = zeros(Float64,dimHilb)

    Threads.@threads for idx in 1:dimHilb
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
                pair_counts[idx] += 0.5

            end
        end
    end

    # build the diagonal interaction Matrix, keeping only states with interacting pairs
    interaction = dropzeros!(spdiagm(0 => ComplexF64.(pair_counts)))

    if_save && saveInteraction(interaction,get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag")),lattice_params; kwargs...)

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

    interaction_cutoff::Float64 = get(hamilt_params,"interaction_cutoff",1e-5)
    lr_dist = sum(abs.(U) .> interaction_cutoff) - 1

    rows,cols,vals = findnz(interaction)

    for idx in 1:length(rows)
        basis = lattice_params["full_basis"][:,rows[idx]]

        # the undressed entry only counts pairs, so rebuild the value as a sum over pairs like applyHam
        dressed_value = 0.0
        for i in 1:length(basis)
            for j in i+1:length(basis)

                # find the distance between the two particles
                dist = coordinate(basis[j],Lx,Ly) .- coordinate(basis[i],Lx,Ly)

                # skip if the particles aren't on the same physical index
                dist[1] != 0 && continue

                # add the interaction term for this pair
                abs(dist[2]) <= lr_dist && abs(U[abs(dist[2])+1]) > interaction_cutoff && (dressed_value += U[abs(dist[2])+1])

            end
        end
        interaction[rows[idx],cols[idx]] = dressed_value
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
    return check_data_exists(metadata_dict,"interaction"; location=dataloc,output_level=output_level,file_type="jld2")
end

function getInteraction(lattice_params::Dict,hamilt_params::Dict; kwargs...)
    if_found,interaction_data = findInteraction(lattice_params; kwargs...)

    !if_found && (interaction = buildInteraction(lattice_params; kwargs...))
    if_found && (interaction = interaction_data[1]["interaction_matrix"])
    
    interaction = dressInteraction(hamilt_params,lattice_params,interaction; kwargs...)
    
    return interaction
end

function buildPeriodicPotential(lattice_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if_save::Bool = get(kwargs,:if_save,true)
    if output_level > 1
        println("Building periodic potential shape")
    end

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    full_basis = lattice_params["full_basis"]
    dimHilb = size(full_basis)[2]

    # staggered shape along the physical direction, one term per particle (strength-independent)
    shape_values = zeros(Float64,dimHilb)
    Threads.@threads for idx in 1:dimHilb
        for s in full_basis[:,idx]
            shape_values[idx] += (-1)^(coordinate(s,Lx,Ly)[1])
        end
    end
    potential = dropzeros!(spdiagm(0 => ComplexF64.(shape_values)))

    if_save && savePeriodicPotential(potential,get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag")),lattice_params; kwargs...)

    return potential
end

function dressPeriodicPotential(hamilt_params::Dict,lattice_params::Dict,potential::SparseMatrixCSC{ComplexF64}; kwargs...)
    potential .*= hamilt_params["periodic_potential_strength"]
    return potential
end

function savePeriodicPotential(potential::SparseMatrixCSC{ComplexF64},dataloc::String,lattice_params::Dict; kwargs...)
    data_dict = Dict([("periodicpotential_matrix",potential)])
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N)])
    filename = "periodicpotential-N-"*string(N)*"-Lx-"*string(Lx)*"-Ly-"*string(Ly)*".jld2"
    full_loc = join([dataloc,filename],"/")
    write_data_jld2(full_loc,data_dict,metadata_dict; kwargs...)
end

function findPeriodicPotential(lattice_params::Dict; kwargs...)
    dataloc = get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag"))
    output_level = get(kwargs,:output_level,1)

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N)])
    return check_data_exists(metadata_dict,"periodicpotential"; location=dataloc,output_level=output_level,file_type="jld2")
end

function getPeriodicPotential(lattice_params::Dict,hamilt_params::Dict; kwargs...)
    if_found,potential_data = findPeriodicPotential(lattice_params; kwargs...)

    !if_found && (potential = buildPeriodicPotential(lattice_params; kwargs...))
    if_found && (potential = potential_data[1]["periodicpotential_matrix"])

    potential = dressPeriodicPotential(hamilt_params,lattice_params,potential; kwargs...)

    return potential
end

function buildDisorder(lattice_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if_save::Bool = get(kwargs,:if_save,true)
    if output_level > 1
        println("Building disorder shape")
    end

    dimHilb = size(lattice_params["full_basis"])[2]

    # one random onsite energy per basis state in [-1,1]; saving the shape makes the
    # disorder realization quenched (reused for every dressing), unlike applyHam
    # which draws a fresh realization on every build
    disorder = spdiagm(0 => ComplexF64.(rand(dimHilb) .* 2 .- 1))

    if_save && saveDisorder(disorder,get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag")),lattice_params; kwargs...)

    return disorder
end

function dressDisorder(hamilt_params::Dict,lattice_params::Dict,disorder::SparseMatrixCSC{ComplexF64}; kwargs...)
    disorder .*= hamilt_params["disorder_strength"]
    return disorder
end

function saveDisorder(disorder::SparseMatrixCSC{ComplexF64},dataloc::String,lattice_params::Dict; kwargs...)
    data_dict = Dict([("disorder_matrix",disorder)])
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N)])
    filename = "disorder-N-"*string(N)*"-Lx-"*string(Lx)*"-Ly-"*string(Ly)*".jld2"
    full_loc = join([dataloc,filename],"/")
    write_data_jld2(full_loc,data_dict,metadata_dict; kwargs...)
end

function findDisorder(lattice_params::Dict; kwargs...)
    dataloc = get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag"))
    output_level = get(kwargs,:output_level,1)

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    N = lattice_params["N"]
    metadata_dict = Dict([("Lx",Lx),("Ly",Ly),("N",N)])
    return check_data_exists(metadata_dict,"disorder"; location=dataloc,output_level=output_level,file_type="jld2")
end

function getDisorder(lattice_params::Dict,hamilt_params::Dict; kwargs...)
    if_found,disorder_data = findDisorder(lattice_params; kwargs...)

    !if_found && (disorder = buildDisorder(lattice_params; kwargs...))
    if_found && (disorder = disorder_data[1]["disorder_matrix"])

    disorder = dressDisorder(hamilt_params,lattice_params,disorder; kwargs...)

    return disorder
end

function getHamiltonian(lattice_params::Dict,hamilt_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
    if output_level > 1
        println("Building Hamiltonian by dressing")
    end

    if lattice_params["twist_angle"] != [0.0,0.0]
        error("Twistings not yet implemented")
    end

    if get(hamilt_params,"which_dir","virt") != "virt"
        error("Interactions only implemented for which_dir = virt")
    end

    hopping::SparseMatrixCSC{ComplexF64} = getHopping(lattice_params,hamilt_params; kwargs...)
    interaction::SparseMatrixCSC{ComplexF64} = getInteraction(lattice_params,hamilt_params; kwargs...)

    ham = hopping + interaction

    # diagonal potential shapes read from file and dressed by their strengths
    (get(hamilt_params,"disorder_strength",0.0) != 0.0) && (ham += getDisorder(lattice_params,hamilt_params; kwargs...))
    (get(hamilt_params,"periodic_potential_strength",0.0) != 0.0) && (ham += getPeriodicPotential(lattice_params,hamilt_params; kwargs...))

    return ham
end




































"fin"
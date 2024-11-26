#####################################################
#=

This file is for testing a way of building more efficient ED Hamiltonians

Depends on:
    exact-diag/execute-ed.jl

=#
######################################################

include("execute-ed.jl")
using Test

function buildHopping(lattice_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
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

function buildInteraction(lattice_params::Dict; kwargs...)
    output_level = get(kwargs,:output_level,1)
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
        return dressInteraction(U[2],lattice_params,interaction; kwargs...)
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

    return interaction
end

function dressInteraction(flat_intstren::Float64,lattice_params::Dict,interaction::SparseMatrixCSC{ComplexF64}; kwargs...)
    return flat_intstren .* interaction
end


do_all::Bool = true

if do_all || false
@testset "Only Hopping with Anis/Fluxes" begin;

    configs = [(4,4,2),(4,3,2),(4,4,3),(3,2,2),(4,3,3)]
    tw1 = 0.0
    tw2 = 0.0
    #tws = range(0.0,1.0,length=5)
    hanises = [1.0,0.5,2.0]
    for if_periodic_x in [true,false]
        for if_periodic_y in [true,false]
            for config in configs
                lx,ly,n = config
                for hanis in hanises
                    #println("Working on Lx = ",lx," Ly = ",ly," N = ",n," tw1 = ",tw1," tw2 = ",tw2," hanis = ",hanis)
                    pdict = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("tw1",tw1),("tw2",tw2),("if_check_fluxes",false),("if_pinning",false),("if_periodic_x",if_periodic_x),("if_periodic_y",if_periodic_y),("hopping_anisotropy",hanis),("interaction_strength",0.0),("lr",0),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
                    lattice_params,hamilt_params,running_args = get_normal_model_params_ed(pdict)    
                    basis_dataloc = running_args.basis_dataloc
                    full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
                    lattice_params["full_basis"] = full_basis
                    ham_correct = buildHam(lattice_params,hamilt_params; output_level=0)

                    ham_new = buildHopping(lattice_params; output_level=0)
                    ham_new = dressHopping(hamilt_params,lattice_params,ham_new; output_level=0)

                    @test isapprox(ham_correct,ham_new,atol=1e-10)
                end
            end
        end
    end

end;
end

if do_all || false
@testset "Interaction Part" begin;
    configs = [(4,4,2),(4,3,2),(4,4,3),(4,3,3)]
    tw1 = 0.0
    tw2 = 0.0
    #tws = range(0.0,1.0,length=5)
    hanises = [1.0,0.5,2.0]
    intstrens = [1.0,2.0,3.0]

    for if_periodic_x in [true,false]
        for if_periodic_y in [true,false]
            for config in configs
            lx,ly,n = config
                for hanis in hanises
                    for intstren in intstrens
                        #println("Working on Lx = ",lx," Ly = ",ly," N = ",n," tw1 = ",tw1," tw2 = ",tw2," hanis = ",hanis)
                        pdict = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("tw1",tw1),("tw2",tw2),("if_check_fluxes",false),("if_pinning",false),("if_periodic_x",if_periodic_x),("if_periodic_y",if_periodic_y),("hopping_anisotropy",hanis),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
                        lattice_params,hamilt_params,running_args = get_normal_model_params_ed(pdict)    
                        basis_dataloc = running_args.basis_dataloc
                        full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
                        lattice_params["full_basis"] = full_basis
                        ham_correct = buildHam(lattice_params,hamilt_params; output_level=0)

                        ham_new = buildHopping(lattice_params; output_level=0)
                        ham_new = dressHopping(hamilt_params,lattice_params,ham_new; output_level=0)
                        intham = buildInteraction(lattice_params; output_level=0)
                        intham = dressInteraction(hamilt_params,lattice_params,intham; output_level=0)
                        ham_new .+= intham

                        @test isapprox(ham_correct,ham_new,atol=1e-10)
                    end
                end
            end
        end
    end
end;
end

#=if do_all || false
@testset "Twistings" begin;
    
end;
end=#


















































"fin"
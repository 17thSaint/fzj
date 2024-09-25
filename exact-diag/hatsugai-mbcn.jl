#####################################################
#=

This file is for calculating Hatsugai-style many body Chern number for ED

Depends on:
    other-funcs/data-storage-funcs.jl
    
=#
######################################################


function get_reference_multiplets(Lx::Int64,Ly::Int64,particles::Int64; kwargs...)
    dataloc = get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag/torus"))

    reference_multiplets::Vector{Vector{ComplexF64}} = [zeros(ComplexF64,2) for i in 1:4]

    params_dict1 = Dict([("Lx",Lx),("Ly",Ly),("N",particles),("if_periodic_x",true),("if_periodic_y",true),("twist_angle1",0.33),("twist_angle2",0.67),("interaction_strength",0.0)])
    params_dict2 = Dict([("Lx",Lx),("Ly",Ly),("N",particles),("if_periodic_x",true),("if_periodic_y",true),("twist_angle1",0.67),("twist_angle2",0.33),("interaction_strength",0.0)])

    if_exists1,found_data1 = check_data_exists(params_dict1,"ed"; location=dataloc,output_level=false)
    if_exists2,found_data2 = check_data_exists(params_dict2,"ed"; location=dataloc,output_level=false)

    if if_exists1 && if_exists2
        reference_multiplets[1:2] = found_data1[1]["state"][1:2]
        reference_multiplets[3:4] = found_data2[1]["state"][1:2]
    else
        error("Data does not exist for reference multiplets")
    end

    return reference_multiplets,found_data1[2]["filename"],found_data2[2]["filename"]
end

function get_gamma(gs1::Vector{ComplexF64},gs2::Vector{ComplexF64},reference_multiplet::Vector{Vector{ComplexF64}})
    gamma_mat = zeros(ComplexF64,2,2)
    for i in 1:2
        for j in 1:2
            gamma_mat[i,j] = dot(conj(reference_multiplet[i]),gs1) * dot(conj(gs1),reference_multiplet[j]) + dot(conj(reference_multiplet[i]),gs2) * dot(conj(gs2),reference_multiplet[j])
        end
    end
    return det(gamma_mat)
end

function get_omega(gs1::Vector{ComplexF64},gs2::Vector{ComplexF64},reference_multiplets::Vector{Vector{ComplexF64}})
    omega_mat = zeros(ComplexF64,2,2)
    for i in 1:2
        for j in 1:2
            omega_mat[i,j] = dot(conj(reference_multiplets[i+2]),gs1) * dot(conj(gs1),reference_multiplets[j]) + dot(conj(reference_multiplets[i+2]),gs2) * dot(conj(gs2),reference_multiplets[j])
        end
    end

    return det(omega_mat)
end














































"fin"
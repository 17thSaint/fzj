#####################################################
#=

This file is for calculating Hatsugai-style many body Chern number for ED

Depends on:
    other-funcs/data-storage-funcs.jl

=#
######################################################


function get_reference_multiplets(Lx::Int64,Ly::Int64,particles::Int64; kwargs...)
    dataloc::String = get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag/torus"))

    reference_multiplets::Vector{Vector{ComplexF64}} = [zeros(ComplexF64,2) for i in 1:4]

    params_dict1::Dict{String,Any} = Dict([("Lx",Lx),("Ly",Ly),("N",particles),("if_periodic_x",true),("if_periodic_y",true),("twist_angle1",0.33),("twist_angle2",0.67),("interaction_strength",0.0)])
    params_dict2::Dict{String,Any} = Dict([("Lx",Lx),("Ly",Ly),("N",particles),("if_periodic_x",true),("if_periodic_y",true),("twist_angle1",0.67),("twist_angle2",0.33),("interaction_strength",0.0)])

    if_exists1::Bool,found_data1::Union{Vector{Dict},Nothing} = check_data_exists(params_dict1,"ed"; location=dataloc,output_level=false)
    if_exists2::Bool,found_data2::Union{Vector{Dict},Nothing} = check_data_exists(params_dict2,"ed"; location=dataloc,output_level=false)

    if if_exists1 && if_exists2
        reference_multiplets[1:2] = found_data1[1]["state"][1:2]
        reference_multiplets[3:4] = found_data2[1]["state"][1:2]
    else
        error("Data does not exist for reference multiplets")
    end

    return reference_multiplets,found_data1[2]["filename"],found_data2[2]["filename"]
end

function save_hatsugai_data(data_dict::Dict{String,ComplexF64},filepath::String,ref_multis_filenames::Vector{String})
    # put in some checks on the filepath
    for (idx,rmf) in enumerate(ref_multis_filenames)
        data_dict["rm$(idx)_name"] = rmf
    end
    modify_data_jld2(data_dict,filepath,"metadata"; output_level=1)
end

function save_hatsugai_data(data_dict::Dict{String,ComplexF64},filepath::String,ref_multis_filename::String,which_gamma::Int)
    # put in some checks on the filepath
    data_dict["rm$(which_gamma)_name"] = ref_multis_filename
    modify_data_jld2(data_dict,filepath,"metadata"; output_level=1)
end

function get_gamma(gs1::Vector{ComplexF64},gs2::Vector{ComplexF64},reference_multiplet::Vector{Vector{ComplexF64}}; kwargs...)
    if_save = get(kwargs,:if_save,false)

    gamma_mat::Matrix{ComplexF64} = zeros(ComplexF64,2,2)
    for i in 1:2
        for j in 1:2
            gamma_mat[i,j] = dot(conj(reference_multiplet[i]),gs1) * dot(conj(gs1),reference_multiplet[j]) + dot(conj(reference_multiplet[i]),gs2) * dot(conj(gs2),reference_multiplet[j])
        end
    end

    result = det(gamma_mat)

    if if_save
        which_gamma::Int = kwargs[:which_gamma]
        data_dict::Dict{String,ComplexF64} = Dict([("gamma$(which_gamma)",result)])
        filepath::String = kwargs[:filepath]
        save_hatsugai_data(data_dict,filepath,kwargs[:ref_multis_filename],which_gamma)
    end

    return result
end

function get_omega(gs1::Vector{ComplexF64},gs2::Vector{ComplexF64},reference_multiplets::Vector{Vector{ComplexF64}}; kwargs...)
    if_save = get(kwargs,:if_save,false)

    omega_mat::Matrix{ComplexF64} = zeros(ComplexF64,2,2)
    for i in 1:2
        for j in 1:2
            omega_mat[i,j] = dot(conj(reference_multiplets[i+2]),gs1) * dot(conj(gs1),reference_multiplets[j]) + dot(conj(reference_multiplets[i+2]),gs2) * dot(conj(gs2),reference_multiplets[j])
        end
    end

    results = det(omega_mat)

    if if_save
        data_dict::Dict{String,ComplexF64} = Dict([("omega",result)])
        filepath::String = kwargs[:filepath]
        ref_multis_filenames::Vector{String} = kwargs[:ref_multis_filenames]
        save_hatsugai_data(data_dict,filepath,ref_multis_filenames)
    end

    return result
end

function get_hatsugaifull(gs1::Vector{ComplexF64},gs2::Vector{ComplexF64},reference_multiplets::Vector{Vector{ComplexF64}}; kwargs...)
    if_save = get(kwargs,:if_save,false)

    gamma1::ComplexF64 = get_gamma(gs1,gs2,reference_multiplets[1:2]; if_save=false)
    gamma1::ComplexF64 = get_gamma(gs1,gs2,reference_multiplets[3:4]; if_save=false)
    omega::ComplexF64 = get_omega(gs1,gs2,reference_multiplets; if_save=false)

    if if_save
        data_dict::Dict{String,ComplexF64} = Dict([("gamma1",gamma1),("gamma2",gamma2),("omega",omega)])
        filepath::String = kwargs[:filepath]
        ref_multis_filenames::Vector{String} = kwargs[:ref_multis_filenames]
        save_hatsugai_data(data_dict,filepath,ref_multis_filenames)
    end

    return gamma1,gamma2,omega
end














































"fin"
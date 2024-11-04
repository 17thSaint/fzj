#####################################################
#=

This file is for calculating Hatsugai-style many body Chern number for effective 1D MPS

Depends on:
    other-funcs/data-storage-funcs.jl

=#
######################################################


# needs work with perturbative strength
function make_new_reference_multiplets(Lx::Int64,Ly::Int64,particles::Int64; kwargs...)
    tw1::Float64 = get(kwargs,:tw1,0.33)
    tw2::Float64 = get(kwargs,:tw2,0.67)
    hanis::Float64 = get(kwargs,:hopping_anisotropy,1.0)

    println("Making new reference multiplets for $(Lx)x$(Ly) N=$(particles) with Twists $(tw1) and $(tw2)")
    params_dict1 = Dict([("es_count",1),("if_remapping",false),("Lphys",Lx),("Lsynth",Ly),("particles",particles),("tw1",tw1),("tw2",tw2),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",hanis),("filling",0.5),("nev",4),("if_find_data",true),("if_save_data",true),("output_level",0)])
    params_dict2 = Dict([("es_count",1),("if_remapping",false),("Lphys",Lx),("Lsynth",Ly),("particles",particles),("tw1",tw2),("tw2",tw1),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",hanis),("filling",0.5),("nev",4),("if_find_data",true),("if_save_data",true),("output_level",0)])

    # make this for mps
    rez1 = run_normal_1deffmps(params_dict1)
    rez2 = run_normal_1deffmps(params_dict2)

    return [rez1[1][1],rez1[1][2],rez2[1][1],rez2[1][2]],rez1[4][:name],rez2[4][:name]
end

function get_reference_multiplets(Lx::Int64,Ly::Int64,particles::Int64; kwargs...)
    if_make_new::Bool = get(kwargs,:if_make_new,false)
    dataloc::String = get(kwargs,:dataloc,get_folder_location("cluster-data/synth-dims/twists"))
    tw1::Float64 = get(kwargs,:tw1,0.33)
    tw2::Float64 = get(kwargs,:tw2,0.67)
    hanis::Float64 = get(kwargs,:hopping_anisotropy,1.0)

    reference_multiplets = Vector{MPS}(undef,4)

    params_dict1::Dict{String,Any} = Dict([("Lphys",Lx),("Lsynth",Ly),("nbosons",particles),("if_periodic_phys",true),("if_periodic_synth",true),("tw1",tw1),("tw2",tw2),("hopping_anisotropy",hanis)])
    params_dict2::Dict{String,Any} = Dict([("Lphys",Lx),("Lsynth",Ly),("nbosons",particles),("if_periodic_phys",true),("if_periodic_synth",true),("tw1",tw2),("tw2",tw1),("hopping_anisotropy",hanis)])

    if_exists1::Bool,metadata1::Union{Vector{Dict},Nothing} = check_data_exists(params_dict1,"mps"; location=dataloc,output_level=false)
    if_exists2::Bool,metadata2::Union{Vector{Dict},Nothing} = check_data_exists(params_dict2,"mps"; location=dataloc,output_level=false)

    if if_exists1 && if_exists2
        found_data1::Dict{String,Any} = read_data_jld2("wavefunc"*metadata1[2]["name"],dataloc; output_level=0)
        reference_multiplets[1] = found_data1["mps"]
        reference_multiplets[2] = found_data1["mps_1"]

        found_data2::Dict{String,Any} = read_data_jld2("wavefunc"*metadata2[2]["name"],dataloc; output_level=0)
        reference_multiplets[3] = found_data2["mps"]
        reference_multiplets[4] = found_data2["mps_1"]
    elseif if_make_new
        return make_new_reference_multiplets(Lx,Ly,particles; tw1=tw1,tw2=tw2,hopping_anisotropy=hanis)
    else
        error("Data does not exist for reference multiplets")
    end

    return reference_multiplets,metadata1[2]["name"],metadata2[2]["name"]
end

function save_hatsugai_data(data_dict::Dict{String,Any},filepath::String,ref_multis_filenames::Vector{String})
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

# needs work
function get_gamma(gs1::MPS,gs2::MPS,reference_multiplet::Vector{MPS}; kwargs...)
    if_save = get(kwargs,:if_save,false)

    gamma_mat::Matrix{ComplexF64} = zeros(ComplexF64,2,2)
    for i in 1:2
        for j in 1:2
            #gamma_mat[i,j] = (transpose(conj(reference_multiplet[i])) * gs1) * (transpose(conj(gs1)) * reference_multiplet[j]) + (transpose(conj(reference_multiplet[i])) * gs2) * (transpose(conj(gs2)) * reference_multiplet[j])
            gamma_mat[i,j] = inner(reference_multiplet[i],gs1) * inner(gs1,reference_multiplet[j]) + inner(reference_multiplet[i],gs2) * inner(gs2,reference_multiplet[j])
        end
    end

    result = det(gamma_mat)

    if if_save
        which_gamma::Int = kwargs[:which_gamma]
        data_dict::Dict{String,Any} = Dict([("gamma$(which_gamma)",result)])
        filepath::String = kwargs[:filepath]
        save_hatsugai_data(data_dict,filepath,kwargs[:ref_multis_filename],which_gamma)
    end

    return result
end

# needs work
function get_omega(gs1::MPS,gs2::MPS,reference_multiplets::Vector{MPS}; kwargs...)
    if_save = get(kwargs,:if_save,false)

    omega_mat::Matrix{ComplexF64} = zeros(ComplexF64,2,2)
    for i in 1:2
        for j in 1:2
            #omega_mat[i,j] = (transpose(conj(reference_multiplets[i+2])) * gs1) * (transpose(conj(gs1)) * reference_multiplets[j]) + (transpose(conj(reference_multiplets[i+2])) * gs2) * (transpose(conj(gs2)) * reference_multiplets[j])
            omega_mat[i,j] = inner(reference_multiplets[i+2],gs1) * inner(gs1,reference_multiplets[j]) + inner(reference_multiplets[i+2],gs2) * inner(gs2,reference_multiplets[j])
        end
    end

    result = det(omega_mat)

    if if_save
        data_dict::Dict{String,Any} = Dict([("omega",result)])
        filepath::String = kwargs[:filepath]
        ref_multis_filenames::Vector{String} = kwargs[:ref_multis_filenames]
        save_hatsugai_data(data_dict,filepath,ref_multis_filenames)
    end

    return result
end

# needs work
function get_hatsugaifull(gs1::MPS,gs2::MPS,reference_multiplets::Vector{MPS}; kwargs...)
    if_save = get(kwargs,:if_save,false)

    gamma1::ComplexF64 = get_gamma(gs1,gs2,reference_multiplets[1:2]; if_save=false)
    gamma2::ComplexF64 = get_gamma(gs1,gs2,reference_multiplets[3:4]; if_save=false)
    omega::ComplexF64 = get_omega(gs1,gs2,reference_multiplets; if_save=false)

    if if_save
        data_dict::Dict{String,Any} = Dict([("gamma1",gamma1),("gamma2",gamma2),("omega",omega)])
        filepath::String = kwargs[:filepath]
        ref_multis_filenames::Vector{String} = kwargs[:ref_multis_filenames]
        save_hatsugai_data(data_dict,filepath,ref_multis_filenames)
    end

    return gamma1,gamma2,omega
end










































"fin"
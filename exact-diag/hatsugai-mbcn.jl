#####################################################
#=

This file is for calculating Hatsugai-style many body Chern number for ED

Depends on:
    other-funcs/data-storage-funcs.jl

=#
######################################################


function make_new_reference_multiplets(Lx::Int64,Ly::Int64,particles::Int64; kwargs...)
    tw1::Float64 = get(kwargs,:tw1,0.33)
    tw2::Float64 = get(kwargs,:tw2,0.67)
    hanis::Float64 = get(kwargs,:hopping_anisotropy,1.0)
    intstren::Float64 = get(kwargs,:interaction_strength,0.0)
    lr_dist::Int = get(kwargs,:lr,0)
    (intstren != 0.0 && lr_dist == 0) ? lr_dist = Ly-1 : nothing
    if_pinning::Bool = get(kwargs,:if_pinning,false)

    params_dict1 = Dict([("Lx",Lx),("Ly",Ly),("N",particles),("tw1",tw1),("tw2",tw2),("if_pinning",if_pinning),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis),("interaction_strength",intstren),("lr",lr_dist),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",true)])
    params_dict2 = Dict([("Lx",Lx),("Ly",Ly),("N",particles),("tw1",tw2),("tw2",tw1),("if_pinning",if_pinning),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis),("interaction_strength",intstren),("lr",lr_dist),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",true)])

    rez1 = run_normal_ed(params_dict1)
    rez2 = run_normal_ed(params_dict2)

    return [rez1[1][1],rez1[1][2],rez2[1][1],rez2[1][2]],rez1[4],rez2[4]
end

function get_reference_multiplets(Lx::Int64,Ly::Int64,particles::Int64; kwargs...)
    if_make_new::Bool = get(kwargs,:if_make_new,true)

    dataloc::String = get(kwargs,:dataloc,get_folder_location("cluster-data/exact-diag/torus"))
    tw1::Float64 = get(kwargs,:tw1,0.33)
    tw2::Float64 = get(kwargs,:tw2,0.67)
    hanis::Float64 = get(kwargs,:hopping_anisotropy,1.0)
    intstren::Float64 = get(kwargs,:interaction_strength,0.0)
    if_pinning::Bool = get(kwargs,:if_pinning,false)


    reference_multiplets::Vector{Vector{ComplexF64}} = [zeros(ComplexF64,2) for i in 1:4]

    params_dict1::Dict{String,Any} = Dict([("Lx",Lx),("Ly",Ly),("N",particles),("if_periodic_x",true),("if_periodic_y",true),("twist_angle1",tw1),("twist_angle2",tw2),("interaction_strength",intstren),("hopping_anisotropy",hanis),("if_pinning",if_pinning)])
    params_dict2::Dict{String,Any} = Dict([("Lx",Lx),("Ly",Ly),("N",particles),("if_periodic_x",true),("if_periodic_y",true),("twist_angle1",tw2),("twist_angle2",tw1),("interaction_strength",intstren),("hopping_anisotropy",hanis),("if_pinning",if_pinning)])

    if_exists1::Bool,found_data1::Union{Vector{Dict},Nothing} = check_data_exists(params_dict1,"ed"; location=dataloc,output_level=false)
    if_exists2::Bool,found_data2::Union{Vector{Dict},Nothing} = check_data_exists(params_dict2,"ed"; location=dataloc,output_level=false)

    if if_exists1 && if_exists2
        reference_multiplets[1:2] = found_data1[1]["state"][1:2]
        reference_multiplets[3:4] = found_data2[1]["state"][1:2]
    else
        if if_make_new
            println("Need to make new Hatsugai reference multiplets")
            return make_new_reference_multiplets(Lx,Ly,particles; tw1=tw1,tw2=tw2,hopping_anisotropy=hanis,interaction_strength=intstren,if_pinning=if_pinning)
        else
            error("Reference multiplets not found and if_make_new is false")
        end
    end

    return reference_multiplets,found_data1[2]["filename"],found_data2[2]["filename"]
end

function collate_hatsugai_data(Lx::Int64,Ly::Int64,N::Int64; kwargs...)
    intstren::Float64 = get(kwargs,:intstren,0.0)
    hanis::Float64 = get(kwargs,:hanis,1.0)
    if_pinning::Bool = get(kwargs,:if_pinning,false)

    dataloc::String = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis),("interaction_strength",intstren)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    #display(all_files)
    storing_file::String = ""

    tw1s::Vector{Float64} = Float64[]
    tw2s::Vector{Float64} = Float64[]
    omegas::Vector{ComplexF64} = ComplexF64[]
    lambda1s::Vector{ComplexF64} = ComplexF64[]
    lambda2s::Vector{ComplexF64} = ComplexF64[]
    for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        m["twist_angle"] == [0.0,0.0] && (storing_file = f)
        if haskey(m,"omega")
            append!(tw1s,m["twist_angle"][1])
            append!(tw2s,m["twist_angle"][2])
            append!(omegas,m["omega"])
            append!(lambda1s,m["gamma1"])
            append!(lambda2s,m["gamma2"])
        else
            println("No omega data found")
        end
    end
    display(tw1s)
    datadict = Dict([("tw1s",tw1s),("tw2s",tw2s),("omegas",omegas),("lambda1s",lambda1s),("lambda2s",lambda2s)])
    modify_data(datadict,dataloc * "/" * storing_file,"metadata"; output_level=0)

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

function get_gamma(gs1::Vector{ComplexF64},gs2::Vector{ComplexF64},reference_multiplet::Vector{Vector{ComplexF64}}; kwargs...)
    if_save = get(kwargs,:if_save,false)

    gamma_mat::Matrix{ComplexF64} = zeros(ComplexF64,2,2)
    for i in 1:2
        for j in 1:2
            gamma_mat[i,j] = (transpose(conj(reference_multiplet[i])) * gs1) * (transpose(conj(gs1)) * reference_multiplet[j]) + (transpose(conj(reference_multiplet[i])) * gs2) * (transpose(conj(gs2)) * reference_multiplet[j])
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

function get_omega(gs1::Vector{ComplexF64},gs2::Vector{ComplexF64},reference_multiplets::Vector{Vector{ComplexF64}}; kwargs...)
    if_save = get(kwargs,:if_save,false)

    omega_mat::Matrix{ComplexF64} = zeros(ComplexF64,2,2)
    for i in 1:2
        for j in 1:2
            omega_mat[i,j] = (transpose(conj(reference_multiplets[i+2])) * gs1) * (transpose(conj(gs1)) * reference_multiplets[j]) + (transpose(conj(reference_multiplets[i+2])) * gs2) * (transpose(conj(gs2)) * reference_multiplets[j])
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

function get_hatsugaifull(gs1::Vector{ComplexF64},gs2::Vector{ComplexF64},reference_multiplets::Vector{Vector{ComplexF64}}; kwargs...)
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
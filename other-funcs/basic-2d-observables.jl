#####################################################
#=

This file contains useful observable functions for 2D physics

=#
######################################################

function get_cdwsf(qvec::Vector{Float64},dens_corr_mat::Array{Float64}; kwargs...)
    Lphys::Int64,Lsynth::Int64 = size(dens_corr_mat)[1],size(dens_corr_mat)[3]

    result::Union{Float64,ComplexF64} = 0.0
    
    for j in 1:Lphys
        for s in 1:Lsynth
            for jj in 1:Lphys
                for ss in 1:Lsynth
                    dist_vect::Vector{Int64} = [minimum([(j-jj),Lphys-(j-jj)]),minimum([(s-ss),Lsynth-(s-ss)])]
                    dotprod::Float64 = dot(qvec,dist_vect)
                    result += exp(-im*2*pi*dotprod)*dens_corr_mat[j,jj,s,ss]
                end
            end
        end
    end

    return result / ((Lphys*Lsynth)^2)
end

function get_cdwsf(angle::Float64,dens_corr_mat::Array{Float64},radius::Int64=1; kwargs...)
    qvec::Vector{Float64} = [radius*cos(pi*angle),radius*sin(pi*angle)]
    return get_cdwsf(qvec,dens_corr_mat; kwargs...)
end

function save_ft_dd(ft_dd_val::ComplexF64,angle::Float64,filepath::String)
    data_dict = Dict([("ft_dd_$angle",ft_dd_val)])
    modify_data_jld2(data_dict,filepath,"metadata"; output_level=1)
end


































"fin"
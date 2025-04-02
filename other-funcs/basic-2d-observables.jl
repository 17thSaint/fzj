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

function save_dd_correlation(dens_corr_mat::Array{Float64},filepath::String)
    data_dict::Dict{String,Array} = Dict([("dens_corr_mat",dens_corr_mat)])
    modify_data_jld2(data_dict,filepath,"metadata"; output_level=1)
end

function integrate_2d_matrix(matrix::Matrix{Float64})
    """
    Performs 2D integration over the values in a 2D matrix.
    
    Parameters:
    matrix (Array{Float64, 2}): The 2D matrix to integrate over.
    
    Returns:
    Float64: The result of the 2D integration.
    """
    
    # Assuming the matrix has dimensions m x n
    m::Int64, n::Int64 = size(matrix)
    
    # Initialize the integral
    integral::Float64 = 0.0
    
    # Perform the 2D integration using trapezoidal rule
    for i in 1:(m-1)
        for j in 1:(n-1)
            integral += (matrix[i,j] + matrix[i,j+1] + matrix[i+1,j] + matrix[i+1,j+1]) / 4
        end
    end
    
    # Account for the step sizes in the x and y directions
    dx::Float64 = 1 / (m - 1)
    dy::Float64 = 1 / (n - 1)
    integral *= dx * dy
    
    return integral
end

function ft_density(momentum::Vector{Float64},occs::Matrix{Float64}; kwargs...)
    ly,lx = size(occs)

    all_positions::Array{Float64,2} = zeros(Float64,ly,lx)
    for j in 1:lx
        for s in 1:ly
            all_positions[s,j] = dot(momentum,[j,s])
        end
    end

    result::ComplexF64 = sum(occs .* exp.(im .* all_positions)) / (lx * ly)

    return result
end

function ft_densitydensity(momentum::Vector{Float64},dds::Matrix{Float64}; kwargs...)
    ly,lx = size(dds)

    center_site::Vector{Int64} = [Int64(ceil(lx/2)),Int64(ceil(ly/2))]

    all_distances::Array{Float64,2} = zeros(Float64,ly,lx)
    for j in 1:lx
        for s in 1:ly
            all_distances[s,j] = dot(momentum,([j,s] - center_site))
        end
    end

    result::ComplexF64 = sum(dds .* exp.(im .* all_distances)) / (lx * ly)

    return result
end

function pairdistribution(densitydensity::Matrix{Float64},occupancy_matrix; kwargs...)
    if_plot::Bool = get(kwargs, :if_plot, false)

    if size(densitydensity) != size(occupancy_matrix)
        println("Sizes of Occs and DDs are not the same. Trimming Occs.")
        occupancy_matrix = occupancy_matrix[1:size(densitydensity)[1],1:size(densitydensity)[2]]
    end

    Lsynth::Int64,Lphys::Int64 = size(densitydensity)

    centersite::Vector{Int64} = [Int64(ceil(Lphys/2)),Int64(ceil(Lsynth/2))]
    pairdist::Matrix{Float64} = densitydensity ./ (occupancy_matrix[centersite[2],centersite[1]] .* occupancy_matrix)    

    if_plot ? plot_pairdistribution(pairdist; kwargs...) : nothing

    return pairdist
end

function entanglement_entropy(spec::Vector{Float64}; kwargs...)
    return -sum((spec .^ 1) .* ( 1 .* log.(spec)))
end

function pairdist_ellipticalness(pairdist::Matrix{Float64}; kwargs...)

    centersite = [Int64(ceil(size(pairdist,2)/2)),Int64(ceil(size(pairdist,1)/2))]
    prob_pairdist = 1 .- pairdist
    prob_pairdist = prob_pairdist ./ sum(prob_pairdist)
    x_avg::Float64 = 0.0
    y_avg::Float64 = 0.0
    x_var::Float64 = 0.0
    y_var::Float64 = 0.0
    xy_var::Float64 = 0.0
    for i in 1:size(prob_pairdist,2)
        for j in 1:size(prob_pairdist,1)
            x_avg += (i - centersite[1]) * prob_pairdist[j,i]
            y_avg += (j - centersite[2]) * prob_pairdist[j,i]
            x_var += (i - centersite[1])^2 * prob_pairdist[j,i]
            y_var += (j - centersite[2])^2 * prob_pairdist[j,i]
            xy_var += (i - centersite[1]) * (j - centersite[2]) * prob_pairdist[j,i]
        end
    end

    ellipticalness::Float64 = x_var / y_var

    return ellipticalness,x_avg,y_avg,x_var,y_var,xy_var
end

function ft_coeff(phys_site::Tuple{Int,Int},momentum::Vector{Float64},op_type::String; kwargs...)
    dag_sign::Int = op_type == "Adag" ? 1 : -1
    return exp(2*pi*im*dag_sign*dot(momentum,phys_site))
end

function diocane(phys_site::Tuple,momentum::Vector,op_type::String; kwargs...)
    Ly::Int = kwargs[:Ly]
    mval = Int(momentum[2] * Ly)

    if phys_site[1]-1 == mval
        dag_sign::Int = op_type == "Adag" ? 1 : -1
        return exp(dag_sign * 2*pi*im*momentum[2]*(phys_site[2]-1)) / sqrt(Ly)
    else
        return 0.0
    end
end

































"fin"
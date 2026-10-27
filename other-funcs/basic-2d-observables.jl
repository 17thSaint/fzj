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

function get_cdwsf_occs(qvec::Vector{Float64},occs::Matrix{Float64}; kwargs...)
    Lphys::Int64,Lsynth::Int64 = size(occs)

    result::Union{Float64,ComplexF64} = 0.0
    
    for j in 1:Lphys
        for s in 1:Lsynth
            dotprod::Float64 = dot(qvec,[j,s])
            result += exp(-im*2*pi*dotprod)*occs[j,s]
        end
    end

    return result
end

function get_cdwsf_occs(angle::Float64,occs::Matrix{Float64},radius::Int64=1; kwargs...)
    qvec::Vector{Float64} = [radius*cos(pi*angle),radius*sin(pi*angle)]
    return get_cdwsf_occs(qvec,occs; kwargs...)
end

# column statistics for the CDW order parameters from configurations i with probability probs[i] (normalized here)
# and column occupations ncols[x,i] = sum_y n_{x,y} (x physical), see column_occupations in exact-diag/observables.jl
# (ED basis) and fpeps/observables.jl (samples): column_corr[x,x'] = <n_x n_x'> and evencolumn_dist[k+1] =
# P(k particles on even columns), k = 0..N. Small enough to save in metadata and correlation based, so works for
# translation invariant states
function column_cdw_statistics(ncols::AbstractMatrix{<:Integer},probs::AbstractVector{Float64})
    Lx::Int64 = size(ncols,1)
    N::Int64 = sum(view(ncols,:,1))
    all(sum(ncols,dims=1) .== N) || error("configurations have different particle numbers")

    w = probs ./ sum(probs)
    ncols_float = Float64.(ncols)
    column_corr = (ncols_float .* transpose(w)) * transpose(ncols_float)

    evencolumn_dist = zeros(Float64,N+1)
    for (i,k) in enumerate(sum(view(ncols,2:2:Lx,:),dims=1))
        evencolumn_dist[k+1] += w[i]
    end

    return column_corr, evencolumn_dist
end

# CDW order parameters from column_cdw_statistics (pass averages over a manifold to make them basis independent):
# S(pi) = <|sum_x e^{i pi x} n_x|^2> / N, correlation ratio 1 - S(pi - 2pi/Lx)/S(pi), and Binder cumulant
# 1 - <m^4>/(3<m^2>^2) of the staggered column magnetization m = sum_x (-1)^x n_x / N
function cdw_order_parameters(column_corr::Matrix{Float64},evencolumn_dist::Vector{Float64})
    Lx::Int64 = size(column_corr,1)
    N::Int64 = length(evencolumn_dist) - 1

    structure_factor(q) = real(sum(exp(im*q*(x-xx)) * column_corr[x,xx] for x in 1:Lx, xx in 1:Lx)) / N
    s_pi = structure_factor(pi)
    corr_ratio = 1 - structure_factor(pi - 2pi/Lx) / s_pi

    stag_mags = [(2k - N) / N for k in 0:N]
    m2,m4 = dot(evencolumn_dist,stag_mags.^2),dot(evencolumn_dist,stag_mags.^4)
    binder = 1 - m4 / (3*m2^2)

    return s_pi, corr_ratio, binder
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
        return 0.00000000000000001 + im*0.0
    end
end

function visibility_fourpt(four_point_moms::Matrix{Float64}; kwargs...)
    @assert size(four_point_moms,2) > 6 "Lattice must be larger than 6 to see visibility"

    restricted_fourpts = four_point_moms[1,4:end-2]

    if maximum(restricted_fourpts) < 1e-7
        restricted_fourpts = four_point_moms[2,4:end-2]
    end

    maxval = maximum(restricted_fourpts)
    minval = minimum(restricted_fourpts)
    maxminratio = minval / maxval

    return maxminratio
end

# gauge-invariant particle current on the bond s -> s + e_dir (dir = 1: x/phys, 2: y/synth), from the density matrix
# rho[i,j] = <b^dag_i b_j>. For the hopping term coeff b^dag_s b_{s+e} + h.c. in H, J = d<n_{s+e}>/dt = -2 Im(coeff rho[s,s+e]).
# hopping_coeff(site,dir) is coeff in the gauge of the code that built H, kwarg site_index(site) the row of rho
# (default linear_index on the Lx x Ly lattice)
function bond_current(densmat::Matrix{ComplexF64},hopping_coeff::Function,Lx::Int64,Ly::Int64,site::Tuple{Int64,Int64},dir::Int64; kwargs...)
    site_index = get(kwargs,:site_index,s -> linear_index(s,Lx,Ly))

    next_site = dir == 1 ? (mod1(site[1]+1,Lx),site[2]) : (site[1],mod1(site[2]+1,Ly))

    return -2 * imag(hopping_coeff(site,dir) * densmat[site_index(site),site_index(next_site)])
end

# currents[x,y] on the bond (x,y) -> (x+1,y), (Lx-1) x Ly for open x
function physical_current(densmat::Matrix{ComplexF64},hopping_coeff::Function,Lx::Int64,Ly::Int64,if_periodic_x::Bool; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    xlen = if_periodic_x ? Lx : Lx-1
    currents = [bond_current(densmat,hopping_coeff,Lx,Ly,(x,y),1; kwargs...) for x in 1:xlen, y in 1:Ly]

    if_plot ? plot_physical_current(currents; kwargs...) : nothing

    return currents
end

# currents[x,y] on the bond (x,y) -> (x,y+1), Lx x (Ly-1) for open y
function synthetic_current(densmat::Matrix{ComplexF64},hopping_coeff::Function,Lx::Int64,Ly::Int64,if_periodic_y::Bool; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    ylen = if_periodic_y ? Ly : Ly-1
    currents = [bond_current(densmat,hopping_coeff,Lx,Ly,(x,y),2; kwargs...) for x in 1:Lx, y in 1:ylen]

    if_plot ? plot_synthetic_current(currents; kwargs...) : nothing

    return currents
end

































"fin"
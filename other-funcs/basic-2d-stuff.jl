#####################################################
#=

This file contains basic useful functions for 2D physics

=#
######################################################

using LinearAlgebra


function find_dist(p1::Tuple{Int,Int}, p2::Tuple{Int,Int}, size::Tuple{Int,Int}, periodic::Tuple{Bool,Bool}=(false, false))
    dx::Float64 = abs(p1[1] - p2[1])
    dy::Float64 = abs(p1[2] - p2[2])

    if periodic[1]
        dx = min(dx, size[1] - dx)
    end

    if periodic[2]
        dy = min(dy, size[2] - dy)
    end

    return sqrt(dx^2 + dy^2),(dx,dy),(p1[1]-p2[1],p1[2]-p2[2])
end

function bin_values(vector::Vector{T}, num_bins::Int) where T
    min_val::T = minimum(vector)
    max_val::T = maximum(vector)
    bin_width::Float64 = (max_val - min_val) / num_bins

    bins::Dict{T,Int} = Dict()
    for val in vector
        bin_index::Int = min(floor(Int, (val - min_val) / bin_width) + 1, num_bins)
        bins[val] = bin_index
    end

    return bins
end

# checks integer number of fluxes and then if periodic checks maintaining magnetic translations
function check_fluxes(alpha::Float64,Lx::Int64,Ly::Int64,if_periodic_x::Bool,if_periodic_y::Bool,flux_direction::String; kwargs...)
    if_error = get(kwargs,:if_error,true)
    if_ed = get(kwargs,:if_ed,true)
    output_level = kwargs[:output_level]

    # needs to be general to ttn and ed
    if if_ed
        if flux_direction == "x"
            flux_direction = "phys"
        elseif flux_direction == "y"
            flux_direction = "synth"
        end
    end

    if alpha == 0.0
        return nothing
    end
	if alpha > 0.4
        error("Alpha is too large: ",alpha)
    end
    x_shift,y_shift = !if_periodic_x, !if_periodic_y
    num_fluxes::Float64 = round(alpha*(Lx - x_shift) * (Ly - y_shift),digits=5)
    if_error && output_level > 0 ? println("Number of Fluxes = ",num_fluxes," for Lx = ",Lx," and Ly = ",Ly) : nothing
    if !isinteger(num_fluxes)
        if_error ? error("Number of fluxes is not an integer") : return false
    end

	if_error && output_level > 0 ? println("Checking fluxes only along Gauge Direction") : nothing
    if flux_direction == "synth"
        if if_periodic_x && !isinteger(num_fluxes/Ly)
            if if_periodic_y && isinteger(num_fluxes/Lx)
                flux_direction = "phys"
                if_error && output_level > 0 ? println("Fluxes don't fit, changing to X direction") : nothing
            else
                if_error ? error("Number of fluxes is not an integer multiple of Lx") : return false
            end
        end
    elseif flux_direction == "phys"
        if if_periodic_y && !isinteger(num_fluxes/Lx)
            if if_periodic_x && isinteger(num_fluxes/Ly)
                flux_direction = "synth"
                if_error && output_level > 0 ? println("Fluxes don't fit, changing to Y direction") : nothing
            else
                if_error ? error("Number of fluxes is not an integer multiple of Ly") : return false
            end
        end
    else
        error("Flux direction is not valid")
    end


	#=
    if if_periodic_x && !isinteger(num_fluxes/Lx)
        error("Number of fluxes is not an integer multiple of Lx")
    end

    if if_periodic_y && !isinteger(num_fluxes/Ly)
        error("Number of fluxes is not an integer multiple of Ly")
    end=#

    # convert flux_direction back to x or y if ED
    if if_ed
        if flux_direction == "phys"
            flux_direction = "x"
        elseif flux_direction == "synth"
            flux_direction = "y"
        end
    end

    return flux_direction
end























"fin"
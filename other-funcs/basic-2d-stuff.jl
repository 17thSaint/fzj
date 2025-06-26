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

# checks flux counting for fluxes in both directions
function check_fluxes(alpha::Vector{Float64},Lx::Int64,Ly::Int64,if_periodic_x::Bool,if_periodic_y::Bool,flux_direction::Vector{String}; kwargs...)
    if_error = get(kwargs,:if_error,true)
    if_ed = get(kwargs,:if_ed,true)
    output_level = kwargs[:output_level]

    # needs to be general to ttn and ed
    if if_ed
        if flux_direction[1] == "x"
            flux_direction[1] = "phys"
        elseif flux_direction[2] == "y"
            flux_direction[2] = "synth"
        end
    end

    total_alpha::Float64 = alpha[2] - alpha[1]

    if alpha[1] == 0.0 && alpha[2] == 0.0
        return nothing
    end
	if total_alpha > 0.4
        error("Alpha is too large: ",alpha)
    end
    x_shift,y_shift = !if_periodic_x, !if_periodic_y
    total_num_fluxes::Float64 = round(total_alpha*(Lx - x_shift) * (Ly - y_shift),digits=5)
    if_error && output_level > 0 ? println("Total Number of Fluxes = ",total_num_fluxes," for Lx = ",Lx," and Ly = ",Ly) : nothing
    if !isinteger(total_num_fluxes)
        if_error ? error("Total Number of fluxes is not an integer") : return false
    end

    results = ["phys","synth"]
    for i in 1:2
        num_fluxes::Float64 = round(alpha[i]*(Lx - x_shift) * (Ly - y_shift),digits=5)
        if_error && output_level > 0 ? println("Checking fluxes only along Gauge Direction") : nothing
        if flux_direction[1] == "synth"

            if if_periodic_x && !isinteger(num_fluxes/Ly)
                if_error ? error("Number of fluxes is not an integer multiple of Lx") : return false
            end

        elseif flux_direction[1] == "phys"

            if if_periodic_y && !isinteger(num_fluxes/Lx)
                if_error ? error("Number of fluxes is not an integer multiple of Ly") : return false
            end

        else
            error("Flux direction is not valid")
        end

    end

    # convert flux_direction back to x or y if ED
    if if_ed
        if results[1] == "phys"
            results[1] = "x"
        elseif results == "synth"
            results = "y"
        end
    end

    return results
end

# finds the linear index assuming jump snake mapping with site 1 at bottom left corner
function linear_index(site::Tuple{Int64,Int64},Lx::Int64,Ly::Int64)
    return (site[2] - 1)*Lx + site[1]
end

function linear_index(site::Vector{Int64},Lx::Int64,Ly::Int64)
    return (site[2] - 1)*Lx + site[1]
end

# finds the site assuming jump snake mapping with site 1 at bottom left corner
function coordinate(site::Int64,Lx::Int64,Ly::Int64)
    x = mod1(site,Lx)
    y = Int((site - x) / Lx + 1)
    return (x,y)
end

# find how many Binary TTN layers need for given lattice size
function get_layers_from_latticesize(Lx::Int64,Ly::Int64)
    if Lx == Ly || Lx == 2*Ly || Ly == 2*Lx
        layers = log(2,Lx*Ly)
        if isinteger(layers)
            return Int(layers)
        else
            error("Lattice size does not fit in Binary Tree")
        end
    else
        error("Lattice size does not fit in Binary Tree")
    end
end

# old version from TTN code
#=function long_range_scaling(x_final,virt_edge_length,initial_strength; kwargs...)
	if_plot = get(kwargs, :if_plot, false)
	if_save_data = get(kwargs, :if_save_data, false)
	if_save_fig = get(kwargs, :if_save_fig, false)
	if_hard_cutoff = get(kwargs, :cliff, false)
	if_rounding = get(kwargs, :rounding, true)
	if if_hard_cutoff
		if_rounding = false
	end
	final_minimum = get(kwargs, :limit, 10^-3)
	trunc_minimum = get(kwargs, :trunc_min, 10^-6)
	trunc = get(kwargs, :trunc, trunc_minimum*initial_strength)
	scaling_func = get(kwargs, :scaling, "flat")
	
	strengths = zeros(virt_edge_length)
	
	if scaling_func == "flat"
		strengths[1:x_final+1] .= initial_strength
	elseif scaling_func == "exp"
		strengths = map(1:virt_edge_length) do x
			initial_strength * exp(-log(1/final_minimum)*(x-1)/x_final)	
		end
		strengths[1] = initial_strength
	elseif scaling_func == "lr_flat"
		strengths[1] = initial_strength
		strengths[2:x_final+1] .= final_minimum
	elseif scaling_func == "rydberg"
		blockade_radius = initial_strength
		strengths = map(0:virt_edge_length-1) do x
			1.0 * (blockade_radius^6) / (blockade_radius^6 + x^6)
		end
		x_final < length(strengths) ? strengths[x_final+2:end] .= 0.0 : nothing
	end
	
	if if_hard_cutoff
		strengths[x_final + 2:end] .= 0.0
	elseif if_rounding
		final_index = findfirst(x -> abs(x) .<= trunc,strengths)
		if !isnothing(final_index)
			strengths[final_index:end] .= 0.0
		end
	end
	
	if_plot || if_save_fig ? plot_long_range_scaling(strengths,virt_edge_length; kwargs...) : nothing
	#if_save_data ? save_long_range_scaling(strengths,virt_edge_length; kwargs...) : nothing

	return strengths
end=#

function long_range_scaling(x_final::Int64,virt_edge_length::Int64,initial_strength::Float64; kwargs...)

	trunc = get(kwargs, :trunc_digits, 5)
	scaling_func = get(kwargs, :scaling, "flat")
	
	strengths = zeros(virt_edge_length)

    if x_final == 0.0 || initial_strength == 0.0
        strengths[1] = initial_strength != 0.0 ? initial_strength : 1.0
        return strengths
    end
	
	if scaling_func == "flat"
		strengths[1:x_final+1] .= initial_strength
    elseif scaling_func == "gaussian"
        sigma = get(kwargs, :sigma, 1.0)
        strengths = map(0:virt_edge_length-1) do x
            #(initial_strength/(sqrt(2*pi)*sigma)) * exp(-x^2/(2*sigma^2))
            initial_strength * exp(-x^2/(2*sigma^2))
        end
	elseif scaling_func == "exp"
        corr_length = get(kwargs, :corr_length, virt_edge_length)
		strengths = map(1:virt_edge_length) do x
			initial_strength * exp(-(x-1)/corr_length)	
		end
	elseif scaling_func == "rydberg"
		blockade_radius = get(kwargs, :blockade_radius, 1.0)
		strengths = map(0:virt_edge_length-1) do x
			initial_strength * (blockade_radius^6) / (blockade_radius^6 + x^6)
		end
	end
	
	strengths = round.(strengths,digits=trunc)

	return strengths
end

function get_corr_length_from_Us(us::Vector{Float64})
    init_stren = us[1]
    all_cs = -1 ./ (log.(us[2:end] ./ init_stren) ./ collect(1:length(us)-1))
    return sum(all_cs) / (length(all_cs))
end


function get_hopping_strengths(t_strength::Float64,hopping_anisotropy::Float64)
	if hopping_anisotropy < 1.0
		t_strength_synth::Float64 = t_strength / hopping_anisotropy
		t_strength_phys::Float64 = t_strength
	else
		t_strength_phys = t_strength * hopping_anisotropy
		t_strength_synth = t_strength
	end

    #t_strength_phys = hopping_anisotropy
    #t_strength_synth = t_strength

	return t_strength_phys,t_strength_synth
end





















"fin"
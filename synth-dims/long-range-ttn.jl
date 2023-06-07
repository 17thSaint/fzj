using LinearAlgebra,PyPlot
include("../review-practice-codes/ttn.jl")

function long_range_scaling(x_final,edge_length,initial_strength; kwargs...)
	if_hard_cutoff = get(kwargs, :cliff, false)
	if_rounding = get(kwargs, :rounding, true)
	final_minimum = get(kwargs, :limit, 10^-3)
	trunc_minimum = get(kwargs, :trunc_min, 10^-6)
	trunc = get(kwargs, :trunc, trunc_minimum*initial_strength)
	scaling_func = get(kwargs, :scaling, "flat")
	
	strengths = zeros(edge_length)
	
	if scaling_func == "flat"
		strengths[1:x_final+1] .= initial_strength
	elseif scaling_func == "exp"
		strengths = map(1:edge_length) do x
			initial_strength * exp(-log(1/final_minimum)*(x-1)/x_final)	
		end
	elseif scaling_func == "power"
		println("Not done power")
	end
	
	if if_hard_cutoff
		strengths[x_final + 1:end] .= 0.0
	elseif if_rounding
		final_index = findfirst(x -> x .<= trunc,strengths)
		if !isnothing(final_index)
			strengths[final_index:end] .= 0.0
		end
	end

	return strengths
end































"fin"

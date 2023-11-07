function start_rand_config(num_parts::Int, m::Int)
    # Calculate the filling fraction
    filling = 1 / m

    # Calculate the characteristic length scale
    rm = sqrt(2 * num_parts / filling)

    # Generate random real and imaginary parts in one step
    real_parts = rand(Float64, num_parts) .* rand(-1:2:1, num_parts) .* rm
    imag_parts = rand(Float64, num_parts) .* rand(-1:2:1, num_parts) .* rm

    # Combine real and imaginary parts to create complex numbers
    config = real_parts .- im .* imag_parts

    return config
end

function jastrow(z,which; kwargs...)	
	if_log = get(kwargs, :if_log, false)
	
	if if_log
		result = 0.0*im
		for i in 1:length(z)
			if z[i] != which
				result += log(Complex(which - z[i]))
			end
		end
	else
		result = 1.0
		for i in 1:length(z)
			if z[i] != which
				result *= (which - z[i])^1
			end
		end
	end
	
	return result
end

function jastrow_squared(z,which; kwargs...)
	if_log = get(kwargs, :if_log, false)
		
	if if_log
		result = 0.0*im
		for i in 1:length(z)
			if z[i] != which
				result += 2*log(Complex(which - z[i]))
			end
		end
	else
		result = 1.0
		for i in 1:length(z)
			if z[i] != which
				result *= (which - z[i])^2
			end
		end
	end
	return result
end

function get_log_add(a,b)
	if real(a) > real(b)
		ordered::Vector{typeof(a)} = [b,a]
	else
		ordered = [a,b]
	end
	result::ComplexF64 = ordered[2] + log(Complex(1 + exp(ordered[1] - ordered[2])))
	return Complex(result)
end

function get_log_subtract(a,b)
	ordered = [a,b]
	result::ComplexF64 = ordered[1] + log(Complex(1 - exp(ordered[2] - ordered[1])))
	return Complex(result)
end

function get_log_det(matrix::Matrix{ComplexF64},reg_input=false)
	num_parts::Int64 = size(matrix)[1]
	maxes::Vector{ComplexF64} = [0.0+0.0*im for i in 1:num_parts]
	rejected_indices::Vector{Int64} = []
	
	changed::Matrix{ComplexF64} = matrix + fill(0.0,(num_parts,num_parts))
	for i in 1:num_parts
		row::Vector{Float64} = real(changed[i,:])
		validation::Bool = true
		j::Int64 = 0
		index::Int64 = 0
		while validation
			val::ComplexF64 = sort(row)[end - j]
			guess_index::Int64 = findfirst(real.(changed[i,:]) .== val )
			if any(rejected_indices .== index)
				j += 1
			else
				index = guess_index
				validation = false
			end
		end
		maxes[i] = changed[i,index]
		changed[i,:] .-= maxes[i]
		append!(rejected_indices,index)
	end
	reduced_logdet::ComplexF64 = sum(maxes) + log(Complex(det(exp.(changed))))

	return reduced_logdet,maxes,changed
end


function dist_btw(z,i,j,pow)
	return (z[i] - z[j])^pow
end


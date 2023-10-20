using SpecialMatrices

function trace_part(config,which_particle=1)
	vdm_deriv = [(i-1)*config[which_particle]^(i-2) for i in 1:length(config)]
	return ( conj(Vandermonde(config)') \ vdm_deriv )[which_particle]
	#return  -1 * det(Vandermonde(config)) * (( conj(Vandermonde(config)') \ vdm_deriv )[which_particle])
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

function reverse_flux_wavefunction(z,m=3)
	p = Int((m+1)/2)
	num_parts = length(z)
	nt = sum([i for i in 1:num_parts-1])
	wavefunc = 0.0*im
	
	c_start = time()
	vdm_logdet = get_log_det(log.(Vandermonde(z)))[1]
	const_part = nt * (log(4*p) + 2*p*vdm_logdet)#log(Complex(det(Vandermonde(z)))))
	c_end = time()
	
	d_start = time()
	deriv_part = 0.0*im
	#all_derivs = [deriv_of_slater(z,i) for i in 1:num_parts]
	for j in 1:num_parts
		for i in 1:j-1
			new_term = trace_part(z,i) - trace_part(z,j)
			deriv_part += log(Complex(new_term))
		end
	end
	d_end = time()
	
	log_exponent = -0.25 * sum(abs2.(z))
	
	times = [100*(c_end-c_start),100*(d_end-d_start)]
	
	#println("Deriv $deriv_part, Const $const_part")
	
	return deriv_part + const_part + log_exponent,times
end

function dist_btw(z,i,j,pow)
	return (z[i] - z[j])^pow
end

function test_3parts(z,m=3)
	num_parts = length(z)
	p = Int((m+1)/2)
	
	dz = im .* zeros(3)
	dz[1] = 2*p*dist_btw(z,2,3,2*p) * (dist_btw(z,1,2,2*p-1)*dist_btw(z,1,3,2*p) + dist_btw(z,1,2,2*p)*dist_btw(z,1,3,2*p-1))
	dz[2] = 2*p*dist_btw(z,1,3,2*p) * (-dist_btw(z,1,2,2*p-1)*dist_btw(z,2,3,2*p) + dist_btw(z,1,2,2*p)*dist_btw(z,2,3,2*p-1))
	dz[3] = 2*p*dist_btw(z,1,2,2*p) * (-dist_btw(z,1,3,2*p-1)*dist_btw(z,2,3,2*p) - dist_btw(z,1,3,2*p)*dist_btw(z,2,3,2*p-1))
	
	result = 1.0
	for j in 1:num_parts
		for i in 1:j-1
			result *= -2*(dz[i] - dz[j])
		end
	end
	
	for i in 1:num_parts
		result *= exp(-0.25*abs2(z[i]))
	end
	
	return result
end








































"fin"

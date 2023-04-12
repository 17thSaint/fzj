using ITensors,LinearAlgebra,PyPlot,Statistics,SparseArrays

x = [0 1;1 0]
xx = [0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0]
z = [1 0;0 -1]
zz = [1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1]
id = [1 0;0 1]
zers = [0 0;0 0]

function get_flattened_index(b_list)
	return sum(b_list .* [2^(length(b_list) - i) for i in 1:length(b_list)]) + 1
end

function int_to_binary(given::Int,len::Int)
    # Initialize an empty array to store the binary digits
    binary = []
    # Loop until the integer is reduced to zero
    while given > 0
        # Get the remainder of n when divided by 2
        digit = given % 2
        # Prepend the digit to the binary array
        pushfirst!(binary, digit)
        # Integer divide n by 2 to reduce its value
        given = div(given, 2)
    end
    # Return the binary array
    if length(binary) < len
    	return append!([0 for i in 1:(len-length(binary))],binary)
    else
    	return binary
    end
end


function get_single_qubit_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	mat,which_qubits = get(kwargs, :arguments, 3)
	prod_parts = []
	for j in 1:length(bs)
		if j != which_qubits
			if bs[j] != bps[j]
				return 0.0
			else
				append!(prod_parts,[bs[j] == bps[j]])
			end
		end
	end
	if length(prod_parts) < 1
		prod_rez = 1
	else
		prod_rez = prod(prod_parts)
	end
	return prod_rez * mat[bs[which_qubits]+1,bps[which_qubits]+1] 
end

function get_exp_single_qubit_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	mat,which_qubits,strength,dt = get(kwargs, :arguments, 3)
	coeff = -im * strength * dt / 2
	exp_mat = exp(coeff.*mat)
	return get_single_qubit_elem(bs=bs,bps=bps,arguments=(exp_mat,which_qubits))
end

function get_two_qubit_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	mat,which_qubits = get(kwargs, :arguments, 3)
	prod_parts = []
	for j in 1:length(bs)
		if j != which_qubits[1] && j != which_qubits[2]
			if bs[j] != bps[j]
				return 0.0
			else
				append!(prod_parts,[bs[j] == bps[j]])
			end
		end
	end
	if length(prod_parts) < 1
		prod_rez = 1
	else
		prod_rez = prod(prod_parts)
	end
	mat_b_elem_index = get_flattened_index([bs[which_qubits[1]],bs[which_qubits[2]]])
	mat_bp_elem_index = get_flattened_index([bps[which_qubits[1]],bps[which_qubits[2]]])
	return prod_rez * mat[mat_b_elem_index,mat_bp_elem_index]
end

function get_exp_two_qubit_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	mat,which_qubits,strength,dt = get(kwargs, :arguments, 3)
	coeff = -im * strength * dt
	exp_mat = exp(coeff.*mat)
	return get_two_qubit_elem(bs=bs,bps=bps,arguments=(exp_mat,which_qubits))
end

function get_mat_type(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	mat,which_qubits,strength,dt = get(kwargs, :arguments, 3)
	if size(mat)[1] > 2
		return get_exp_two_qubit_elem(bs=bs,bps=bps,arguments=(mat,which_qubits,strength,dt))
	else
		return get_exp_single_qubit_elem(bs=bs,bps=bps,arguments=(mat,which_qubits,strength,dt))
	end
end

#=
function elem_multiply_exp_matrices(mat1,mat2,bs,bps,which_qubits1,which_qubits2,strength1,strength2,dt)
	site_count = length(bs)
	final_val = im*0.0
	for i in 1:2^site_count
		summed_index = int_to_binary(i-1,site_count)
		
		left = get_mat_type(mat1,which_qubits1,bs,summed_index,strength1,dt)
		right = get_mat_type(mat2,which_qubits2,summed_index,bps,strength2,dt)
		
		final_val += left * right
	end
	return final_val
end
=#
function build_matrix_from_elements(func,arguments,site_count)
	mat = im.*zeros(2^site_count,2^site_count)
	for i in 1:2^site_count
		for j in 1:2^site_count
			bs = int_to_binary(i-1,site_count)
			bps = int_to_binary(j-1,site_count)
			mat[i,j] = func(bs=bs,bps=bps,arguments=arguments)
		end
	end
	return mat
end

function elem_mult_matrices(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	left_func,left_args,right_func,right_args = get(kwargs, :arguments, 3)
	site_count = length(bs)
	final_val = im*0.0
	for i in 1:2^site_count
		summed_index = int_to_binary(i-1,site_count)

		left = left_func(bs=bs,bps=summed_index,arguments=left_args)

		right = right_func(bs=summed_index,bps=bps,arguments=right_args)
		
		final_val += left * right
	end
	return final_val
end

function get_zz_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	which_qubits,hz_strength,dt = get(kwargs, :arguments, 3)
	left_args = (z,which_qubits[1],hz_strength,dt)
	right_args = (z,which_qubits[2],hz_strength,dt)
	return elem_mult_matrices(bs=bs,bps=bps,arguments=(get_exp_single_qubit_elem,left_args,get_exp_single_qubit_elem,right_args))
end

function get_xx_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	which_qubits,j_strength,dt = get(kwargs, :arguments, 3)
	local_args = (xx,which_qubits,j_strength,dt)
	return get_exp_two_qubit_elem(bs=bs,bps=bps,arguments=local_args)
end

function get_x_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	which_qubit,hx_strength,dt = get(kwargs, :arguments, 3)
	local_args = (x,which_qubit,hx_strength,dt)
	return get_exp_single_qubit_elem(bs=bs,bps=bps,arguments=local_args)
end

function get_fullz_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	which_qubits,j_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	return elem_mult_matrices(bs=bs,bps=bps,arguments=(get_xx_elem,(which_qubits,j_strength,dt),get_zz_elem,(which_qubits,hz_strength,dt)))
end

function get_left_localham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	which_qubits,j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	return elem_mult_matrices(bs=bs,bps=bps,arguments=(get_x_elem,(which_qubits[1],hx_strength,dt),get_fullz_elem,(which_qubits,j_strength,hz_strength,dt)))
end

function get_localham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	which_qubits,j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	return elem_mult_matrices(bs=bs,bps=bps,arguments=(get_left_localham_elem,(which_qubits,j_strength,hx_strength,hz_strength,dt),get_x_elem,(which_qubits[1],hx_strength,dt)))
end

function get_3ham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	left_args = ([1,2],j_strength,hx_strength,hz_strength,dt)
	right_args = ([2,3],j_strength,hx_strength,hz_strength,dt)
	return  elem_mult_matrices(bs=bs,bps=bps,arguments=(get_localham_elem,left_args,get_localham_elem,right_args))
end

function get_4ham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	left_args = (j_strength,hx_strength,hz_strength,dt)
	right_args = ([3,4],j_strength,hx_strength,hz_strength,dt)
	return  elem_mult_matrices(bs=bs,bps=bps,arguments=(get_3ham_elem,left_args,get_localham_elem,right_args))
end

function get_5ham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	left_args = (j_strength,hx_strength,hz_strength,dt)
	right_args = ([4,5],j_strength,hx_strength,hz_strength,dt)
	return  elem_mult_matrices(bs=bs,bps=bps,arguments=(get_4ham_elem,left_args,get_localham_elem,right_args))
end

function get_6ham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	left_args = (j_strength,hx_strength,hz_strength,dt)
	right_args = ([5,6],j_strength,hx_strength,hz_strength,dt)
	return  elem_mult_matrices(bs=bs,bps=bps,arguments=(get_5ham_elem,left_args,get_localham_elem,right_args))
end

function get_7ham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	left_args = (j_strength,hx_strength,hz_strength,dt)
	right_args = ([6,7],j_strength,hx_strength,hz_strength,dt)
	return  elem_mult_matrices(bs=bs,bps=bps,arguments=(get_6ham_elem,left_args,get_localham_elem,right_args))
end

function get_8ham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	left_args = (j_strength,hx_strength,hz_strength,dt)
	right_args = ([7,8],j_strength,hx_strength,hz_strength,dt)
	return  elem_mult_matrices(bs=bs,bps=bps,arguments=(get_7ham_elem,left_args,get_localham_elem,right_args))
end

function get_9ham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	left_args = (j_strength,hx_strength,hz_strength,dt)
	right_args = ([8,9],j_strength,hx_strength,hz_strength,dt)
	return  elem_mult_matrices(bs=bs,bps=bps,arguments=(get_8ham_elem,left_args,get_localham_elem,right_args))
end

function get_10ham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	left_args = (j_strength,hx_strength,hz_strength,dt)
	right_args = ([9,10],j_strength,hx_strength,hz_strength,dt)
	return  elem_mult_matrices(bs=bs,bps=bps,arguments=(get_9ham_elem,left_args,get_localham_elem,right_args))
end

function get_count_ham_elem(; kwargs...)
	bs = get(kwargs, :bs, 1)
	bps = get(kwargs, :bps, 2)
	site_count,j_strength,hx_strength,hz_strength,dt = get(kwargs, :arguments, 3)
	local_args = j_strength,hx_strength,hz_strength,dt
	if site_count == 3
		return get_3ham_elem(bs=bs,bps=bps,arguments=local_args)
	elseif site_count == 4
		return get_4ham_elem(bs=bs,bps=bps,arguments=local_args)
	elseif site_count == 5
		return get_5ham_elem(bs=bs,bps=bps,arguments=local_args)
	elseif site_count == 6
		return get_6ham_elem(bs=bs,bps=bps,arguments=local_args)
	elseif site_count == 7
		return get_7ham_elem(bs=bs,bps=bps,arguments=local_args)
	elseif site_count == 8
		return get_8ham_elem(bs=bs,bps=bps,arguments=local_args)
	elseif site_count == 9
		return get_9ham_elem(bs=bs,bps=bps,arguments=local_args)
	elseif site_count == 10
		return get_10ham_elem(bs=bs,bps=bps,arguments=local_args)
	else
		println("Don't have func for $site_count Sites")
	end
	return
end

function get_expectation(operator_func,operator_args,wavefunc_func)

end

#=
function get_inter_elem(bs,bps,which_qubits,strengths,dt)
	site_count = length(bs)
	final_val = im*0.0
	j_strength,hz_strength = strengths
	for i in 1:2^site_count
		summed_index = int_to_binary(i-1,site_count)
		
		xx_part = get_exp_two_qubit_elem(xx,bs,summed_index,which_qubits,j_strength,dt)
		zz_part = get_exp_single_qubit_elem(z,summed_index,bps,which_qubits[1],hz_strength,dt) + get_exp_single_qubit_elem(z,summed_index,bps,which_qubits[2],hz_strength,dt)
		
		final_val += zz_part * xx_part
	end
	return final_val
end


function get_leftx_on_inter_elem(bs,bps,which_qubits,j_strength,hx_strength,hz_strength,dt)
	site_count = length(bs)
	final_val = im*0.0
	for i in 1:2^site_count
		summed_index = int_to_binary(i-1,site_count)
		
		x_part = get_exp_single_qubit_elem(x,bs,summed_index,which_qubits[1],hx_strength,dt)
		inter_part = get_inter_elem(summed_index,bps,which_qubits,j_strength,hz_strength,dt)

		final_val += x_part * inter_part
	end
	return final_val
end

function get_localham_elem(bs,bps,which_qubits,j_strength,hx_strength,hz_strength,dt)
	site_count = length(bs)
	final_val = im*0.0
	for i in 1:2^site_count
		summed_index = int_to_binary(i-1,site_count)
		
		left_part = get_leftx_on_inter_elem(bs,summed_index,which_qubits,j_strength,hx_strength,hz_strength,dt)
		x_part = get_exp_single_qubit_elem(x,summed_index,bs,which_qubits[1],hx_strength,dt)

		final_val += left_part * x_part
	end
	return final_val
end
=#

function make_site(state)
	site_tensor = ITensor(Index(2))
	site_tensor[1] = 0
	if state == 0
		site_tensor[:] = [0,1]
	elseif state == 1
		site_tensor[:] = [1,0]
	end
	return site_tensor
end

function get_wavefunc_givenorg(local_org)
	site1 = make_site(local_org[1])
	seq_wavefunc_localorg = [site1]
	for j in 2:length(local_org)
		next_site = make_site(local_org[j])
		append!(seq_wavefunc_localorg,[next_site])
	end
	wavefunc_localorg = prod(seq_wavefunc_localorg)
	return wavefunc_localorg
end

function turn_matrix_into_tensor(mat)
	a1 = Index(size(mat)[1])
	a2 = Index(size(mat)[2])
	tensor_version = ITensor(eltype(mat),mat,a1,a2)
	return tensor_version
end

function normalize_wavefunc(wavefunc)
	norm_factor = real((conj(wavefunc) * wavefunc)[1])
	normed_wavefunc = (1/sqrt(norm_factor)) .* wavefunc
	return normed_wavefunc
end

function make_rand_wavefunc(site_count)
	base_wavefunc = randomITensor([Index(2) for i in 1:site_count])
	normed_wavefunc = normalize_wavefunc(base_wavefunc)
	return normed_wavefunc
end

#=
function find_nonzero_exp_elems(mat,site_count,which_qubit,strength,dt)
	data = []
	for i in 1:2^site_count
		for j in 1:2^site_count
			bs = int_to_binary(i-1,site_count)
			bps = int_to_binary(j-1,site_count)
			if length(which_qubit) > 1
				local_value = get_exp_two_qubit_elem(bs,bps,mat,which_qubit,strength,dt)
			else
				local_value = get_exp_single_qubit_elem(bs,bps,mat,which_qubit,strength,dt)
			end
			if local_value != 0.0
				append!(data,[[i,j,local_value]])
			end
		end
	end
	return data
end

function find_any_nonzero_elems(mat)
	data = []
	for i in 1:size(mat)[1]
		for j in 1:size(mat)[2]
			if mat[i,j] != 0.0
				append!(data,[[i,j,mat[i,j]]])
			end
		end
	end
	return data
end

function make_mat_from_nonzeros(data,dims)
	mat = im.*zeros(dims[1],dims[2])
	for i in 1:length(data)
		mat[Int(real(data[i][1])),Int(real(data[i][2]))] = data[i][3]
	end
	return mat
end

function multiply_matrices(A,B)
	vals = Dict{Tuple{Int,Int},ComplexF64}()
    	for i = 1:length(A)
        	rowA, colA, valA = A[i]
		for j = 1:length(B)
		    	rowB, colB, valB = B[j]
		    	if colA == rowB
		        	key = (Int(rowA), Int(colB))
		        	if haskey(vals, key)
		            		vals[key] += valA * valB
		        	else
		           		vals[key] = valA * valB
		        	end
		    	end
		end
	end
    
    # Convert the dictionary of non-zero values to an array of tuples
    non_zero_C = [[key[1], key[2], vals[key]] for key in keys(vals)]
    
    return non_zero_C
end

function multiply_matrices2(mat1,mat2,site_count)
	locs_1,vals_1 = mat1
	locs_2,vals_2 = mat2
	locs_rez = []
	vals_rez = []
	all_possibilities = [k for k in 1:2^site_count]
	for i in 1:2^site_count
		for j in 1:2^site_count
			inner_sum_indices = findall(v1 -> ([i,v1] in locs_1) && ([v1,j] in locs_2),all_possibilities)
			if length(inner_sum_indices) != 0
				found_vals_1 = [vals_1[findfirst(x->x==[i,inner_sum_indices[b]],locs_1)] for b in 1:length(inner_sum_indices)]
				found_vals_2 = [vals_2[findfirst(x->x==[inner_sum_indices[b],j],locs_2)] for b in 1:length(inner_sum_indices)]
				result = sum(found_vals_1 .* found_vals_2)
				append!(locs_rez,[[i,j]])
				append!(vals_rez,[result])
			end
		end
	end
	return locs_rez,vals_rez
end

function get_locvals_local_ham(which_qubits,site_count,j_strength,hz_strength,hx_strength,dt)
	x_part = find_nonzero_exp_elems(x,site_count,which_qubits[1],hx_strength,dt)
	xx_part = find_nonzero_exp_elems(xx,site_count,which_qubits,j_strength,dt)
	z_part1 = find_nonzero_exp_elems(z,site_count,which_qubits[1],hz_strength,dt)
	z_part2 = find_nonzero_exp_elems(z,site_count,which_qubits[2],hz_strength,dt)
	
	comb_z_ons = multiply_matrices(z_part1,z_part2)
	comb_zpart = multiply_matrices(xx_part,comb_z_ons)
	comb_xleft = multiply_matrices(x_part,comb_zpart)
	final_form = multiply_matrices(comb_xleft,x_part)
	return final_form
end

function get_lv_fullham(site_count,j_strength,hz_strength,hx_strength,dt)
	seq_ham = [get_locvals_local_ham([1,2],site_count,j_strength,hz_strength,hx_strength,dt)]
	for i in 2:site_count
		println(i)
		next_sites = [i,i+1]
		if i == site_count
			next_sites[2] = 1
		end
		next_local_ham = get_locvals_local_ham(next_sites,site_count,j_strength,hz_strength,hx_strength,dt)
		if i == 2
			append!(seq_ham,[next_local_ham])
		else
			seq_ham[2] = next_local_ham
		end
		seq_ham[1] = multiply_matrices(seq_ham[1],seq_ham[2])
	end
	return seq_ham[1]
end
=#
function make_manybody_form(mat,site_count,which_qubit)
	full_mat = im.*zeros(2^site_count,2^site_count)
	for i in 1:2^site_count
		for j in 1:2^site_count
			bs = int_to_binary(i-1,site_count)
			bps = int_to_binary(j-1,site_count)
			if length(which_qubit) > 1
				full_mat[i,j] = get_two_qubit_elem(bs=bs,bps=bps,arguments=(mat,which_qubit))
			else
				full_mat[i,j] = get_single_qubit_elem(bs=bs,bps=bps,arguments=(mat,which_qubit))
			end
		end
	end
	return full_mat
end

function get_exp_xpart(which_site,site_count,dt,hx_strength,order=2)
	coeff = -im * dt * hx_strength
	if order == 2
		coeff *= 0.5
	end
	fin_x = exp(coeff .* make_manybody_form(x,site_count,which_site))
	return fin_x
end

function get_exp_zpart(which_sites,site_count,j_strength,hz_strength,dt)
	coeff_int = -im * j_strength * dt
	coeff_ons = -im * 0.5 * dt * hz_strength
	int_part = exp(coeff_int .* make_manybody_form(xx,site_count,which_sites))
	ons_part = exp(coeff_ons .* (make_manybody_form(z,site_count,which_sites[1]) + make_manybody_form(z,site_count,which_sites[2])))
	fin_z = int_part * ons_part
	return fin_z,int_part,ons_part
end

function get_mbham_local(site_count,which_sites,j_strength,hz_strength,hx_strength,dt,order=2)
	zpart = get_exp_zpart(which_sites,site_count,j_strength,hz_strength,dt)[1]
	xpart = get_exp_xpart(which_sites[1],site_count,dt,hx_strength,order)
	if order == 2
		fin_ham = xpart * zpart * xpart
	else
		fin_ham = zpart * xpart
	end
	return fin_ham
end

function get_full_ham(site_count,j_strength,hz_strength,hx_strength,dt,order=2)
	seq_ham = [get_mbham_local(site_count,[1,2],j_strength,hz_strength,hx_strength,dt,order),Matrix{ComplexF64}(undef,2^site_count,2^site_count)]
	for i in 2:site_count-1
		next_sites = [i,i+1]
		#if i == site_count
		#	next_sites[2] = 1
		#end
		next_contrib = get_mbham_local(site_count,next_sites,j_strength,hz_strength,hx_strength,dt,order)
		seq_ham[2] = next_contrib
		seq_ham[1] = prod(seq_ham)
	end
	return seq_ham[1]
end

function make_tens_wavefunc_vec(input_wavefunc)
	comb = combiner(inds(input_wavefunc))
	flat_tens = comb * input_wavefunc
	return Vector(flat_tens)
end

function get_local_magnetization(which_site,wavefunc)
	site_count = Int(log(2,length(wavefunc)))
	full_z = make_manybody_form(z,site_count,which_site)
	magn_val = (transpose(wavefunc) * full_z * wavefunc) / (transpose(wavefunc) * wavefunc)
	return magn_val
end

function get_arrow_shifts(loc_bit)
	yp = cos(angle(loc_bit))/10
	xp = sin(angle(loc_bit))/10
	return xp,yp
end

function get_average_magnetization(all_magns,time_steps)
	return [mean(all_magns[i][1:end]) for i in 1:time_steps+1]
end

function get_twosite_correlation(sites,wavefunc)
	site_count = Int(log(2,length(wavefunc)))
	first_elem = (transpose(wavefunc) * make_manybody_form(zz,site_count,sites) * wavefunc) / (transpose(wavefunc) * wavefunc)
	second_elem = get_local_magnetization(sites[1],wavefunc) * get_local_magnetization(sites[2],wavefunc)
	return first_elem - second_elem
end

function get_avg_correl_wdists(wavefunc)
	sites_count = Int(log(2,length(wavefunc)))
	all_corrs = [0.0*im for i in 1:sites_count]
	for i in 1:sites_count
		local_val = get_twosite_correlation([Int(ceil(sites_count/2)),i],wavefunc)
		all_corrs[i] = local_val
	end
	return all_corrs
end

function plot_site_mag_time_ev(site_vals,all_magns,all_times,site_count,time_steps)
	for i in 1:n
		for j in 1:time_steps+1
			xcomp,ycomp = get_arrow_shifts(magns[j][i])
			plot([site_vals[i],site_vals[i] + xcomp],[all_times[j],all_times[j] + ycomp],"-b")
		end
		plot([site_vals[i],site_vals[i]],[all_times[1],all_times[end]],"-k")
	end
	xlabel("Sites")
	ylabel("Time (x10)")
	return
end

function do_trotter_step(input_wavefunc,hamilt)
	return hamilt * input_wavefunc
end

#=
count = 3
org = [0 for i in 1:count]
n = length(org)
first_wavefunc = sparse(make_tens_wavefunc_vec(get_wavefunc_givenorg(org)))

steps = [5,10,25,50]
final_time = 0.05
hx = 0.0
js = 2.0
hz = 0.25
for i in 1:length(steps)
	time_steps = steps[i]
	dt = final_time/time_steps
	ham = get_full_ham(n,js,hz,hx,dt)
	
	next_wavefunc = do_trotter_step(first_wavefunc,ham)
	for i in 2:time_steps
		next_wavefunc = do_trotter_step(next_wavefunc,ham)
	end
	
end
=#












"fin"

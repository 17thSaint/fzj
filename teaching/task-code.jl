using Einsum,ITensors

x = [0 1;1 0]
xx = [0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0]
z = [1 0;0 -1]
zz = [1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1]
id = [1 0;0 1]

function get_flattened_index(b_list)
	return sum(b_list .* [2^(length(b_list) - i) for i in 1:length(b_list)]) + 1
end

function get_single_qubit_elem(mat,bs,bps,which_qubit)
	prod_parts = []
	for j in 1:length(bs)
		if j != which_qubit
			append!(prod_parts,[bs[j] == bps[j]])
		end
	end
	return prod(prod_parts) * mat[bs[which_qubit]+1,bps[which_qubit]+1] 
end

function get_two_qubit_elem(mat,bs,bps,which_qubits)
	prod_parts = []
	for j in 1:length(bs)
		if j != which_qubits[1] && j != which_qubits[2]
			append!(prod_parts,[bs[j] == bps[j]])
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

n = 2
downup_wavefunc = get_wavefunc_givenorg([0,1])
x_tens = turn_matrix_into_tensor(x)
phi_form = randomITensor([Index(2) for i in 1:n])

@einsum phi_form[b1,b2] := x_tens[b1,b1p] * downup_wavefunc[b1p,b2]











"fin"

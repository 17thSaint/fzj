using LinearAlgebra

x = [0 1;1 0]
xx = [0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0]
z = [1 0;0 -1]
id = [1 0;0 1]

n = 3
qubit = 2
full_mat = zeros(2^n,2^n)

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

#
for b1 in 0:1
	for b1p in 0:1
		for b2 in 0:1
			for b2p in 0:1
				for b3 in 0:1
					for b3p in 0:1
						bs_local = [b1,b2,b3]
						b_elem = get_flattened_index(bs_local)
						bps_local = [b1p,b2p,b3p]
						bp_elem = get_flattened_index(bps_local)
						full_mat[b_elem,bp_elem] = get_two_qubit_elem(xx,bs_local,bps_local,[1,2])
					end
				end
			end
		end
	end
end
#

display(full_mat)








"fin"

using LinearAlgebra

x = [0 1;1 0]
z = [1 0;0 -1]
id = [1 0;0 1]

n = 2
qubit = 2
full_mat = zeros(2^n,2^n)
full_mat2 = zeros(2^n,2^n)

function get_full_mat_elem(mat,bs,bps,which_qubit)
	prod_parts = []
	for j in 1:length(bs)
		if j != which_qubit
			append!(prod_parts,[bs[j] == bps[j]])
		end
	end
	return prod(prod_parts) * mat[bs[which_qubit]+1,bps[which_qubit]+1] 
end

for b1 in 0:1
	for b1p in 0:1
		for b2 in 0:1
			for b2p in 0:1
				bs_local = [b1,b2]
				bps_local = [b1p,b2p]
				full_mat[b1*(2^1) + b2*(2^0) + 1,b1p*(2^1) + b2p*(2^0) + 1] = get_full_mat_elem(x,bs_local,bps_local,qubit)
			end
		end
	end
end

display(full_mat)








"fin"

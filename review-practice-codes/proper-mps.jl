using ITensors

function make_site(state,bond_dim=2)
	local_index = Index(bond_dim)
	site_tensor = ITensor(local_index)
	site_tensor[1] = 0
	if state == 0
		site_tensor[:] = [0 1]
	elseif state == 1
		site_tensor[:] = [1 0]
	end
	return site_tensor,local_index	
end

function next_step_organizations(all_orgs)
	new_all_orgs = vcat(all_orgs,all_orgs) .* 1
	current_count = Int(length(new_all_orgs)/2)
	for i in 1:current_count
		append!(new_all_orgs[i],[0])
		append!(new_all_orgs[i+current_count],[1])
	end
	return new_all_orgs
end

function get_organizations(site_count)
	site_current = [0]
	seq_orgs = [[],[]]
	while length(site_current) < site_count
		if length(site_current) == 1
			seq_orgs[1] = [[0],[1]]
		end
		seq_orgs[2] = next_step_organizations(seq_orgs[1]) .* 1
		append!(site_current,[0])
		seq_orgs[1] = seq_orgs[2] .* 1
	end
	return seq_orgs[2]
end

function get_wavefunc_givenorg(local_org,amplitude=1.0)
	site1, s1_index = make_site(local_org[1])
	seq_wavefunc_localorg = [site1]
	all_indices = [s1_index]
	for j in 2:length(local_org)
		next_site, next_index = make_site(local_org[j])
		append!(all_indices,[next_index])
		if j == 2
			append!(seq_wavefunc_localorg,[next_site])
		else
			seq_wavefunc_localorg[2] = next_site
		end
		seq_wavefunc_localorg[1] = prod(seq_wavefunc_localorg)
	end
	wavefunc_localorg = seq_wavefunc_localorg[1] .* amplitude
	return wavefunc_localorg,all_indices
end

function normalize_wavefunc(wavefunc)
	norm_factor = real((conj(wavefunc) * wavefunc)[1])
	normed_wavefunc = (1/sqrt(norm_factor)) .* wavefunc
	return normed_wavefunc,norm_factor
end

function get_full_wavefunc(amp_tensor,site_count,normed=true)
	all_orgs = get_organizations(site_count)
	first_wavefunc,first_indices = get_wavefunc_givenorg(all_orgs[1],amp_tensor[1])
	all_wavefunc_parts = [first_wavefunc]
	for i in 2:length(all_orgs)
		local_org = all_orgs[i]
		local_wavefunc,local_indices = get_wavefunc_givenorg(local_org,amp_tensor[i])
		replaceinds!(local_wavefunc,local_indices,first_indices)
		append!(all_wavefunc_parts,[local_wavefunc])
	end
	wavefunc = sum(all_wavefunc_parts)
	if normed
		final_wavefunc, coeff = normalize_wavefunc(wavefunc)
		return final_wavefunc,first_indices,coeff
	else
		return wavefunc,first_indices
	end
end

function get_creat_annih_mats(num_parts=1)
	annih = zeros(Float64,(num_parts+1,num_parts+1))
	for i in 1:num_parts
		annih[i,i+1] = sqrt(i)
	end
	creat = transpose(annih)
	return creat, annih
end

function turn_matrix_into_tensor(mat)
	a1 = Index(size(mat)[1])
	a2 = Index(size(mat)[2])
	tensor_version = ITensor(eltype(mat),mat,a1,a2)
	return tensor_version,a1,a2
end

function make_identity_tensor(bond_dim=2)
	a1 = Index(bond_dim)
	a2 = Index(bond_dim)
	iden_tensor = ITensor(a1,a2)
	iden_tensor[1] = 1.0
	for i in 1:bond_dim
		iden_tensor[i,i] = 1.0
	end
	return iden_tensor,a1,a2
end

function assign_new_indices(tensor,indices)
	new_indices = []
	for i in 1:length(indices)
		append!(new_indices,[Index(dim(indices[i]))])
	end
	new_tensor = replaceinds(tensor,indices,new_indices)
	return new_tensor,new_indices
end

function make_onsite_ham(site_count,onsite_strength=1.0)
	iden,d5,d6 = make_identity_tensor()
	creat_mat,annih_mat = get_creat_annih_mats()
	creat_ten,d1,d2 = turn_matrix_into_tensor(creat_mat)
	annih_ten,d3,d4 = turn_matrix_into_tensor(annih_mat)
	ham_indices = []
	seq_ham = Array{ITensor}(undef,2)
	for i in 1:site_count
		local_ham_parts = []
		local_part_indices = []
		if i == 1  
			first_local_annih = replaceind(annih_ten,d3,d2)
			append!(local_ham_parts,[onsite_strength.*(creat_ten*first_local_annih)])
			append!(local_part_indices,[d1,d4])
		else
			append!(local_ham_parts,[iden])
			append!(local_part_indices,[d5,d6])
		end
		for j in 2:site_count
			if j == i
				local_annih = replaceind(annih_ten,d3,d2)
				resulting_comp,local_indices = assign_new_indices(onsite_strength.*(creat_ten*local_annih),[d1,d4])
				append!(local_ham_parts,[resulting_comp])
				append!(local_part_indices,local_indices)
			else
				resulting_comp,local_indices = assign_new_indices(iden,[d5,d6])
				append!(local_ham_parts,[resulting_comp])
				append!(local_part_indices,local_indices)
			end
		end
		local_ham_contrib = prod(local_ham_parts)
		# add svd dim reduction step here? or at each tensor step therefore diff line above
		if i == 1
			append!(ham_indices,local_part_indices)
			seq_ham[1] = local_ham_contrib
		else
			for k in 1:length(ham_indices)
				replaceind!(local_ham_contrib,local_part_indices[k],ham_indices[k])
			end
			seq_ham[2] = local_ham_contrib
			seq_ham[1] = sum(seq_ham)
		end
	end
	final_ham = seq_ham[1]
	return final_ham,ham_indices
end



bond_dim = 2
num_states = 2
num_sites = 2
d1 = Index(1)
d2 = Index(bond_dim)
n_L = Index(num_states)
n_R = Index(num_states)

left_tensor = ITensor(d1,d2,n_L)#randomITensor(d1,d2,n_L) + im .* randomITensor(d1,d2,n_L)
left_tensor[1,1,1] = 0
right_tensor = ITensor(d2,d1,n_R)#randomITensor(d2,d1,n_R) + im .* randomITensor(d2,d1,n_R)
right_tensor[1,1,1] = 0
left_tensor[:,:,1] = [0; 0]
left_tensor[:,:,2] = [0; 1]
right_tensor[:,:,2] = [0; 1]
right_tensor[:,:,1] = [0; 0]
rez_amp_tensor = left_tensor * right_tensor

#rez = get_full_wavefunc(rez_amp_tensor,num_sites)
#show rez[1]

site_count = 4#num_sites
onsite_str = 1.0












"fin"

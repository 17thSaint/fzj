using ITensors

num_states = 2 # b/c spin-1/2 particles
bond_dim = 2
num_sites = 2

function get_creat_annih_mats(num_parts=1)
	annih = zeros(Float64,(num_parts+1,num_parts+1))
	for i in 1:num_parts
		annih[i,i+1] = sqrt(i)
	end
	creat = transpose(annih)
	return creat, annih
end

function turn_matrix_into_tensor(mat,index_1="row",index_2="col")
	a1 = Index(size(mat)[1],index_1)
	a2 = Index(size(mat)[2],index_2)
	tensor_version = ITensor(eltype(mat),mat,a1,a2)
	return tensor_version,a1,a2
end

function get_up_down_state(state)
	if state == 0
		return [0,1]
	elseif state == 1
		return [1,0]
	end
end

function get_organization_state(organization)
	d1 = Index(2,"d1")
	sites = Index(length(organization),"sites")
	organization_state = ITensor(d1,sites)
	for i in 1:length(organization)
		for j in 1:2
			organization_state[d1=>j,sites=>i] = get_up_down_state(organization[i])[j]
		end
	end
	return organization_state,d1,sites
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

function normalize_wavefunc(wavefunc)
	norm_factor = real((conj(wavefunc) * wavefunc)[1])
	normed_wavefunc = (1/sqrt(norm_factor)) .* wavefunc
	return normed_wavefunc,norm_factor
end

function get_wavefunc(amp_tensor,normed=true)
	all_orgs = get_organizations(length(size(amp_tensor)))
	first_org, d1, sites = get_organization_state(all_orgs[1])
	all_org_states = [first_org]
	all_amps = [amp_tensor[1]]
	for i in 2:prod(size(amp_tensor))
		local_org_state, d3, d4 = get_organization_state(all_orgs[i])
		replaceind!(local_org_state,d3,d1)
		replaceind!(local_org_state,d4,sites)
		append!(all_org_states,[local_org_state])
		append!(all_amps,[amp_tensor[i]])
	end
	wavefunc = sum(all_org_states .* all_amps)
	if normed
		final_wavefunc,coeff = normalize_wavefunc(wavefunc)
		return final_wavefunc,d1,sites,all_org_states,all_amps,coeff
	else
		return wavefunc,d1,sites,all_org_states,all_amps
	end
end

function make_hamilt(site_count,onsite_strength=1.0,inter_strength=1.0)
	cre_mat,ani_mat = get_creat_annih_mats()
	onsite_index = Index(site_count)
	intersite_index = Index(site_count)
	hd1 = Index(bond_dim)
	hd2 = Index(bond_dim)
	ham = ITensor(hd1,hd2,onsite_index,intersite_index)
	ham[1,1,1,1] = 0.0
	for i in 1:site_count
		for j in 1:site_count
			if abs(i-j) == 1
				ham[:,:,i,j] = inter_strength.* (cre_mat * ani_mat)
			elseif abs(i-j) == 0
				ham[:,:,i,j] = onsite_strength.* (cre_mat * ani_mat)
			end
		end
	end
	return ham,hd1,hd2,onsite_index,intersite_index
end

function get_ham_expect(amp_tensor,site_count)
	rez_wavefunc,inner_dim,site_dim,orgs,amps,normcoef = get_wavefunc(amp_tensor)
	hamilt,local_dim1,local_dim2,onsite,intersite = make_hamilt(site_count)
	
	replaceind!(hamilt,local_dim2,inner_dim)
	replaceind!(hamilt,intersite,site_dim)
	
	operate_ham_on_wavefunc = hamilt * rez_wavefunc
	
	dag_wavefunc_for_ham = conj(rez_wavefunc) .* 1.0
	
	replaceind!(dag_wavefunc_for_ham,inner_dim,local_dim1)
	replaceind!(dag_wavefunc_for_ham,site_dim,onsite)
	
	expect_val = (dag_wavefunc_for_ham * operate_ham_on_wavefunc)[1]
	return expect_val
end

d1 = Index(1)
d2 = Index(bond_dim)
n_L = Index(num_states)
n_R = Index(num_states)

left_tensor = randomITensor(d1,d2,n_L) + im .* randomITensor(d1,d2,n_L)
right_tensor = randomITensor(d2,d1,n_R) + im .* randomITensor(d2,d1,n_R)
rez_amp_tensor = left_tensor * right_tensor











"fin"

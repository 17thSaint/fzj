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

function make_random_site(dim=2)
	local_index = Index(dim)
	site_tensor = rand([-1,1]) .* randomITensor(local_index) + (im*rand([-1,1])) .* randomITensor(local_index)
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
	creat = zeros(Float64,(num_parts+1,num_parts+1))
	for i in 1:num_parts
		creat[i,i+1] = sqrt(i)
	end
	annih = transpose(creat)
	return creat, annih
end

function turn_matrix_into_tensor(mat)
	a1 = Index(size(mat)[1])
	a2 = Index(size(mat)[2])
	tensor_version = ITensor(eltype(mat),mat,a1,a2)
	return tensor_version,[a1,a2]
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
		# add svd dim reduction step here? or at each tensor step thus diff line above
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

function make_nearest_neighbor_ham(site_count,periodic=false,neighbor_strength=1.0)
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
			append!(local_ham_parts,[annih_ten])
			append!(local_ham_parts,[neighbor_strength.*creat_ten]) # annihilate at size i, create at i+1
			append!(local_part_indices,[d3,d4,d1,d2])
		elseif i == site_count
			if periodic
				new_creat,new_creat_indices = assign_new_indices(creat_ten,[d1,d2])
				append!(local_ham_parts,[new_creat])
				append!(local_part_indices,new_creat_indices)
			else
				continue
			end
		else
			append!(local_ham_parts,[iden])
			append!(local_part_indices,[d5,d6])
		end
		for j in 2:site_count
			if j == i
				if j == site_count
					if periodic
						new_annih,new_annih_indices = assign_new_indices(annih_ten,[d3,d4])
						append!(local_ham_parts,[new_annih])
						append!(local_part_indices,new_annih_indices)
					else
						resulting_comp,local_indices = assign_new_indices(iden,[d5,d6])
						append!(local_ham_parts,[resulting_comp])
						append!(local_part_indices,local_indices)
					end
				else
					new_creat,new_creat_indices = assign_new_indices(creat_ten,[d1,d2])
					new_annih,new_annih_indices = assign_new_indices(annih_ten,[d3,d4])
					append!(local_ham_parts,[new_annih])   # annihilate at size i, create at i+1
					append!(local_ham_parts,[neighbor_strength.*new_creat])
					append!(local_part_indices,new_annih_indices)
					append!(local_part_indices,new_creat_indices)
				end
			elseif j == i+1
				continue
			else
				resulting_comp,local_indices = assign_new_indices(iden,[d5,d6])
				append!(local_ham_parts,[resulting_comp])
				append!(local_part_indices,local_indices)
			end
		end
		local_ham_contrib = prod(local_ham_parts)
		# add svd dim reduction step here? or at each tensor step thus diff line above
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

function get_expect_ham_val(hamilt,hamilt_indices,wavefunc,wavefunc_indices,site_count)
	#replaceind!(hamilt,hamilt_indices[3],wavefunc_indices[1]) # replaces indices to operate ham on wavefunc
	#replaceind!(hamilt,hamilt_indices[1],wavefunc_indices[2])
	for i in 1:site_count
		replaceind!(hamilt,hamilt_indices[2*i-1],wavefunc_indices[end-i+1])
	end	

	ham_on_wavefunc = hamilt * wavefunc

	replaceinds!(ham_on_wavefunc,(hamilt_indices[2],hamilt_indices[4]),wavefunc_indices) # get indices ready to operate bra wavefunc
	expect_val = (conj.(wavefunc) * ham_on_wavefunc)[1]
	
	return expect_val
end

function get_random_mps_coeffs(site_count,bond_dim,possible_states=2)
	d1 = Index(1)
	states_indices = [Index(possible_states) for i in 1:site_count]
	inner_indices = [Index(bond_dim) for i in 1:site_count-1]
	left_tensor = rand([-1,1]).*randomITensor(d1,inner_indices[1],states_indices[1]) + (im * rand([-1,1])) .* randomITensor(d1,inner_indices[1],states_indices[1])
	right_tensor = rand([-1,1]).*randomITensor(inner_indices[end],d1,states_indices[end]) + (im * rand([-1,1])) .* randomITensor(inner_indices[end],d1,states_indices[end])
	all_tensors = [left_tensor]
	for i in 1:site_count-2
		local_tensor = rand([-1,1]).*randomITensor(inner_indices[i],inner_indices[i+1],states_indices[i+1]) + (im * rand([-1,1])) .* randomITensor(inner_indices[i],inner_indices[i+1],states_indices[i+1])
		append!(all_tensors,[local_tensor])
	end
	append!(all_tensors,[right_tensor])
	return prod(all_tensors),all_tensors,d1,states_indices,inner_indices
end

function make_random_wavefunc(site_count,possible_states=2,normed=true)
	site1, s1_index = make_random_site(possible_states)
	seq_wavefunc = [site1]
	all_indices = [s1_index]
	for j in 2:site_count
		next_site, next_index = make_random_site(possible_states)
		append!(all_indices,[next_index])
		if j == 2
			append!(seq_wavefunc,[next_site])
		else
			seq_wavefunc[2] = next_site
		end
		seq_wavefunc[1] = prod(seq_wavefunc)
	end
	wavefunc = seq_wavefunc[1]
	if normed
		return normalize_wavefunc(wavefunc)[1],all_indices
	else
		return wavefunc,all_indices	
	end
end

function make_reshaped_wavefunc(input_wavefunc,site_count,possible_states)
	reshaped_wavefunc = im .* zeros(possible_states,possible_states^(site_count-1))
	for i in 1:possible_states^(site_count-1)
		for j in 1:possible_states
			reshaped_wavefunc[j,i] = input_wavefunc[possible_states*(i-1) + j][1]
		end
	end
	tens_reshaped_wavefunc, new_indices = turn_matrix_into_tensor(reshaped_wavefunc)
	return tens_reshaped_wavefunc, new_indices
end

# to start input matrix is reshaped wavefunc
function do_element_svd(input_matrix,input_indices,keep_all=false,trash_percent=0.001)
	u,s,vt = svd(input_matrix,input_indices[1])
	function while_loop()
		keep = true
		which_dim = 1
		while keep
			if s[which_dim,which_dim] < trash_percent * s[1,1] && !keep_all
				#println("Throwing Out ",size(s)[1]-which_dim+1,"/",size(s)[1])
				keep = false
			elseif which_dim == size(s)[1]
				#println("Keeping All")
				which_dim += 1
				keep = false
			else
				which_dim += 1
			end
		end
		return which_dim
	end
	which_dim = while_loop()
	new_s_indices = [Index(which_dim-1),Index(which_dim-1)]
	new_s = ITensor(new_s_indices)
	new_s[1] = 0.0*im
	for i in 1:which_dim-1
		new_s[i,i] = s[i,i]
	end
	new_u_indices = [inds(u)[1],Index(which_dim-1)]
	new_u = ITensor(new_u_indices)
	new_u[1] = 0.0*im
	for i in 1:size(u)[1]
		for j in 1:which_dim-1
			new_u[i,j] = u[i,j]
		end
	end
	new_vt_indices = [inds(vt)[1],Index(which_dim-1)]
	new_vt = ITensor(new_vt_indices)
	new_vt[1] = 0.0*im
	for i in 1:size(vt)[1]
		for j in 1:which_dim-1
			new_vt[i,j] = vt[i,j]
		end
	end
	replaceind!(new_s,new_s_indices[1],new_u_indices[2])
	replaceind!(new_s,new_s_indices[2],new_vt_indices[2])
	next_beginning_matrix = new_s * new_vt
	next_beginning_indices = inds(next_beginning_matrix)
	return new_u,new_u_indices,next_beginning_matrix,next_beginning_indices,size(s)[1]-which_dim+1
	#return u,s,vt,new_u,new_s,new_vt,which_dim
end

function make_As(input_wavefunc,site_count,possible_states,keep_all=false,trash_percent=0.001)
	all_as = []
	as_indices = []
	all_cs = []
	cs_indices = []
	throwouts = []
	for i in 1:site_count-1
		#println("Doing Site $i")
		if i == 1
			local_matrix, local_indices = make_reshaped_wavefunc(input_wavefunc,site_count,possible_states)
			append!(all_cs,[local_matrix])
			append!(cs_indices,[local_indices])
		else
			local_matrix = all_cs[end]
			local_indices = cs_indices[end]
		end
		local_a,local_a_indices,next_c,next_c_indices,throwout_count = do_element_svd(local_matrix,local_indices,keep_all,trash_percent)
		append!(throwouts,[throwout_count])
		append!(all_as,[local_a])
		append!(as_indices,[local_a_indices])
		append!(all_cs,[next_c])
		append!(cs_indices,[next_c_indices])
	end
	return all_as,as_indices,all_cs,cs_indices,throwouts
end



#=d1 = Index(1)
d2 = Index(bond_dim)
d3 = Index(bond_dim)
d4 = Index(bond_dim)
n_L = Index(num_states)
n_R = Index(num_states)
n_C = Index(num_states)

left_tensor = rand([-1,1]).*randomITensor(d1,d2,n_L) + (im * rand([-1,1])) .* randomITensor(d1,d2,n_L)#ITensor(d1,d2,n_L)
#left_tensor[1,1,1] = 0
right_tensor = rand([-1,1]).*randomITensor(d3,d1,n_R) + (im * rand([-1,1])) .* randomITensor(d3,d1,n_R)#ITensor(d3,d1,n_R)
#right_tensor[1,1,1] = 0
center_tensor = rand([-1,1]).*randomITensor(d2,d3,n_C) + (im * rand([-1,1])) .* randomITensor(d2,d3,n_C)
center_tensor[1] = 0
left_tensor[:,:,1] = [0; 0]
left_tensor[:,:,2] = [0; 1]
center_tensor[:,:,1] = [0 0; 0 0]
center_tensor[:,:,2] = [0 0; 0 1]
right_tensor[:,:,2] = [0; 1]
right_tensor[:,:,1] = [0; 0]

rez_amp_tensor = left_tensor * center_tensor * right_tensor
=#

bond_dim = 2
#local_site_count = num_sites
#nn_ham,nn_ham_indices = make_nearest_neighbor_ham(local_site_count,true)
#onsite_ham,onsite_ham_indices = make_onsite_ham(local_site_count)
#for i in 1:length(nn_ham_indices)
#	replaceind!(nn_ham,nn_ham_indices[i],onsite_ham_indices[i])
#end
#full_ham = nn_ham + onsite_ham
#rez_amp_tensor = get_random_mps_coeffs(local_site_count,bond_dim)[1]
#chosen_wavefunc,wavefunc_indices = get_full_wavefunc(rez_amp_tensor,local_site_count)
#u1,s1,v1 = svd(chosen_wavefunc,wavefunc_indices)
#nrg = get_expect_ham_val(full_ham,onsite_ham_indices,chosen_wavefunc,wavefunc_indices,local_site_count)
#println(nrg)

for i in 2:10
	num_states = i
	for j in 2:10
		num_sites = j
		rand_wavefunc, wavefunc_indices = make_random_wavefunc(num_sites,num_states)
		components = make_As(rand_wavefunc,num_sites,num_states,false,10^-10)
		println("Sites $j: ",components[end][1],"/",i)
	end
end










"fin"

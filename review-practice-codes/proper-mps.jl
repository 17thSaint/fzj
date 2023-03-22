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
	return site_tensor	
end

function make_random_site(dim=2)
	local_index = Index(dim)
	site_tensor = rand([-1,1]) .* randomITensor(local_index) + (im*rand([-1,1])) .* randomITensor(local_index)
	return site_tensor
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
	site1 = make_site(local_org[1])
	seq_wavefunc_localorg = [site1]
	for j in 2:length(local_org)
		next_site = make_site(local_org[j])
		if j == 2
			append!(seq_wavefunc_localorg,[next_site])
		else
			seq_wavefunc_localorg[2] = next_site
		end
		seq_wavefunc_localorg[1] = prod(seq_wavefunc_localorg)
	end
	wavefunc_localorg = seq_wavefunc_localorg[1] .* amplitude
	return wavefunc_localorg
end

function normalize_wavefunc(wavefunc)
	norm_factor = real((conj(wavefunc) * wavefunc)[1])
	normed_wavefunc = (1/sqrt(norm_factor)) .* wavefunc
	return normed_wavefunc,norm_factor
end

function get_full_wavefunc(amp_tensor,site_count,normed=true)
	all_orgs = get_organizations(site_count)
	first_wavefunc = get_wavefunc_givenorg(all_orgs[1],amp_tensor[1])
	all_wavefunc_parts = [first_wavefunc]
	for i in 2:length(all_orgs)
		local_org = all_orgs[i]
		local_wavefunc = get_wavefunc_givenorg(local_org,amp_tensor[i])
		replaceinds!(local_wavefunc,inds(local_wavefunc),inds(first_wavefunc))
		append!(all_wavefunc_parts,[local_wavefunc])
	end
	wavefunc = sum(all_wavefunc_parts)
	if normed
		final_wavefunc, coeff = normalize_wavefunc(wavefunc)
		return final_wavefunc,coeff
	else
		return wavefunc
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
	return tensor_version
end

function get_creat_annih_tens(num_parts=1)
	creat_mat, annih_mat = get_creat_annih_mats(num_parts)
	creat_tens = turn_matrix_into_tensor(creat_mat)
	annih_tens = turn_matrix_into_tensor(annih_mat)
	return creat_tens,annih_tens
end

function make_identity_tensor(bond_dim=2)
	a1 = Index(bond_dim)
	a2 = Index(bond_dim)
	iden_tensor = ITensor(a1,a2)
	iden_tensor[1] = 1.0
	for i in 1:bond_dim
		iden_tensor[i,i] = 1.0
	end
	return iden_tensor
end

# this in principle can be replaced by use of priming
function assign_new_indices(tensor,indices)
	new_indices = []
	for i in 1:length(indices)
		append!(new_indices,[Index(dim(indices[i]))])
	end
	new_tensor = replaceinds(tensor,indices,new_indices)
	return new_tensor,new_indices
end

function make_onsite_ham(site_count,onsite_strength=1.0)
	iden = make_identity_tensor()
	d5,d6 = inds(iden)
	creat_ten,annih_ten = get_creat_annih_tens()
	d1,d2 = inds(creat_ten)
	d3,d4 = inds(annih_ten)
	seq_ham = Array{ITensor}(undef,2)
	for i in 1:site_count
		local_ham_parts = []
		if i == 1  
			first_local_annih = replaceind(annih_ten,d3,d2)
			append!(local_ham_parts,[onsite_strength.*(creat_ten*first_local_annih)])
		else
			append!(local_ham_parts,[iden])
		end
		for j in 2:site_count
			if j == i
				local_annih = replaceind(annih_ten,d3,d2)
				resulting_comp,local_indices = assign_new_indices(onsite_strength.*(creat_ten*local_annih),[d1,d4])
				append!(local_ham_parts,[resulting_comp])
			else
				resulting_comp,local_indices = assign_new_indices(iden,[d5,d6])
				append!(local_ham_parts,[resulting_comp])
			end
		end
		local_ham_contrib = prod(local_ham_parts)
		# add svd dim reduction step here? or at each tensor step thus diff line above
		if i == 1
			seq_ham[1] = local_ham_contrib
		else
			replaceinds!(local_ham_contrib,inds(local_ham_contrib),inds(seq_ham[1]))
			seq_ham[2] = local_ham_contrib
			seq_ham[1] = sum(seq_ham)
		end
	end
	final_ham = seq_ham[1]
	return final_ham
end

function make_nearest_neighbor_ham(site_count,periodic=false,neighbor_strength=1.0)
	iden = make_identity_tensor()
	d5,d6 = inds(iden)
	creat_mat,annih_mat = get_creat_annih_mats()
	creat_ten = turn_matrix_into_tensor(creat_mat)
	d1,d2 = inds(creat_ten)
	annih_ten = turn_matrix_into_tensor(annih_mat)
	d3,d4 = inds(annih_ten)
	seq_ham = Array{ITensor}(undef,2)
	for i in 1:site_count
		local_ham_parts = []
		if i == 1
			append!(local_ham_parts,[annih_ten])
			append!(local_ham_parts,[neighbor_strength.*creat_ten]) # annihilate at size i, create at i+1
		elseif i == site_count
			if periodic
				new_creat,new_creat_indices = assign_new_indices(creat_ten,[d1,d2])
				append!(local_ham_parts,[new_creat])
			else
				continue
			end
		else
			append!(local_ham_parts,[iden])
		end
		for j in 2:site_count
			if j == i
				if j == site_count
					if periodic
						new_annih,new_annih_indices = assign_new_indices(annih_ten,[d3,d4])
						append!(local_ham_parts,[new_annih])
					else
						resulting_comp,local_indices = assign_new_indices(iden,[d5,d6])
						append!(local_ham_parts,[resulting_comp])
					end
				else
					new_creat,new_creat_indices = assign_new_indices(creat_ten,[d1,d2])
					new_annih,new_annih_indices = assign_new_indices(annih_ten,[d3,d4])
					append!(local_ham_parts,[new_annih])   # annihilate at size i, create at i+1
					append!(local_ham_parts,[neighbor_strength.*new_creat])
				end
			elseif j == i+1
				continue
			else
				resulting_comp,local_indices = assign_new_indices(iden,[d5,d6])
				append!(local_ham_parts,[resulting_comp])
			end
		end
		local_ham_contrib = prod(local_ham_parts)
		# add svd dim reduction step here? or at each tensor step thus diff line above
		if i == 1
			seq_ham[1] = local_ham_contrib
		else
			replaceinds!(local_ham_contrib,inds(local_ham_contrib),inds(seq_ham[1]))
			seq_ham[2] = local_ham_contrib
			seq_ham[1] = sum(seq_ham)
		end
	end
	final_ham = seq_ham[1]
	return final_ham
end

function get_expect_ham_val(hamilt,wavefunc,site_count)
	#replaceind!(hamilt,hamilt_indices[3],wavefunc_indices[1]) # replaces indices to operate ham on wavefunc
	#replaceind!(hamilt,hamilt_indices[1],wavefunc_indices[2])
	hamilt_indices = inds(hamilt)
	wavefunc_indices = inds(wavefunc)
	for i in 1:site_count
		replaceind!(hamilt,hamilt_indices[2*i-1],wavefunc_indices[end-i+1])
	end	

	ham_on_wavefunc = hamilt * wavefunc

	replaceinds!(ham_on_wavefunc,(hamilt_indices[2],hamilt_indices[4]),wavefunc_indices) # get indices ready to operate bra wavefunc
	expect_val = (conj.(wavefunc) * ham_on_wavefunc)[1]
	
	return expect_val
end

function get_random_mps_coeffs(site_count,possible_states=2)
	d1 = Index(1)
	states_indices = [Index(possible_states) for i in 1:site_count]
	inner_indices = [Index(possible_states) for i in 1:site_count-1]
	left_tensor = rand([-1,1]).*randomITensor(d1,inner_indices[1],states_indices[1]) + (im * rand([-1,1])) .* randomITensor(d1,inner_indices[1],states_indices[1])
	right_tensor = rand([-1,1]).*randomITensor(inner_indices[end],d1,states_indices[end]) + (im * rand([-1,1])) .* randomITensor(inner_indices[end],d1,states_indices[end])
	all_tensors = [left_tensor]
	for i in 1:site_count-2
		local_tensor = rand([-1,1]).*randomITensor(inner_indices[i],inner_indices[i+1],states_indices[i+1]) + (im * rand([-1,1])) .* randomITensor(inner_indices[i],inner_indices[i+1],states_indices[i+1])
		append!(all_tensors,[local_tensor])
	end
	append!(all_tensors,[right_tensor])
	return prod(all_tensors),all_tensors
end

function make_random_wavefunc(site_count,possible_states=2,normed=true)
	site1 = make_random_site(possible_states)
	seq_wavefunc = [site1]
	for j in 2:site_count
		next_site = make_random_site(possible_states)
		if j == 2
			append!(seq_wavefunc,[next_site])
		else
			seq_wavefunc[2] = next_site
		end
		seq_wavefunc[1] = prod(seq_wavefunc)
	end
	wavefunc = seq_wavefunc[1]
	if normed
		return normalize_wavefunc(wavefunc)[1]
	else
		return wavefunc
	end
end

function make_reshaped_wavefunc(input_wavefunc,site_count,possible_states)
	reshaped_wavefunc = im .* zeros(possible_states,possible_states^(site_count-1))
	for i in 1:possible_states^(site_count-1)
		for j in 1:possible_states
			reshaped_wavefunc[j,i] = input_wavefunc[possible_states*(i-1) + j][1]
		end
	end
	tens_reshaped_wavefunc = turn_matrix_into_tensor(reshaped_wavefunc)
	return tens_reshaped_wavefunc
end

function get_keeping_type(how_keeping,keeping_val=0)
	if how_keeping == "all"
		return true
	elseif how_keeping == "count" || how_keeping == "percent"
		return keeping_val
	end
end


function reshape_Cs(current_c,next_bond_dim)
	next_right_dim = Int(size(current_c)[2]/next_bond_dim)
	reshaped = ITensor(Index(next_bond_dim),Index(size(current_c)[1]),Index(next_right_dim))
	reshaped[1] = 0.0*im
	for i in 1:next_bond_dim
		reshaped[i,:,:] = current_c[:,(i-1)*next_right_dim+1:i*next_right_dim]
	end
	c1 = combiner(inds(reshaped)[1:2])
	flattened = c1 * reshaped
	return flattened
end

# to start input matrix is reshaped wavefunc
# keeping type default is throw away <1%
function do_element_svd(input_matrix,keeping_type=0.01)
	u,s,vt = svd(input_matrix,inds(input_matrix)[1])
	max_dim = size(s)[1]
	if typeof(keeping_type) == Bool
		keep_all = true
	else
		keep_all = false
		if typeof(keeping_type) == Float64
			trash_percent = keeping_type
		elseif typeof(keeping_type) == Int64
			max_dim = keeping_type
		end
	end
	function while_loop()
		keep = true
		which_dim = 1
		while keep
			if !keep_all && typeof(keeping_type) == Float64 && s[which_dim,which_dim] < trash_percent * s[1,1]
				println("Throwing Out ",size(s)[1]-which_dim+1,"/",size(s)[1])
				keep = false
			elseif which_dim == max_dim
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
	
	return new_u,next_beginning_matrix,size(s)[1]-which_dim+1
end

function check_A_sym(input_A)
	transp_A = setprime(transpose(input_A),0)
	return transp_A == input_A
end

function check_A_diag(input_A)
	if size(input_A)[1] == size(input_A)[2]
		mat_form = Matrix(input_A,inds(input_A))
		return !(false in isapprox.(transpose(mat_form),mat_form,atol=10^-7))
	else
		println("Not Square")
		return true
	end
end

function make_As(input_wavefunc,site_count,possible_states,keeping_type=0.01,labels=false)
	all_as = []
	all_cs = []
	throwouts = []
	
	for i in 1:site_count-1
		#println("Doing Site $i")
		if i == 1
			local_matrix = make_reshaped_wavefunc(input_wavefunc,site_count,possible_states)
			append!(all_cs,[local_matrix])
		else
			local_matrix = reshape_Cs(all_cs[end],possible_states)
		end

		local_a,next_c,throwout_count = do_element_svd(local_matrix,keeping_type)
		#=
		if labels
			if i == 1
				replaceind!(local_a,inds(local_a)[1],addtags(inds(local_a)[1],"s1"))
				replaceind!(local_a,inds(local_a)[2],addtags(inds(local_a)[2],"a1"))
			else
				ind1 = i - 1
				replaceind!(local_a,inds(local_a)[1],addtags(inds(local_a)[1],"a$ind1"))
				replaceind!(local_a,inds(local_a)[2],addtags(inds(local_a)[2],"a$i"))
			end
		end
		=#
		append!(throwouts,[throwout_count])
		append!(all_as,[local_a])
		append!(all_cs,[next_c])
		#=
		if i != 1
			if !check_A_sym(local_a)
				println("Middle A at site $i NOT SYM")
				return all_as,all_cs
			elseif !check_A_diag(local_a)
				println("Middle A at site $i NOT diagonal")
				return all_as,all_cs
			end
		end
		=#
	end
	if typeof(keeping_type) == Float64
		return all_as,all_cs,throwouts
	else
		return all_as,all_cs
	end
end

function get_eigvec_tensor(state,possible_states)
	eigvec_tensor = ITensor(Index(possible_states))
	eigvec_tensor[state] = 1.0
	return eigvec_tensor
end

# the 1 state is (1 0) and the 2 state is (0 1)
# j is s'
function get_first_rho(a1)
	dim1 = size(a1)[1]
	rho1_coeffs = conj(transpose(a1)) * a1
	mat_local_dens = zeros(dim1,dim1) .+ im*0.0
	for i in 1:dim1
		for j in 1:dim1
			coeff = sum(rho1_coeffs[j,k1,i,k2] for k2 in 1:size(a1)[2] for k1 in 1:size(a1)[2])
			mat_local_dens[i,j] = coeff
		end
	end
	println(tr(mat_local_dens))
	local_dens_mat = turn_matrix_into_tensor(mat_local_dens)
	return local_dens_mat
end

function get_inner_rho(ai)
	coeff = sum((ai * conj(transpose(ai)))[i,j] for i in 1:size(ai)[1] for j in 1:size(ai)[1])
	local_dens_mat = ITensor(Index(size(ai)[1]),Index(size(ai)[1]))
	local_dens_mat = ones(size(ai)[1],size(ai)[1]) .* coeff
	return local_dens_mat
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

#bond_dim = 2
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



num_sites = 5
num_states = 3
keeping = "percent"
keep_count = 0.01
keep_type = get_keeping_type(keeping,keep_count)
rand_wavefunc = make_random_wavefunc(num_sites,num_states)
as,cs = make_As(rand_wavefunc,num_sites,num_states,keep_type)









"fin"

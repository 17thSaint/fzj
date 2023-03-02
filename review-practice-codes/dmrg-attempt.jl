using ITensors

num_states = 2 # b/c spin-1/2 particles
bond_dim = 2
num_sites = 3

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

function get_wavefunc(amp_tensor)
	all_orgs = get_organizations(length(size(f)))
	for i in 2:prod(size(amp_tensor))
		
	end
end


f, d11, d12 = get_organization_state([0,1])
g, d21, d22 = get_organization_state([1,0])
replaceind!(g,d21,d11)
replaceind!(g,d22,d12)
z = f + g


d1 = Index(1)
d2 = Index(bond_dim)
n_L = Index(num_states)
n_R = Index(num_states)

left_tensor = randomITensor(d1,d2,n_L)
right_tensor = randomITensor(d2,d1,n_R)













"fin"

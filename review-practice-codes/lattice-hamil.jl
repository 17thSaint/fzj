using LinearAlgebra,PyPlot

annih = [0 1; 0 0]
creat = [0 0; 1 0]
id = [1 0; 0 1]


function get_position(which,edge_count,spacing)
	ycomp = spacing * (ceil(which/edge_count) - 1)
	xcomp = ((which-1) % edge_count) * spacing
	return xcomp + im*ycomp
end

function get_j(j,k,t,phi,edge_count,spacing)
	z = get_position(k,edge_count,spacing) - get_position(j,edge_count,spacing)
	gz = (-1)^(real(z) + imag(z) + imag(z)*real(z))
	exp_phi_part = (get_position(j,edge_count,spacing) * conj(z) - conj(get_position(j,edge_count,spacing)) * z + abs2(z) ) * pi * 0.5
	exp_other_part = -pi * 05 * abs2(z)
	coef = gz * exp(exp_other_part)
	final_val = t * gz * exp(exp_other_part + phi * exp_phi_part)
	return final_val
end

function which_mat(choice)
	if choice == "id"
		return id
	elseif choice == "c"
		return creat
	elseif choice == "a"
		return annih
	else
		return "Error"
	end
end

function make_component_matrix(j,k,edge_count)
	which_mats = ["id" for i in 1:edge_count^2]
	c_site = j
	an_site = k
	which_mats[c_site] = "c"
	which_mats[an_site] = "a"

	series_of_prods = Vector{Matrix{Int64}}(undef,edge_count^2-1)
	for i1 in 2:edge_count^2
		if i1 == 2
			this_element = which_mat(which_mats[1])
		else
			this_element = series_of_prods[i1-2]
		end
		new_mat = kron(this_element,which_mat(which_mats[i1]))
		series_of_prods[i1-1] = new_mat
	end
	return series_of_prods[end]
end

function make_ham_mat(edge_count,t,phi,spacing)
	ham_mat = zeros(Int64,2^(edge_count^2),2^(edge_count^2))
	for i2 in 1:edge_count^2
		#println("Site $i2")
		for j2 in 1:edge_count^2
			if j2 == i2
				continue
			end
			
			coeff = get_j(j2,i2,t,phi,edge_count,spacing)
			matrix_comp = make_component_matrix(j2,i2,edge_count)
			ham_mat += coeff.*matrix_comp
		end
	end
	
	return eigvals(ham_mat),eigvecs(ham_mat),ham_mat
end

lat_sep = 1.0
t_val = 1.0
edges_count = [2,3]
phis = [1/3]
#phi_val = 1/3
#
percent_thru_eigvals = [[] for i in 1:length(edges_count)]
energies = [[] for i in 1:length(edges_count)]
for i in 1:length(edges_count)
	num_edge_sites = edges_count[i]
	for j in 1:length(phis)
		phi_val = phis[j]
		en_vals,eivecs,mat = make_ham_mat(num_edge_sites,t_val,phi_val,lat_sep)
		display(eivecs)
		energies[i] = en_vals./maximum(en_vals)
		percent_thru_eigvals[i] = [(k-1)/length(en_vals) for k in 1:length(en_vals)]
	end
end
#
for i in 1:length(edges_count)
	edges = edges_count[i]
	plot(percent_thru_eigvals[i],energies[i],"-p",label="$edges")
end
legend()


#=
lat_sep = 1.0
t_val = 1.0
site_counts = [10]
phis_count = 10
phis = []

#eivec_alls = [[] for i in 1:phis_count+1]
energy_vals = [[] for i in 1:length(site_counts)]

for j in 1:phis_count+1
	phi_val = 0.1 + (j-1)*0.9/phis_count
	append!(phis,phi_val)
	for i in 1:length(site_counts)
		num_sites = site_counts[i]
		println("Working on $num_sites Sites")
		en_vals,eivecs,mat = make_ham_mat(num_sites,t_val,phi_val,lat_sep)
		append!(energy_vals[i],minimum(en_vals))
		#display(eivecs)
		#println(en_vals)
	end
end
#
for p in 1:length(site_counts)
	sites = site_counts[p]
	plot(phis,energy_vals[p],label="$sites")
end
legend()
=#










"fin"

using LinearAlgebra,PyPlot

function get_creat_annih_ops(num_parts=1)
	annih = zeros(Float64,(num_parts+1,num_parts+1))
	for i in 1:num_parts
		annih[i,i+1] = sqrt(i)
	end
	creat = transpose(annih)
	return creat, annih
end

function get_position(which,edge_count,spacing)
	ycomp = spacing * (ceil(which/edge_count) - 1)
	xcomp = ((which-1) % edge_count) * spacing
	return xcomp + im*ycomp
end

function get_torus_distance_btw(pos_1,pos_2,edge_count,spacing)	#enforces periodic boundary conditions of square lattice
	x1, y1 = real(pos_1), imag(pos_1)
	x2, y2 = real(pos_2), imag(pos_2)
	
	abs_dists_x = [abs(x1 - x2), spacing*edge_count - abs(x1 - x2)]
	if sort(abs_dists_x) == abs_dists_x
		dist_x = x1 - x2
	else
		dist_x = sign(x1-x2) * (-1) * abs_dists_x[2]
	end
	
	abs_dists_y = [abs(y1 - y2), spacing*edge_count - abs(y1 - y2)]
	if sort(abs_dists_y) == abs_dists_y
		dist_y = y1 - y2
	else
		dist_y = sign(y1-y2) * (-1) * abs_dists_y[2]
	end
	
	z_val = dist_x + im*dist_y
	return z_val
end

function get_j(j,k,t,phi,edge_count,spacing=1.0,periodic=true)
	z = get_position(k,edge_count,spacing) - get_position(j,edge_count,spacing)
	if periodic
		z = get_torus_distance_btw(get_position(k,edge_count,spacing),get_position(j,edge_count,spacing),edge_count,spacing) # periodic boundary conditions
	end
	gz = (-1)^(real(z) + imag(z) + imag(z)*real(z))
	exp_phi_part = (get_position(j,edge_count,spacing) * conj(z) - conj(get_position(j,edge_count,spacing)) * z + abs2(z) ) * pi * 0.5
	exp_other_part = -pi * 05 * abs2(z)
	coef = gz * exp(exp_other_part)
	final_val = t * gz * exp(exp_other_part + phi * exp_phi_part)
	return final_val
end

function which_mat(choice,num_parts=1)
	creat, annih = get_creat_annih_ops(num_parts)
	if choice == "id"
		return Matrix(I,num_parts+1,num_parts+1)
	elseif choice == "c"
		return creat
	elseif choice == "a"
		return annih
	else
		return "Error"
	end
end

function make_component_matrix(j,k,edge_count,num_parts=1)
	which_mats = ["id" for i in 1:edge_count^2]
	c_site = j
	an_site = k
	which_mats[c_site] = "c"
	which_mats[an_site] = "a"

	series_of_prods = Vector{Matrix{Float64}}(undef,edge_count^2-1)
	for i1 in 2:edge_count^2
		if i1 == 2
			this_element = which_mat(which_mats[1],num_parts)
		else
			this_element = series_of_prods[i1-2]
		end
		new_mat = kron(this_element,which_mat(which_mats[i1],num_parts))
		series_of_prods[i1-1] = new_mat
		if i1 > 2
			series_of_prods[i1-2] = Matrix{Float64}(undef,(2,2))
		end
	end
	return series_of_prods[end]
end

function make_ham_mat(edge_count,t,phi,num_parts=1,spacing=1.0,periodic=true)
	ham_mat = zeros(Int64,(num_parts+1)^(edge_count^2),(num_parts+1)^(edge_count^2))
	for i2 in 1:edge_count^2
		#println("Site $i2 / $edge_count")
		for j2 in 1:edge_count^2
			if j2 == i2
				continue
			end
			
			coeff = get_j(j2,i2,t,phi,edge_count,spacing,periodic)
			matrix_comp = make_component_matrix(j2,i2,edge_count,num_parts)
			ham_mat += coeff.*matrix_comp
			#println(i2,", ",j2,", ",coeff)
			#display(ham_mat)
		end
	end
	
	if ham_mat != conj(transpose(ham_mat))
		println("Not Hermitian, Need to Stop")
	end
	
	return eigvals(ham_mat),eigvecs(ham_mat),ham_mat
end

function get_total_number_particles(edge_count,wavefunc)
	count = 0
	for i in 1:edge_count^2
		site_number_operator = make_component_matrix(i,i,edge_count)
		multiplied_wavefunc = real.( (site_number_operator * wavefunc ) ./ wavefunc )
		found_number = false
		display(multiplied_wavefunc)
		display(wavefunc)
		element = 0
		while !found_number
			element += 1
			local_number = multiplied_wavefunc[element] 
			if !isnan(local_number) 
				count += Int(round(local_number))
				println("Found at Element $element: ",Int(round(local_number)),", Count is now $count")
				found_number = true
			else
				println("Not Found at Element $element")
			end
		end
	end
	println("Final Count is $count")
end

function plot_nrgs_range_phi(phis_count,edge_count,num_parts=1,plot_ground_state=true,periodic=true,phi_start=0.1,phi_end=1.0,spacing=1.0)
	phis = [phi_start + (i-1)*(phi_end-phi_start)/phis_count for i in 1:phis_count+1]
	energies = [[] for i in 1:phis_count+1]
	for i in 1:phis_count+1
		phi_val = phis[i]
		en_vals = make_ham_mat(edge_count,1.0,phi_val,num_parts)[1]
		energies[i] = en_vals
	end
	if plot_ground_state
		title_here = "Fermion"
		if num_parts > 1
			title_here = "Boson Np=$num_parts"
		end
		plot(phis,[energies[i][1] for i in 1:phis_count+1],"-p",label="Ne=$edge_count, $title_here")
		xlabel("Phi")
		ylabel("Lowest NRG")
		legend()
	end
	
	return phis, energies
end

count = 30
sites = 3
for particles in [2]
	plot_nrgs_range_phi(count,sites,particles)
end


#=
lat_sep = 1.0
t_val = 1.0
edges_count = [2]
phis_count = 30
phis = [0.1 + (j-1)*(1.0-0.1)/phis_count for j in 1:phis_count]
#phi_val = 1/3
#
#eivecs_here = Matrix{ComplexF64}(undef,(2^(edges_count[1]^2),2^(edges_count[1]^2)))
percent_thru_eigvals = [[[] for j in 1:length(phis)] for i in 1:length(edges_count)]
energies = [[[] for j in 1:length(phis)] for i in 1:length(edges_count)]
for i in 1:length(edges_count)
	num_edge_sites = edges_count[i]
	for j in 1:length(phis)
		#println(100*j/length(phis))
		phi_val = phis[j]
		en_vals, eivecs, mat = make_ham_mat(num_edge_sites,t_val,phi_val,lat_sep,true)
		println(mat == conj(transpose(mat)))
		#for k in 1:1#size(eivecs)[1]
		#	get_total_number_particles(num_edge_sites,eivecs[:,k])
		#end
		#ham_mats = mat
		energies[i][j] = en_vals#./maximum(en_vals)
		percent_thru_eigvals[i][j] = [(k-1)/length(en_vals) for k in 1:length(en_vals)]
	end
end
=#
#=
for i in 1:length(edges_count)
	edges = edges_count[i]
	plot(phis,[energies[i][j][1] for j in 1:length(phis)],"-p",label="$edges")
end
legend()
xlabel("Phi")
ylabel("Lowest Energy")
=#

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

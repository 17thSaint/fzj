using LinearAlgebra,PyPlot

function get_creat_annih_ops(num_parts=1)
	annih = zeros(Float64,(num_parts+1,num_parts+1))
	for i in 1:num_parts
		annih[i,i+1] = sqrt(i)
	end
	creat = transpose(annih)
	return creat, annih
end

function get_position(which,edge_count,spacing=1.0)
	ycomp = spacing * (ceil(which/edge_count) - 1) - 0.5*(edge_count-1)*spacing
	xcomp = (((which-1) % edge_count) * spacing) - 0.5*(edge_count-1)*spacing
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

function make_single_particle_component_matrix(j,k,edge_count)
	j_vec = zeros(Int64,edge_count^2)
	j_vec[j] = 1
	k_vec = zeros(Int64,edge_count^2)
	k_vec[k] = 1
	comp_mat = j_vec * transpose(k_vec)
	return comp_mat
end

function make_single_part_ham(edge_count,phi,periodic=true,spacing=1.0,t=1.0)
	ham_mat = zeros(Float64,(edge_count^2,edge_count^2))
	for j in 1:edge_count^2
		for k in 1:edge_count^2
			if j == k
				continue
			end
			coeff = round(get_j(j,k,t,phi,edge_count,spacing,periodic),digits=10)
			matrix_comp = make_single_particle_component_matrix(j,k,edge_count)
			ham_mat += coeff.*matrix_comp
		end
	end
	
	if ham_mat != conj(transpose(ham_mat))
		println("Not Hermitian, Need to Stop: Edges=$edge_count, Phi=$phi")
	end
	
	return eigvals(ham_mat),eigvecs(ham_mat),ham_mat
end

function make_ham_mat(edge_count,phi,num_parts=1,spacing=1.0,periodic=true,t=1.0)
	ham_mat = zeros(Int64,(num_parts+1)^(edge_count^2),(num_parts+1)^(edge_count^2))
	for i2 in 1:edge_count^2
		#println("Site $i2 / $edge_count")
		for j2 in 1:edge_count^2
			if j2 == i2
				continue
			end
			
			coeff = round(get_j(j2,i2,t,phi,edge_count,spacing,periodic),digits=10)
			matrix_comp = make_component_matrix(j2,i2,edge_count,num_parts)
			ham_mat += coeff.*matrix_comp
			#println(i2,", ",j2,", ",coeff)
			#display(ham_mat)
		end
	end
	
	if ham_mat != conj(transpose(ham_mat))
		println("Not Hermitian, Need to Stop: Edges=$edge_count, Phi=$phi")
	end
	
	return eigvals(ham_mat),eigvecs(ham_mat),ham_mat
end
#=
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
=#

function get_expected_position(wavefunc,edge_count,spacing=1.0) # gets expectation value of position from single particle bases wavefunction
	wavefunc_with_xpos = [wavefunc[i]*real(get_position(i,edge_count,spacing)) for i in 1:length(wavefunc)]
	wavefunc_with_ypos = [wavefunc[i]*imag(get_position(i,edge_count,spacing)) for i in 1:length(wavefunc)]
	final_val_x = abs2(transpose(wavefunc) * wavefunc_with_xpos)
	final_val_y = abs2(transpose(wavefunc) * wavefunc_with_ypos)
	return [final_val_x,final_val_y]
end

function plot_nrgs_range_phi(phis_count,edge_count,single_part=true,num_parts=1,plot_spectrum=false,plot_ground_state=true,periodic=true,phi_start=0.01,phi_end=1.0,spacing=1.0)
	phis = [phi_start + (i-1)*(phi_end-phi_start)/phis_count for i in 1:phis_count+1]
	energies = [[] for i in 1:phis_count+1]
	ground_state_wavefunc = [[] for i in 1:phis_count+1]
	for i in 1:phis_count+1
		#println(round(100*i/(phis_count+1),digits=1),"%")
		phi_val = phis[i]
		if single_part
			en_vals,eigvals,hamil = make_single_part_ham(edge_count,phi_val,periodic,spacing)
		else
			en_vals,eigvals,hamil = make_ham_mat(edge_count,phi_val,num_parts)
		end
		energies[i] = en_vals
		ground_state_wavefunc[i] = eigvals[:,1]
	end
	if plot_ground_state
		plot(phis,[energies[i][1] for i in 1:phis_count+1],"-p",label="Ne=$edge_count")
		xlabel("Phi")
		ylabel("Lowest NRG")
		legend()
	elseif plot_spectrum
		#fig, axs = plt.subplots(1,2)
		#colors = ["b","g","r","c","m","y","k"]
		for j in 1:edge_count^2
			plot(phis,[energies[i][j] for i in 1:phis_count+1],"-b",label="Ne=$edge_count, $j")
			xlabel("Phi")
			ylabel("Lowest NRG")
			legend()
		end
	end
	
	return phis, energies, ground_state_wavefunc
end

function plot_spectra_range_sites(phi_count,start_sites,end_sites)
	colors = ["b","g","r","c","m","y","k"]
	if start_sites == end_sites
		edge_sites = start_sites
		count_color = colors[1]
		phis_here, nrgs_here = plot_nrgs_range_phi(phi_count,edge_sites,true,1,false,false)[1:2]
		for j in 1:edge_sites^2
			if j == 1
				plot(phis_here,[nrgs_here[i][j] for i in 1:phi_count+1],"-$count_color",label="Ne=$edge_sites, $j")
			else
				plot(phis_here,[nrgs_here[i][j] for i in 1:phi_count+1],"-$count_color")
			end
		end
		legend()
		xlabel("Phi")
		title("Energy Spectrum for $edge_sites Edge Sites")
	else
		fig, axs = plt.subplots(1,end_sites-start_sites+1)
		for edge_sites in start_sites:end_sites
			println("Working on $edge_sites Sites")
			count_color = colors[edge_sites-start_sites+1]
			phis_here, nrgs_here = plot_nrgs_range_phi(phi_count,edge_sites,true,1,false,false)[1:2]
			for j in 1:edge_sites^2
				if j == 1
					axs[edge_sites-start_sites+1].plot(phis_here,[nrgs_here[i][j] for i in 1:phi_count+1],"-$count_color",label="Ne=$edge_sites, $j")
				else
					axs[edge_sites-start_sites+1].plot(phis_here,[nrgs_here[i][j] for i in 1:phi_count+1],"-$count_color")
				end
			end
			#axs[edge_sites-start_sites+1].xlabel("Phi")
			#ylabel("Lowest NRG")
			axs[edge_sites-start_sites+1].legend()
		end
		xlabel("Phi")
		title("Energy Spectrum for range Lattice Sites")
	end
end

function bin_nrg_spectra(phi_count,edge_sites,bin_count=10,phi_first=0.01,phi_last=1.0)
	phis_here, nrgs_here = plot_nrgs_range_phi(phi_count,edge_sites,true,1,false,false,true,phi_first,phi_last)[1:2]
	num_nrg_vals = zeros(Int64,phi_count+1)
	for i in 1:phi_count+1
		binned_nrgs = hist(nrgs_here[i],bin_count)
		number_of_nrg_vals_local = 0
		for j in 1:bin_count
			if binned_nrgs[1][j] != 0.0
				number_of_nrg_vals_local += 1
			end
		end
		num_nrg_vals[i] = number_of_nrg_vals_local
	end
	plot(phis_here,num_nrg_vals,"-p",label="$edge_sites")
	xlabel("Phi")
	title("Number of Distinct Energies")
	legend()
	
	return phis_here, num_nrg_vals
end

edges = 10
bins = 50
phis = 50
bin_nrg_spectra(phis,edges,bins)

#=
count = 100
start_site = 3
end_site = 5
for edge_sites in start_site:end_site
	figure()
	phi_vals, nrgs, ground_states = plot_nrgs_range_phi(count,edge_sites,true,1,false,false)
	scatter3D([real(get_position(i,edge_sites))/(edge_sites-1) for i in 1:edge_sites^2],[imag(get_position(i,edge_sites))/(edge_sites-1) for i in 1:edge_sites^2],[0.0 for i in 1:edge_sites^2])
	expected_locations = [get_expected_position(ground_states[i],edge_sites) for i in 1:count+1]
	plot3D([expected_locations[i][1]/(edge_sites-1) for i in 1:count+1],[expected_locations[i][2]/(edge_sites-1) for i in 1:count+1],phi_vals,"-pr")
	zlabel("Phi")
	title("Expected GS Position for Edge Sites=$edge_sites")
end
=#
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

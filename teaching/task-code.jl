using Einsum,ITensors,LinearAlgebra,PyPlot,Statistics

x = [0 1;1 0]
xx = [0 0 0 1;0 0 1 0;0 1 0 0;1 0 0 0]
z = [1 0;0 -1]
zz = [1 0 0 0;0 -1 0 0;0 0 -1 0;0 0 0 1]
id = [1 0;0 1]

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
	if length(prod_parts) < 1
		prod_rez = 1
	else
		prod_rez = prod(prod_parts)
	end
	return prod_rez * mat[bs[which_qubit]+1,bps[which_qubit]+1] 
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

function make_site(state)
	site_tensor = ITensor(Index(2))
	site_tensor[1] = 0
	if state == 0
		site_tensor[:] = [0,1]
	elseif state == 1
		site_tensor[:] = [1,0]
	end
	return site_tensor	
end

function get_wavefunc_givenorg(local_org)
	site1 = make_site(local_org[1])
	seq_wavefunc_localorg = [site1]
	for j in 2:length(local_org)
		next_site = make_site(local_org[j])
		append!(seq_wavefunc_localorg,[next_site])
	end
	wavefunc_localorg = prod(seq_wavefunc_localorg)
	return wavefunc_localorg
end

function turn_matrix_into_tensor(mat)
	a1 = Index(size(mat)[1])
	a2 = Index(size(mat)[2])
	tensor_version = ITensor(eltype(mat),mat,a1,a2)
	return tensor_version
end

function normalize_wavefunc(wavefunc)
	norm_factor = real((conj(wavefunc) * wavefunc)[1])
	normed_wavefunc = (1/sqrt(norm_factor)) .* wavefunc
	return normed_wavefunc
end

function make_rand_wavefunc(site_count)
	base_wavefunc = randomITensor([Index(2) for i in 1:site_count])
	normed_wavefunc = normalize_wavefunc(base_wavefunc)
	return normed_wavefunc
end

function int_to_binary(given::Int,len::Int)
    # Initialize an empty array to store the binary digits
    binary = []
    # Loop until the integer is reduced to zero
    while given > 0
        # Get the remainder of n when divided by 2
        digit = given % 2
        # Prepend the digit to the binary array
        pushfirst!(binary, digit)
        # Integer divide n by 2 to reduce its value
        given = div(given, 2)
    end
    # Return the binary array
    if length(binary) < len
    	return append!([0 for i in 1:(len-length(binary))],binary)
    else
    	return binary
    end
end

function make_manybody_form(mat,site_count,which_qubit)
	full_mat = zeros(2^site_count,2^site_count)
	for i in 1:2^site_count
		for j in 1:2^site_count
			bs = int_to_binary(i-1,site_count)
			bps = int_to_binary(j-1,site_count)
			if length(which_qubit) > 1
				full_mat[i,j] = get_two_qubit_elem(mat,bs,bps,which_qubit)
			else
				full_mat[i,j] = get_single_qubit_elem(mat,bs,bps,which_qubit)
			end
		end
	end
	return full_mat
end

function get_exp_xpart(which_site,site_count,dt,hx_strength,order=2)
	coeff = -im * dt * hx_strength
	if order == 2
		coeff *= 0.5
	end
	fin_x = exp(coeff .* make_manybody_form(x,site_count,which_site))
	return fin_x
end

function get_exp_zpart(which_sites,site_count,j_strength,hz_strength,dt)
	coeff_int = -im * j_strength * dt
	coeff_ons = -im * 0.5 * dt * hz_strength
	int_part = coeff_int .* make_manybody_form(xx,site_count,which_sites)
	ons_part = coeff_ons .* (make_manybody_form(z,site_count,which_sites[1]) + make_manybody_form(z,site_count,which_sites[2]))
	fin_z = exp(int_part + ons_part)
	return fin_z
end

function get_mbham_local(site_count,which_sites,j_strength,hz_strength,hx_strength,dt,order=2)
	zpart = get_exp_zpart(which_sites,site_count,j_strength,hz_strength,dt)
	xpart = get_exp_xpart(which_sites[1],site_count,dt,hx_strength,order)
	if order == 2
		fin_ham = xpart * zpart * xpart
	else
		fin_ham = zpart * xpart
	end
	return fin_ham
end

function get_full_ham(site_count,j_strength,hz_strength,hx_strength,dt,order=2)
	seq_ham = [get_mbham_local(site_count,[1,2],j_strength,hz_strength,hx_strength,dt,order),Matrix{ComplexF64}(undef,2^site_count,2^site_count)]
	for i in 2:site_count
		next_sites = [i,i+1]
		if i == site_count
			next_sites[2] = 1
		end
		next_contrib = get_mbham_local(site_count,next_sites,j_strength,hz_strength,hx_strength,dt,order)
		seq_ham[2] = next_contrib
		seq_ham[1] = prod(seq_ham)
	end
	return seq_ham[1]
end

function make_tens_wavefunc_vec(input_wavefunc)
	comb = combiner(inds(input_wavefunc))
	flat_tens = comb * input_wavefunc
	return Vector(flat_tens)
end

function get_local_magnetization(which_site,wavefunc)
	site_count = Int(log(2,length(wavefunc)))
	full_z = make_manybody_form(z,site_count,which_site)
	magn_val = (transpose(wavefunc) * full_z * wavefunc) / (transpose(wavefunc) * wavefunc)
	return magn_val
end

function get_arrow_shifts(loc_bit)
	yp = cos(angle(loc_bit))/10
	xp = sin(angle(loc_bit))/10
	return xp,yp
end

function get_average_magnetization(all_magns,time_steps)
	return [mean(all_magns[i][1:end]) for i in 1:time_steps+1]
end

function get_twosite_correlation(sites,wavefunc)
	site_count = Int(log(2,length(wavefunc)))
	first_elem = (transpose(wavefunc) * make_manybody_form(zz,site_count,sites) * wavefunc) / (transpose(wavefunc) * wavefunc)
	second_elem = get_local_magnetization(sites[1],wavefunc) * get_local_magnetization(sites[2],wavefunc)
	return first_elem - second_elem
end

function get_avg_correl_wdists(wavefunc)
	sites_count = Int(log(2,length(wavefunc)))
	all_corrs = [0.0*im for i in 1:sites_count]
	for i in 1:sites_count
		local_val = get_twosite_correlation([Int(ceil(sites_count/2)),i],wavefunc)
		all_corrs[i] = local_val
	end
	return all_corrs
end

function plot_site_mag_time_ev(site_vals,all_magns,all_times,site_count,time_steps)
	for i in 1:n
		for j in 1:time_steps+1
			xcomp,ycomp = get_arrow_shifts(magns[j][i])
			plot([site_vals[i],site_vals[i] + xcomp],[all_times[j],all_times[j] + ycomp],"-b")
		end
		plot([site_vals[i],site_vals[i]],[all_times[1],all_times[end]],"-k")
	end
	xlabel("Sites")
	ylabel("Time (x10)")
	return
end


org = [0,0,1,0,0]
n = length(org)
starting_wavefunc = make_tens_wavefunc_vec(get_wavefunc_givenorg(org))
#x_tens = turn_matrix_into_tensor(x)
#phi_form = randomITensor([Index(2) for i in 1:n])

#@einsum phi_form[b1,b2] := x_tens[b1,b1p] * given_wavefunc[b1p,b2]
#
final_time = 2.0
time_steps = 100
dt = final_time/time_steps
hx = 0.0
hz = 1.0
js = 2.0
ham = get_full_ham(n,js,hz,hx,dt)
#magns = [[0.0*im for i in 1:n] for j in 1:time_steps+1]
times = [10*(i-1)*final_time/time_steps for i in 1:time_steps+1]
twosite_correls = im.*zeros(n,time_steps+1)
#sites = [i for i in 1:n]
for i in 1:time_steps+1
	if i == 1
		global loc_wavefunc = starting_wavefunc
	end
	println(i)
	#=
	for j in 1:n
		next_magn = get_local_magnetization(j,loc_wavefunc)
		magns[i][j] = next_magn
	end
	=#
	rez = get_avg_correl_wdists(loc_wavefunc)
	twosite_correls[:,i] = rez
	global loc_wavefunc = ham * loc_wavefunc
end
#
imshow(real.(twosite_correls))
colorbar()




"fin"

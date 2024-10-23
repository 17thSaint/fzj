#####################################################
#=

This file is for calculating observables for effective 1D MPS

Depends on:
    

=#
######################################################

function make_density_correlations(wavefunc::MPS; kwargs...)
    Lphys::Int64,Lsynth::Int64 = get_mps_dims(wavefunc)

    density_correlations::Array{ComplexF64,4} = zeros(Lphys,Lphys,Lsynth,Lsynth)

    for s in 1:Lsynth
        for ss in 1:Lsynth
            corr_mat::Matrix{ComplexF64} = correlation_matrix(wavefunc,"Ns$(s)","Ns$(ss)")
            density_correlations[:,:,s,ss] = corr_mat
        end
    end

    return density_correlations
    
end

function get_occupancy(wavefunc::MPS; kwargs...)
	L,nflavors = get_mps_dims(wavefunc)
	if_squared = get(kwargs, :if_squared, false)
	remapping = kwargs[:remapping]
	
	if_plot = get(kwargs, :if_plot, true)
	if_3d = get(kwargs, :if_3d, false)
	if_save_data = get(kwargs, :if_save_data, false)
	if if_save_data
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "occs")
		filename = check_plot_label(filename,"occs")
		metadata = get(kwargs, :metadata, nothing)
	end
		
	occ_mat = zeros(L,nflavors)
	for s in 1:nflavors
		loc_op = "Ns$(s)"
		if if_squared
			loc_op = "Ns$(s) * Ns$(s)"
		end
		occ_mat[:, s] = expect(wavefunc, loc_op)[remapping]
	end
	if if_plot
		if if_3d
			plot_occupancy_3d(occ_mat; kwargs...)
		else
			plot_occupancy(occ_mat; kwargs...)
		end
	end
	data_dict = Dict([("vals",occ_mat)])
	if_save_data ? write_data_jld2(filename,data_dict,location,metadata) : nothing
	
	return occ_mat
end

function occupancy_variance(wavefunc::MPS; kwargs...)
	L,nflavors = get_mps_dims(wavefunc)
	fo_occ = get_occupancy(wavefunc; if_plot=false)
	occ_var = sqrt.(abs.(fo_occ.^2 .- get_occupancy(wavefunc; if_plot=false,if_squared=true)))
	if_plot = get(kwargs, :if_plot, true)
	if_plot ? plot_occupancy(occ_var; kwargs...) : nothing
	
	return occ_var
end

function entanglement_spectrum(psi::MPS, bond::Int)
    # Split the MPS at the bond using svd
    l, s, r = svd(psi[bond], linkind(psi, bond))

    # Form the reduced density matrix of the left part
    rho = l * s * dag(l)

    # Diagonalize the reduced density matrix to get the eigenvalues
    evals = eigvals(matrix(rho))

    # The entanglement spectrum is given by the negative logarithm of these eigenvalues
    spectrum = evals

    return spectrum
end

function entanglement_entropy(psi::MPS)
	# Get the entanglement spectrum
	spectrum = entanglement_spectrum(psi, Int(floor(length(psi)/2)))

	# Form the entanglement entropy by summing the spectrum
	entropy = -sum(spectrum .* log.(spectrum))

	return abs.(entropy)
end

function e_entropy(M::MPS, l::Int)
        
	# get length of MPS
	L = length(M)
	
	# catch errors where the cut is "outside of the system"
	if l < 1 || l >= L
		return 0.0
	end
	
	# shifting the center of orthogonality to site l
	ITensors.orthogonalize!(M, l)
	
	# perform a singular value decomposition / Schmidt decomposition
	U,S,V = svd(M[l], (linkind(M, l)))#, siteind(M, l)))
	
	# set entropy to zero
	res = 0.0
	
	# loop over all diagonal elements of S
	for i in 1:size(S, 1)
		res += -2.0 * S[i,i]^2 * log(S[i,i]) # = - S² * log(S²)
	end

	return res
end

function density_matrix(wavefunc::MPS; kwargs...)
	L::Int64,nflavors::Int64 = get_mps_dims(wavefunc)
	densmat::Matrix{ComplexF64} = zeros(ComplexF64,L*nflavors,L*nflavors)
	for s in 1:nflavors
		for sp in 1:nflavors
			densmat[L*(s-1)+1:L*s,L*(sp-1)+1:L*sp] = correlation_matrix(wavefunc,"Cr$(s)","Anh$(sp)")
		end
	end
	return densmat
end

# not sure this is actually being used
function normalize_densmat(dens_mat::Matrix,part_count::Int; kwargs...)
    println("Using normalize_densmat function")
	if_log = get(kwargs, :if_log, false)
	L = size(dens_mat)[1]
	current_trace = if_log ? log_sum(diag(dens_mat)) : tr(dens_mat)
	if if_log
		shift_mat = Diagonal([log(part_count) - current_trace for i in 1:L])
		norm_densmat = dens_mat + shift_mat
	else
		shift_mat = ones(L,L)
		for i in 1:L
			shift_mat[i,i] = part_count/current_trace
		end
		norm_densmat = dens_mat .* shift_mat
	end
	return norm_densmat
end

function hopping_correlation(wavefunc::MPS, hopping_direction::String; kwargs...)
    Lx::Int64,Ly::Int64 = get_mps_dims(wavefunc)
	start_point::Int = 1
	
	all_greens::Matrix{ComplexF64} = zeros(ComplexF64,Ly,Lx)
	for s in 1:Ly
		all_greens[s,:] = hopping_direction == "virt" ? diag(ITensors.correlation_matrix(wavefunc,"Cr$(start_point)","Anh$(s)")) : ITensors.correlation_matrix(wavefunc,"Cr$(s)","Anh$(s)")[start_point,:]
	end
	
	start_norm::Vector{ComplexF64} = hopping_direction == "virt" ? ITensors.expect(wavefunc,"Ns$(start_point)") : [ITensors.expect(wavefunc,"Ns$(i)";sites=start_point) for i in 1:Ly]
	all_norms::Matrix{ComplexF64} = zeros(ComplexF64,Ly,Lx)
	for s in 1:Ly
		const_part::Vector{ComplexF64} = hopping_direction == "virt" ? start_norm : [start_norm[s] for i in 1:Lx]
		all_norms[s,:] = ITensors.expect(wavefunc,"Ns$(s)") .* const_part
	end
	all_greens ./= sqrt.(all_norms)
	all_greens = abs.(all_greens)
	
	if_plot::Bool = get(kwargs, :if_plot, true)
	if_plot ? plot_greenfunc(real.(all_greens),hopping_direction; kwargs...) : nothing
	
    # needs be to checked
	#=if_save_data = get(kwargs, :if_save_data, false)
	if if_save_data
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "$hopping_direction-dir-GF")
		filename = check_plot_label(filename,"$hopping_direction-dir-GF")
		metadata = get(kwargs, :metadata, nothing)
		data_dict = Dict([("vals",all_greens)])
	end
	if_save_data ? write_data_jld2(filename,data_dict,location,metadata) : nothing=#
	
	
	return real.(all_greens)
end

function hopping_correlation(densmat::Matrix{ComplexF64}, hopping_direction::String, Lx::Int64, Ly::Int64; kwargs...)
    error("Not implemented yet")
end

    










































"fin"
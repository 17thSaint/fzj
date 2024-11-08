#####################################################
#=

This file is for calculating observables for effective 1D MPS

Depends on:
    

=#
######################################################

function get_local_nrgs(metadata::Dict)
    number_levels::Int64 = length(filter(x -> occursin("observer",x),keys(metadata)))
    local_nrgs::Vector{Float64} = zeros(Float64,number_levels)
    for i in 1:number_levels
        obs_string::String = i == 1 ? "observer" : "observer_$(i-1)"
        local_nrgs[i] = metadata[obs_string].energies[end]
    end
    return local_nrgs
end

function get_correct_gs_mpsstring(metadata::Dict)
	local_nrgs = get_local_nrgs(metadata)
	if local_nrgs[1] != sort(local_nrgs)[1]
		println("Levels are not sorted")
		which_level_gs = findfirst(x -> minimum(local_nrgs) == local_nrgs[x],1:length(local_nrgs))
		wavefunc_string = "mps_$(which_level_gs-1)"
	else
		which_level_gs = 1
		wavefunc_string = "mps"
	end
	return wavefunc_string
end

# maybe zero shift is needed
function make_density_correlations(wavefunc::MPS; kwargs...)
	if_save::Bool = get(kwargs,:if_save,false)
	if_zero_shift::Bool = get(kwargs,:if_zero_shift,true)
    Lphys::Int64,Lsynth::Int64 = get_mps_dims(wavefunc)

    density_correlations::Array{Float64,4} = zeros(Lphys,Lphys,Lsynth,Lsynth)

    for s in 1:Lsynth
        for ss in 1:Lsynth
            corr_mat::Matrix{Float64} = real.(correlation_matrix(wavefunc,"Ns$(s)","Ns$(ss)"))
            density_correlations[:,:,s,ss] = corr_mat
			if if_zero_shift
				density_correlations[:,:,s,ss] -= real.(expect(wavefunc,"Ns$(s)")) * transpose(real.(expect(wavefunc,"Ns$(ss)")))
			end
        end
    end

	if if_save
		filepath = kwargs[:filepath]
		save_dd_correlation(density_correlations,filepath)
	end

    return density_correlations
    
end

function ft_densitydensity_correlation(momentum_angle::Float64,momentum_radius::Float64,wavefunc::Union{Nothing,MPS}; kwargs...)
    denscorrs = get(kwargs,:denscorrs,nothing)
	if_save::Bool = get(kwargs,:if_save,false)

	if isnothing(denscorrs) && isnothing(wavefunc)
		error("Either denscorrs or wavefunc must be provided")
	end

    if isnothing(denscorrs)
        denscorrs = make_density_correlations(wavefunc; kwargs...)
    end

    momentum = momentum_radius .* [cos(momentum_angle),sin(momentum_angle)]

    all_distances::Array{Float64,4} = zeros(Float64,size(denscorrs))
    for j in 1:size(denscorrs)[1]
        for jj in 1:size(denscorrs)[2]
            for s in 1:size(denscorrs)[3]
                for ss in 1:size(denscorrs)[4]
                    all_distances[j,jj,s,ss] = dot(momentum,[j-jj,s-ss])
                end
            end
        end
    end

    result::ComplexF64 = sum(denscorrs .* exp.(im .* all_distances)) / ((size(denscorrs,1) * size(denscorrs,3))^2)

    if if_save
        filepath = kwargs[:filepath]
        save_ft_dd(result,round(momentum_angle/pi,digits=3),filepath)
    end

    return result
end

function ft_densitydensity_correlation(momentum::Vector{Float64},wavefunc::Union{Nothing,MPS}; kwargs...)
	if momentum == [0.0,0.0]
		return ft_densitydensity_correlation(0.0,0.0,wavefunc; kwargs...)
	else
    	return ft_densitydensity_correlation(atan(momentum[2]/momentum[1]),sqrt(momentum[1]^2 + momentum[2]^2),wavefunc; kwargs...)
	end
end

function findall_ft_dd_1deff(lx::Int64,ly::Int64,n::Int64; kwargs...)
    hanis::Float64 = get(kwargs,:hanis,1.0)

	ks = range(0,stop=2*pi,length=51)

    pdict = Dict([("Lphys",lx),("Lsynth",ly),("nbosons",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    all_files = find_data_file(pdict,"mps",dataloc; output_level=0)

    filter!(x -> !occursin("twist_angle1",x),all_files)
    filter!(x -> !occursin("mk",x),all_files)

	if length(all_files) != 1
		error("Either too many files or none")
	end

	f = all_files[1]

	filepath = dataloc * "/" * f
	d,m = read_data_jld2(filepath; output_level=0)

	if !haskey(m,"dens_corr_mat")
		error("No density correlation matrix found")
	end

	println("Working on $(lx)x$(ly) n=$(n)")

	all_ft_vals::Matrix{Float64} = zeros(Float64,length(ks),length(ks))
	for (idx,kx) in enumerate(ks)
		for (idx2,ky) in enumerate(ks)
			all_ft_vals[idx,idx2] = abs(ft_densitydensity_correlation([kx,ky],nothing; denscorrs=m["dens_corr_mat"]))
		end
	end

	normalization_factor::Float64 = integrate_2d_matrix(all_ft_vals)

	result = maximum(all_ft_vals) / normalization_factor

    return result
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

# measures the flatness parameter which is max E2 - E1 / E3 - E1
function twist_flatness_1deff(lx::Int,ly::Int,n::Int; kwargs...)
    hanis = get(kwargs,:hanis,1.0)
	if_plot_spectrum::Bool = get(kwargs,:if_plot_spectrum,false)
    dataloc = get_folder_location("cluster-data/synth-dims/twists")
    
    all_flatnesses = []

    params_dict = Dict([("Lphys",lx),("Lsynth",ly),("nbosons",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",hanis)])
    all_files = find_data_file(params_dict,"mps",dataloc; output_level=0)

	tw1s = Float64[]
	tw2s = Float64[]
	all_nrgs = Dict([("1",Float64[]),("2",Float64[]),("3",Float64[])])

    for f in all_files

		filename_dict = get_params_dict_from_filename(f)
		if haskey(filename_dict,"tw1") && (filename_dict["tw1"] == 0.67 || filename_dict["tw1"] == 0.33 || filename_dict["tw2"] == 0.67 || filename_dict["tw2"] == 0.33)
			continue
		end

		#=if (lx,ly,n) == (8,3,3) && (filename_dict["tw2"] == 0.0 || filename_dict["tw2"] == 1.0)
			continue
		end=#

		if (lx,ly,n) == (4,6,3) && (filename_dict["tw2"] == 0.0 || filename_dict["tw2"] == 1.0)
            continue
        end

		if (lx,ly,n) == (8,5,5) && (filename_dict["tw2"] == 0.0 || filename_dict["tw2"] == 1.0)
            continue
        end

        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

		if haskey(m,"observer") && haskey(m,"observer_1") && haskey(m,"observer_2")
			flatness = (m["observer_1"].energies[end] - m["observer"].energies[end]) / (m["observer_2"].energies[end] - m["observer"].energies[end])
			append!(all_flatnesses,[flatness])
			append!(tw1s,[filename_dict["tw1"]])
			append!(tw2s,[filename_dict["tw2"]])
			local_nrgs = zeros(Float64,3)
			for i in 1:3
				obs_string = i == 1 ? "observer" : "observer_$(i-1)"
				local_nrgs[i] = m[obs_string].energies[end]
			end
			sort!(local_nrgs)
			for i in 1:3
				append!(all_nrgs[string(i)],[local_nrgs[i]])
			end
		else
			println("File $f doesn't have all states")
			#display(keys(m))
		end
    end

	if_plot_spectrum ? plot_twisting_spectrum(tw1s,tw2s,all_nrgs; kwargs...) : nothing

	filter!(x->x >= 0.0 && x <= 1.0,all_flatnesses)

    return maximum(all_flatnesses)
end










































"fin"
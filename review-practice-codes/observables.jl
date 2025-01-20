#####################################################
#=

This file contains observable calculations for TTNs

Depends on:
    review-practice-codes/ttn.jl

=#
######################################################

include("../other-funcs/include-other-files.jl")

include_other_files(["review-practice-codes/ttn.jl"])

function get_occupancy(ttn::TTNKit.TreeTensorNetwork; kwargs...)
	densmat = get(kwargs, :densmat, nothing)

	if isnothing(densmat)
		exp_occ = abs.(TTNKit.expect(ttn,"N"))
	else
		lat = TTNKit.physical_lattice(TTNKit.network(ttn))
		phys_length,virt_length = get_lattice_dims(ttn)
		exp_occ = zeros(phys_length,virt_length)
		for j in 1:phys_length
			for s in 1:virt_length
				linear_index = TTNKit.linear_ind(lat,(j,s))
				exp_occ[j,s] = abs(densmat[linear_index,linear_index])
			end
		end
	end
	
	if_save_data = get(kwargs, :if_save_data, false)
	if_save_fig = get(kwargs, :if_save_fig, false)
	if_plot = get(kwargs, :if_plot, true)
	
	if_save_data ? save_occupancy(exp_occ; kwargs...) : nothing
	if_plot	|| if_save_fig ? plot_occupancy(exp_occ; kwargs...) : nothing
		
	return exp_occ
end

function get_occupancy(densmat::Matrix; kwargs...)
	if isinteger(sqrt(size(densmat)[1]))
		phys_length = Int(sqrt(size(densmat)[1]))
		virt_length = Int(phys_length)
	else
		phys_length = Int(sqrt(2*size(densmat)[1]))
		virt_length = Int(phys_length/2)
	end

	lat = TTNKit.SimpleLattice((phys_length,virt_length),TTNKit.ITensorNode,"Boson")

	exp_occ = zeros(phys_length,virt_length)
	for j in 1:phys_length
		for s in 1:virt_length
			linear_index = TTNKit.linear_ind(lat,(j,s))
			exp_occ[j,s] = abs(densmat[linear_index,linear_index])
		end
	end

    if_synth_rectangle::Bool = get(kwargs,:if_synth_rectangle,false)
    if_synth_rectangle ? exp_occ = transpose(exp_occ) : nothing

	if_save_data = get(kwargs, :if_save_data, false)
	if_save_fig = get(kwargs, :if_save_fig, false)
	if_plot = get(kwargs, :if_plot, true)

    if_save_data ? save_occupancy(exp_occ; kwargs...) : nothing
	if_plot	|| if_save_fig ? plot_occupancy(exp_occ; kwargs...) : nothing

	return exp_occ
end

function save_occupancy(exp_occ; kwargs...)
	location = get(kwargs, :location, pwd())
	filename = get(kwargs, :name, "occs")
	metadata = get(kwargs, :metadata, nothing)
	occs_data_dict = Dict([("vals",exp_occ)])
	write_data_jld2(filename,occs_data_dict,location,metadata)
	return
end

function physical_correlation(densmat::Matrix{ComplexF64},Lx::Int64,Ly::Int64; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    phys_corrs = Array{Float64,2}(undef,Lx-1,Ly)
    for j in 2:Lx
        for s in 1:Ly
            phys_corr = densmat[linear_index((j,s),Lx,Ly),linear_index((1,s),Lx,Ly)]
            phys_corr /= sqrt(densmat[linear_index((j,s),Lx,Ly),linear_index((j,s),Lx,Ly)] * densmat[linear_index((1,s),Lx,Ly),linear_index((1,s),Lx,Ly)])
            phys_corrs[j-1,s] = abs(phys_corr)
        end
    end

    if_plot ? plot_physical_correlation(phys_corrs; kwargs...) : nothing

    return phys_corrs
end

function physical_correlation(wavefunc::TTNKit.TreeTensorNetwork; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
    densmat = get(kwargs,:densmat,nothing)
	if !isnothing(densmat)
        lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))
	    Lx,Ly = size(lat)
        return physical_correlation(densmat,Lx,Ly; kwargs...)
    else
        densmat = density_matrix(wavefunc)
    end

	lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))
	Lx,Ly = size(lat)

    phys_corrs = Array{Float64,2}(undef,Lx-1,Ly)
    for j in 2:Lx
        for s in 1:Ly
            phys_corr = densmat[TTNKit.linear_ind(lat,(j,s)),TTNKit.linear_ind(lat,(1,s))]
            phys_corr /= sqrt(densmat[TTNKit.linear_ind(lat,(j,s)),TTNKit.linear_ind(lat,(j,s))] * densmat[TTNKit.linear_ind(lat,(1,s)),TTNKit.linear_ind(lat,(1,s))])
            phys_corrs[j-1,s] = abs(phys_corr)
        end
    end

    if_plot ? plot_physical_correlation(phys_corrs; kwargs...) : nothing

    return phys_corrs
end

function synthetic_correlation(densmat::Matrix{ComplexF64},Lx::Int64,Ly::Int64; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    syn_corrs = Array{Float64,2}(undef,Lx,Ly-1)
    for j in 1:Lx
        for s in 2:Ly
			syn_corr = densmat[linear_index((j,s),Lx,Ly),linear_index((j,1),Lx,Ly)]
			syn_corr /= sqrt(densmat[linear_index((j,s),Lx,Ly),linear_index((j,s),Lx,Ly)] * densmat[linear_index((j,1),Lx,Ly),linear_index((j,1),Lx,Ly)])
			syn_corrs[j,s-1] = abs(syn_corr)
        end
    end

    if_plot ? plot_synthetic_correlation(syn_corrs; kwargs...) : nothing

    return syn_corrs
end

function synthetic_correlation(wavefunc::TTNKit.TreeTensorNetwork; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
	densmat = get(kwargs,:densmat,nothing)
	if isnothing(densmat)
		densmat = density_matrix(wavefunc)
    else
        lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))
        Lx,Ly = size(lat)
        return synthetic_correlation(densmat,Lx,Ly; kwargs...)
	end

	lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))
	Lx,Ly = size(lat)

    syn_corrs = Array{Float64,2}(undef,Lx,Ly-1)
    for j in 1:Lx
        for s in 2:Ly
			syn_corr = densmat[TTNKit.linear_ind(lat,(j,s)),TTNKit.linear_ind(lat,(j,1))]
			syn_corr /= sqrt(densmat[TTNKit.linear_ind(lat,(j,s)),TTNKit.linear_ind(lat,(j,s))] * densmat[TTNKit.linear_ind(lat,(j,1)),TTNKit.linear_ind(lat,(j,1))])
			syn_corrs[j,s-1] = abs(syn_corr)
        end
    end

    if_plot ? plot_synthetic_correlation(syn_corrs; kwargs...) : nothing

    return syn_corrs
end

function physical_current(densmat::Matrix{ComplexF64},Lx::Int64,Ly::Int64; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    currents = Array{Float64,2}(undef,Lx,Ly)
    for s in 1:Ly
        for j in 1:Lx
            site1 = linear_index((j,s),Lx,Ly)
            site2 = linear_index((mod1(j+1,Lx),s),Lx,Ly)
            current_val = imag(densmat[site1,site2] - densmat[site2,site1])
            current_normalization = densmat[site1,site1] + densmat[site2,site2]
            current_val /= current_normalization
            currents[j,s] = real(current_val)
        end
    end

    if_plot ? plot_physical_current(currents; kwargs...) : nothing

    return currents
end

function physical_current(wavefunc::TTNKit.TreeTensorNetwork; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
	densmat = get(kwargs,:densmat,nothing)
	if isnothing(densmat)
		densmat = density_matrix(wavefunc)
    else
        lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))
        Lx,Ly = size(lat)
        return physical_current(densmat,Lx,Ly; kwargs...)
	end

    lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))
	Lx,Ly = size(lat)

    currents = Array{Float64,2}(undef,Lx,Ly)
    for s in 1:Ly
        for j in 1:Lx
            site1 = TTNKit.linear_ind(lat,(j,s))
            site2 = TTNKit.linear_ind(lat,(mod1(j+1,Lx),s))
            current_val = imag(densmat[site1,site2] - densmat[site2,site1])
            current_normalization = densmat[site1,site1] + densmat[site2,site2]
            current_val /= current_normalization
            currents[j,s] = real(current_val)
        end
    end

    if_plot ? plot_physical_current(currents; kwargs...) : nothing

    return currents
end

function synthetic_current(densmat::Matrix{ComplexF64},Lx::Int64,Ly::Int64; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    currents = Array{Float64,2}(undef,Lx,Ly)
    for j in 1:Lx
        for s in 1:Lx
            site1 = linear_index((j,s),Lx,Ly)
            site2 = linear_index((j,mod1(s+1,Ly)),Lx,Ly)
            current_val = imag(densmat[site1,site2] - densmat[site2,site1])
            current_normalization = densmat[site1,site1] + densmat[site2,site2]
            current_val /= current_normalization
            currents[j,s] = real(current_val)
        end
    end

    if_plot ? plot_synthetic_current(currents; kwargs...) : nothing

    return currents
end

function synthetic_current(wavefunc::TTNKit.TreeTensorNetwork; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
	densmat = get(kwargs,:densmat,nothing)
	if isnothing(densmat)
		densmat = density_matrix(wavefunc)
    else
        lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))
        Lx,Ly = size(lat)
        return synthetic_current(densmat,Lx,Ly; kwargs...)
	end

	lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))
	Lx,Ly = size(lat)

    currents = Array{Float64,2}(undef,Lx,Ly)
    for j in 1:Lx
        for s in 1:Lx
            site1 = TTNKit.linear_ind(lat,(j,s))
            site2 = TTNKit.linear_ind(lat,(j,mod1(s+1,Ly)))
            current_val = imag(densmat[site1,site2] - densmat[site2,site1])
            current_normalization = densmat[site1,site1] + densmat[site2,site2]
            current_val /= current_normalization
            currents[j,s] = real(current_val)
        end
    end

    if_plot ? plot_synthetic_current(currents; kwargs...) : nothing

    return currents
end

function check_nrg_convergence(metadata::Dict,if_perfect::Bool=true)
    nrg_tol::Float64 = metadata["nrgtol"]

    if_perfect ? nothing : nrg_tol *= 10

    all_nrgs::Dict{String,Array{Float64,1}} = Dict()
    all_nrgs["0"] = metadata["observer"].nrg

    found_excited_nrg::Bool = haskey(metadata,"observer_1")
    next_nrg_level::Int = 1
    while found_excited_nrg
        all_nrgs[string(next_nrg_level)] = metadata["observer_$(next_nrg_level)"].nrg
        next_nrg_level += 1
        found_excited_nrg = haskey(metadata,"observer_$(next_nrg_level)")
    end

    if_converged::Dict{String,Bool} = Dict()
    for (k,v) in all_nrgs
        if_converged[k] = abs(v[end] - v[end-1]) < nrg_tol
    end

    return all(values(if_converged)),if_converged
end

# find some way to speed this up
function make_density_correlations(wavefunc::TTNKit.TreeTensorNetwork; kwargs...)
    Lphys::Int64,Lsynth::Int64 = get_lattice_dims(wavefunc)

    density_correlations::Array{Float64,4} = zeros(Lphys,Lphys,Lsynth,Lsynth)

    for s in 1:Lsynth
        for j in 1:Lphys
            for ss in 1:Lsynth
                for jj in 1:Lphys
                    corr_val::Float64 = real(TTNKit.correlation(wavefunc,"N","N",(j,s),(jj,ss)))
                    density_correlations[j,jj,s,ss] = corr_val
                end
            end
        end
    end

    return density_correlations
end

function get_cdwsf(qvec::Vector{Float64},dens_corr_mat::Array{Float64}; kwargs...)
    Lphys::Int64,Lsynth::Int64 = size(dens_corr_mat)[1],size(dens_corr_mat)[3]

    result::Union{Float64,ComplexF64} = 0.0
    
    for j in 1:Lphys
        for s in 1:Lsynth
            for jj in 1:Lphys
                for ss in 1:Lsynth
                    dist_vect::Vector{Int64} = [minimum([(j-jj),Lphys-(j-jj)]),minimum([(s-ss),Lsynth-(s-ss)])]
                    dotprod::Float64 = dot(qvec,dist_vect)
                    result += exp(-im*2*pi*dotprod)*dens_corr_mat[j,jj,s,ss]
                    #println("Working on sites $(j) and $(jj) with synth indices $(s) and $(ss) with dotprod $(dotprod)")
                end
            end
        end
    end

    return result / ((Lphys*Lsynth)^2)
end

function range_cdwsf_angles(points_count::Int64,dens_corr_mat::Array{Float64},radius::Int64=1.0; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)
    plot_title::String = get(kwargs,:plot_title,"")
    if_multiple_lines::Bool = get(kwargs,:if_multiple_lines,false)
    line_label::String = get(kwargs,:line_label,"")

    angles = range(0.0,2*pi,length=points_count)
    cdwsfs::Vector{ComplexF64} = zeros(ComplexF64,points_count)
    for (idx,angle) in enumerate(angles)
        qvec::Vector{Float64} = [radius*cos(angle),radius*sin(angle)]
        cdwsfs[idx] = get_cdwsf(qvec,denscorrs)
    end

    if if_plot
        if if_multiple_lines
            plot(angles ./ (2*pi),abs.(cdwsfs),"-p",label=line_label)
            legend()
        else
            plot(angles ./ (2*pi),abs.(cdwsfs),"-p")
        end
        xlabel("Angle")
        ylabel("CDW Structure Factor")
        title(plot_title)
    end

    return angles,cdwsfs
end

function pointfourhopping(psi::TTNKit.TreeTensorNetwork,s1::Int64,s2::Int64)
    result = TTNKit.correlation(psi,"N","N",s1,s2)

    if s1 == s2
        result -= TTNKit.expect(psi,"N",s1)
    end

    return result
end

# to be used when using smaller lattice than TTN actual size
function fourpoint_alberto(psi::TTNKit.TreeTensorNetwork,restricted_size::Vector{Int64}; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)

    ttn_lx,ttn_ly = get_lattice_dims(psi)
    r_lx,r_ly = restricted_size

    center_site::Vector{Int64} = [Int64(ceil(r_lx/2)),Int64(ceil(r_ly/2))]
    center_linear::Int64 = linear_index(center_site,ttn_lx,ttn_ly)

    rez::Matrix{Float64} = zeros(Float64,r_ly,r_lx)
    for j in 1:r_lx
        for s in 1:r_ly
            println("Working on site ($(j), $(s))")
            rez[s,j] = real(pointfourhopping(psi,center_linear,linear_index((j,s),ttn_lx,ttn_ly)))
        end
    end

    if_plot ? plot_fourpointcorrelator(rez; kwargs...) : nothing

    return rez
end

function fourpoint_alberto(psi::TTNKit.TreeTensorNetwork; kwargs...)
    if_restricted_size::Bool = get(kwargs,:if_restricted_size,false)
    if_restricted_size && fourpoint_alberto(psi,if_restricted_size; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)

    lx,ly = get_lattice_dims(psi)

    center_site::Vector{Int64} = [Int64(ceil(lx/2)),Int64(ceil(ly/2))]
    center_linear::Int64 = linear_index(center_site,lx,ly)

    rez::Matrix{Float64} = zeros(Float64,ly,lx)
    for j in 1:lx
        for s in 1:ly
            println("Working on site ($(j), $(s))")
            rez[s,j] = real(pointfourhopping(psi,center_linear,linear_index((j,s),lx,ly)))
        end
    end

    if_plot ? plot_fourpointcorrelator(rez; kwargs...) : nothing

    return rez
end

function pairdistribution(psi::TTNKit.TreeTensorNetwork; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)
    occs::Union{Nothing,Matrix{Float64}} = get(kwargs,:occs,nothing)
    if isnothing(occs)
        occs = get_occupancy(psi; kwargs...,if_plot=false)
    end

    Lsynth::Int64,Lphys::Int64 = size(occs)

    fourpoint::Matrix{Float64} = fourpoint_alberto(psi; kwargs...)
    centersite::Vector{Int64} = [Int64(ceil(Lphys/2)),Int64(ceil(Lsynth/2))]
    pairdist::Matrix{Float64} = fourpoint ./ (occs[centersite[2],centersite[1]] .* occs)

    if_plot ? plot_pairdistribution(pairdist; kwargs...) : nothing

    return pairdist
end

function spatial_entanglement_spectrum(psi::TTNKit.TreeTensorNetwork; kwargs...)
    if_save::Bool = get(kwargs,:if_save,false)
    filepath::Union{String,Nothing} = get(kwargs,:filepath,nothing)

    numlayers = TTNKit.number_of_layers(psi)
    TTNKit.move_ortho!(psi,(numlayers,1))
    top_tensor = psi[numlayers,1]
    idx_left = inds(top_tensor; tags = "Link,nl=$(numlayers-1),np=1")
    u,s,v,spec = svd(top_tensor,idx_left)

    if_save && save_spatialentanglementspectrum(spec.eigs[1:100],filepath)

    return spec.eigs[1:100]
end

function save_spatialentanglementspectrum(spec::Vector{Float64},filepath::Nothing)
    error("File Path not provided")
end

function save_spatialentanglementspectrum(spec::Vector{Float64},filepath::String)
    modify_data_jld2(Dict([("entanglement_spectrum",spec)]),filepath,"metadata"; output_level=0)
end








































"fin"
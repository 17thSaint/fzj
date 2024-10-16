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









































"fin"
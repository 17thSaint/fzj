#####################################################
#=

This file contains the simple observable functions for ED

Depends on:

=#
######################################################


# functions for onsite occupancy n_i from either wavefunction or density matrix
function get_occupancy(x::Vector{ComplexF64},lattice_params::Dict{String,Any}; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    full_basis = lattice_params["full_basis"]

    occs = zeros(Float64,Ly,Lx)
    
    for i in 1:length(x)
        basis_state = full_basis[:,i]
        for n in 1:length(basis_state)
            site = coordinate(basis_state[n],Lx,Ly)
            occs[site[2],site[1]] += abs(x[i])^2
        end
    end

    if_plot ? plot_occupancy(occs; kwargs...) : nothing

    return occs
end

function get_occupancy(rho::Array{ComplexF64,2},lattice_params::Dict{String,Any}; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

    occs = zeros(Float64,Ly,Lx)
    all_occs = diag(rho)

    for (idx,val) in enumerate(all_occs)
        site = coordinate(idx,Lx,Ly)
        occs[site[2],site[1]] = real(val)
    end

    if_plot ? plot_occupancy(occs; kwargs...) : nothing

    return occs
end

function plot_occupancy(exp_occ; kwargs...)
	fig = figure()
	imshow(exp_occ)
	colorbar()
	plot_title = get(kwargs, :plot_title, "")
	title_string = "Occupancy, " * plot_title
	title(title_string)
	ylabel("Synthetic")
	xlabel("Physical")

    return nothing
end


# functions for two-site correlation along physical or synthetic dimension, made from density matrix
function physical_correlation(densmat::Array{ComplexF64,2},Lx::Int64,Ly::Int64; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    phys_corrs = Array{Float64,2}(undef,Lx,Ly)
    for j in 1:Lx
        for s in 1:Ly
            phys_corr = densmat[linear_index((j,s),Lx,Ly),linear_index((1,s),Lx,Ly)]
            phys_corr /= sqrt(densmat[linear_index((j,s),Lx,Ly),linear_index((j,s),Lx,Ly)] * densmat[linear_index((1,s),Lx,Ly),linear_index((1,s),Lx,Ly)])
            phys_corrs[j,s] = abs(phys_corr)
        end
    end

    if_plot ? plot_physical_correlation(phys_corrs; kwargs...) : nothing

    return phys_corrs
end

function plot_physical_correlation(phys_corrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(phys_corrs)[2]
        plot(0:size(phys_corrs)[1]-1,phys_corrs[:,i],"-p",label="$i")
    end
    xlabel("Physical Distance")
    ylabel("Correlation")
    title("Physical Correlation "*plot_title)
    legend()
    yscale("log")
    return nothing
end

function synthetic_correlation(densmat::Array{ComplexF64,2},Lx::Int64,Ly::Int64; kwargs...)
    if_plot = get(kwargs,:if_plot,true)

    syn_corrs = Array{Float64,2}(undef,Lx,Ly)
    for j in 1:Lx
        for s in 1:Ly
            syn_corr = densmat[linear_index((j,s),Lx,Ly),linear_index((j,1),Lx,Ly)]
            syn_corr /= sqrt(densmat[linear_index((j,s),Lx,Ly),linear_index((j,s),Lx,Ly)] * densmat[linear_index((j,1),Lx,Ly),linear_index((j,1),Lx,Ly)])
            syn_corrs[j,s] = abs(syn_corr)
        end
    end

    if_plot ? plot_synthetic_correlation(syn_corrs; kwargs...) : nothing

    return syn_corrs
end

function plot_synthetic_correlation(syn_corrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(syn_corrs)[1]
        plot(0:size(syn_corrs)[2]-1,syn_corrs[i,:],"-p",label="$i")
    end
    xlabel("Synthetic Distance")
    ylabel("Correlation")
    title("Synthetic Correlation "*plot_title)
    legend()
    yscale("log")
    return nothing
end


# functions for current along physical or synthetic dimension, made from density matrix
function physical_current(densmat::Array{ComplexF64,2},lattice_params::Dict{String,Any}; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    if_periodic_x = lattice_params["if_periodic_x"]

    xlen = if_periodic_x ? Lx : Lx-1

    currents = Array{Float64,2}(undef,xlen,Ly)
    for s in 1:Ly
        for j in 1:xlen
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

function plot_physical_current(currents::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(currents)[1]
        plot(1:size(currents)[2],currents[i,:],"-p",label="$i")
    end
    xlabel("Synthetic Site")
    ylabel("Current")
    title("Physical Current "*plot_title)
    legend()
    return nothing
end

function synthetic_current(densmat::Array{ComplexF64,2},lattice_params::Dict{String,Any}; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]
    if_periodic_y = lattice_params["if_periodic_y"]

    ylen = if_periodic_y ? Ly : Ly-1

    currents = Array{Float64,2}(undef,Lx,ylen)
    for j in 1:Lx
        for s in 1:ylen
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

function plot_synthetic_current(currents::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    for i in 1:size(currents)[2]
        plot(1:size(currents)[1],currents[:,i],"-p",label="$i")
    end
    xlabel("Physical Site")
    ylabel("Current")
    title("Synthetic Current "*plot_title)
    legend()
    return nothing
end


# gap from Cecille Repellin's paper on Harper-Hofstadter model by hopping anisotropy
function hh_gap_exact(hopping_anisotropy::Float64,alpha::Float64)
    a2 = 1/(4*(sin(pi*alpha)^2))
    a1 = (1+2*(1/hopping_anisotropy^2)*(a2^2))^-0.5
    val = 4*alpha*(1/hopping_anisotropy^2)*(a1^4)*(a2^2)
    if val > 0.5
        return 1/(1+2*(1/hopping_anisotropy^2)*(a2^2))
    else
        return val
    end
end

function hh_gap_fit(x,p)
    return p[1] .* hh_gap_exact.(x,0.25) .+ p[2]
end


# correlations as a function of distance for along a given direction
# needs testing, maybe redundant from physical_correlation and synthetic_correlation
function distance_correlation(rho::Matrix,lattice_params::Dict,direction::String="x")
    if_periodic_x = lattice_params["if_periodic_x"]
    if_periodic_y = lattice_params["if_periodic_y"]
    Lx = lattice_params["Lx"]
    Ly = lattice_params["Ly"]

    if direction == "x"
        len = Lx
        other_len = Ly
    else
        len = Ly
        other_len = Lx
    end
    dist_corrs = zeros(Float64,len-1)
    corr_counts = zeros(Int64,len-1)

    for x1 in 1:len
        for y1 in 1:other_len
            s1 = direction == "x" ? linear_index((x1,y1),Lx,Ly) : linear_index((y1,x1),Lx,Ly)
            for x2 in 1:len-1
                if x1 == x2
                    continue
                end
                s2 = direction == "x" ? linear_index((x2,y1),Lx,Ly) : linear_index((y1,x2),Lx,Ly)
                dist_corrs[Int(abs(x1-x2))] += abs(rho[s1,s2])
                corr_counts[Int(abs(x1-x2))] += 1
            end
        end
    end

    dist_corrs ./= corr_counts

    return dist_corrs
end

distance_correlation(x::Vector{ComplexF64},lattice_params::Dict,direction::String="x") = distance_correlation(density_matrix(x,lattice_params),lattice_params,direction)

# Charge Density Wave structure factor as a function of direction given by qvec
# needs testing
function cdw_structure_factor(rho::Array{ComplexF64,2},qvec::Vector,psi::Vector{ComplexF64},lattice_params::Dict; kwargs...)
    if_periodic_phys = lattice_params["if_periodic_x"]
    if_periodic_synth = lattice_params["if_periodic_y"]
    phys_len = lattice_params["Lx"]
    synth_len = lattice_params["Ly"]


	occs = get_occupancy(psi,lattice_params; if_plot=false,densmat=rho)

	struc_fact = 0.0
	for j in 1:phys_len
		for s in 1:synth_len
			p1 = (j,s)
			p1_linear = linear_index(p1,phys_len,synth_len)
			for jj in 1:phys_len
				for ss in 1:synth_len
					p2 = (jj,ss)
					p2_linear = linear_index(p2,phys_len,synth_len)
					dist = find_dist(p1, p2, (phys_len,synth_len), (if_periodic_phys,if_periodic_synth))[2]
					struc_fact += occs[p1[1],p1[2]] * occs[p2[1],p2[2]] * exp(im * dot(qvec,dist))
				end
			end
		end
	end
	return struc_fact / sum(occs)
end

# needs testing
function cdw_structure_factor(occs::Matrix,qvec::Vector,lattice_params::Dict)
    if_periodic_phys = lattice_params["if_periodic_x"]
    if_periodic_synth = lattice_params["if_periodic_y"]
    phys_len = lattice_params["Lx"]
    synth_len = lattice_params["Ly"]

    struc_fact = 0.0
    for j in 1:phys_len
        for s in 1:synth_len
            p1 = (j,s)
            #p1_linear = linear_index(p1,phys_len,synth_len)
            for jj in 1:phys_len
                for ss in 1:synth_len
                    p2 = (jj,ss)
                    #p2_linear = linear_index(p2,phys_len,synth_len)
                    dist = find_dist(p1, p2, (phys_len,synth_len), (if_periodic_phys,if_periodic_synth))[2]
                    struc_fact += occs[p1[2],p1[1]] * occs[p2[2],p2[1]] * exp(im * dot(qvec,dist))
                end
            end
        end
    end
    return struc_fact / sum(occs)
end

# needs testing
function cdw_sf(rho::Array{ComplexF64,2},psi::Vector{ComplexF64},lattice_params::Dict,qend=(0.0,3.0),howmany=50; kwargs...)
	if_plot = get(kwargs, :if_plot, true)

    if qend[1] == 0.0
        direction = 2
        axis_dir = "Y"
        qmax = qend[2]
    elseif qend[2] == 0.0
        direction = 1
        axis_dir = "X"
        qmax = qend[1]
    else
        direction = 3
    end

    if direction < 3
        qvec = zeros(Float64,2)
        qs = range(0.0,stop=qmax,length=howmany)
        struct_factor = zeros(ComplexF64,howmany)
        for (i,qx) in enumerate(qs)
            qvec[direction] = qx
            struct_factor[i] = cdw_structure_factor(rho,qvec,psi,lattice_params; kwargs...)
        end

        if if_plot
            plot_title = get(kwargs, :plot_title, "")
            plot_label = get(kwargs, :plot_label, nothing)
            isnothing(plot_label) ? fig = figure() : nothing
            isnothing(plot_label) ? plot(qs,abs.(struct_factor)) : plot(qs,abs.(struct_factor),label=plot_label)
            xlabel("q"*axis_dir)
            ylabel("Structure Factor")
            title("CDW Structure Factor, " * plot_title)
            legend()
        end
    else
        struct_factor = zeros(ComplexF64,howmany,howmany)
        qxs = range(-qend[1],stop=qend[1],length=howmany)
        qys = range(-qend[2],stop=qend[2],length=howmany)
        for (i,qx) in enumerate(qxs)
            for (j,qy) in enumerate(qys)
                qvec = [qx,qy]
                struct_factor[i,j] = cdw_structure_factor(rho,qvec,psi,lattice_params; kwargs...)
            end
        end

        if if_plot
            plot_title = get(kwargs, :plot_title, "")
            fig = figure()
            imshow(abs.(struct_factor),extent=[qxs[1],qxs[end],qys[1],qys[end]])
            colorbar()
            xlabel("qx")
            ylabel("qy")
            title("CDW Structure Factor, " * plot_title)
        end
    end

	return struct_factor
end


# Spin Stiffness as a function of twist angles theta
# needs testing
function spin_stiffness(energies::Vector,thetas::Vector; kwargs...)
    if_plot = get(kwargs, :if_plot, true)

    stiffnesses = zeros(Float64,length(thetas)-2)
    for i in 2:length(thetas)-1
        stiffnesses[i-1] = (energies[i+1] - 2*energies[i] + energies[i-1]) / (thetas[i+1] - thetas[i])^2
    end

    if if_plot
        fig = figure()
        plot(thetas[2:end-1],stiffnesses,"-p")
        xlabel("Theta")
        ylabel("Spin Stiffness")
        title("Spin Stiffness")
    end

    return stiffnesses
end

# needs testing, not sure what this is
function charge_polarization(psi::Vector{ComplexF64},lattice_params::Dict; kwargs...)
    occs = get_occupancy(psi,lattice_params; if_plot=false)
    Ly = lattice_params["Ly"]

    cppul = sum([sum(occs[m,:] .* m) for m in 1:Ly])/(Ly) # charge polarization per unit length
    return cppul
end

# gets bulk density from wavefunction given size of edge region
function bulk_density(psi::Vector{ComplexF64},lattice_params::Dict,bulk_width_phys::Int64=1,bulk_width_synth::Int64=1; kwargs...)
    occ_mat::Matrix = get_occupancy(psi,lattice_params; if_plot=false)
    #size(occ_mat)[1] == size(occ_mat)[2] ? bulk_width_virt = bulk_width_phys : nothing
	bulk_occ_mat::Matrix = occ_mat[1+bulk_width_phys:end-bulk_width_phys,1+bulk_width_virt:end-bulk_width_virt]
	bulk_density::Float64 = sum(bulk_occ_mat)/prod(size(bulk_occ_mat))
	return bulk_density
end


# functions for two- and four-point correlators
# needs testing
function twopointcorrelator(densmat::Array{ComplexF64,2},lattice_params::Dict; kwargs...)
    !isodd(lattice_params["Lx"]) && !isodd(lattice_params["Ly"]) ? error("Must have odd sized lattice") : nothing
    lattice_params["Lx"] != lattice_params["Ly"] ? error("Not implemented for non-square lattice") : nothing

    if_plot = get(kwargs,:if_plot,true)

    center_site = (Int((lattice_params["Lx"]+1)/2),Int((lattice_params["Ly"]+1)/2))
    center_site_linear = linear_index(center_site,lattice_params["Lx"],lattice_params["Ly"])

    twopointcorrs = zeros(Float64,lattice_params["Lx"],lattice_params["Ly"])
    for j in 1:lattice_params["Lx"]
        for s in 1:lattice_params["Ly"]
            next_site = (j,s)
            next_site_linear = linear_index(next_site,lattice_params["Lx"],lattice_params["Ly"])
            twopointcorrs[next_site[1],next_site[2]] = abs2(densmat[center_site_linear,next_site_linear])
        end
    end
    twopointcorrs ./= twopointcorrs[center_site[1],center_site[2]]

    if_plot ? plot_twopointcorrelator(twopointcorrs; kwargs...) : nothing

    return twopointcorrs
end

function plot_twopointcorrelator(twopointcorrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    imshow(twopointcorrs)
    colorbar()
    title("Two Point Correlator"*plot_title)
    xlabel("Physical")
    ylabel("Synthetic")
    return nothing
end

function pointfourhopping(psi::Vector{ComplexF64},s1::Int64,s2::Int64,lattice_params::Dict)
    hopmat = buildHopping(lattice_params,s1,s2)

    if s1 != s2
        return conj(transpose(psi)) * (conj(transpose(hopmat)) * hopmat) * psi
    else
        return conj(transpose(psi)) * (conj(transpose(hopmat)) * hopmat - hopmat) * psi
    end
end

function fourpointcorrelator(psi::Vector{ComplexF64},lattice_params::Dict; kwargs...)
    !isodd(lattice_params["Lx"]) && !isodd(lattice_params["Ly"]) ? error("Must have odd sized lattice") : nothing
    lattice_params["Lx"] != lattice_params["Ly"] ? error("Not implemented for non-square lattice") : nothing

    if_plot = get(kwargs,:if_plot,true)

    center_site = (Int((lattice_params["Lx"]+1)/2),Int((lattice_params["Ly"]+1)/2))
    center_site_linear = linear_index(center_site,lattice_params["Lx"],lattice_params["Ly"])

    fourpointcorrs = zeros(Float64,lattice_params["Lx"],lattice_params["Ly"])
    for j in 1:lattice_params["Lx"]
        for s in 1:lattice_params["Ly"]
            next_site = (j,s)
            next_site_linear = linear_index(next_site,lattice_params["Lx"],lattice_params["Ly"])
            fourpointcorrs[next_site[1],next_site[2]] = abs2(pointfourhopping(psi,center_site_linear,next_site_linear,lattice_params))
        end
    end

    #fourpointcorrs ./= fourpointcorrs[center_site[1],center_site[2]]

    if_plot ? plot_fourpointcorrelator(fourpointcorrs; kwargs...) : nothing

    return fourpointcorrs
end

function plot_fourpointcorrelator(fourpointcorrs::Array{Float64,2}; kwargs...)
    plot_title = get(kwargs,:plot_title,"")
    fig = figure()
    imshow(fourpointcorrs)
    colorbar()
    title("Four Point Correlator" * plot_title)
    xlabel("Physical")
    ylabel("Synthetic")
    return nothing
end
























"fin"

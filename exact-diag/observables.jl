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
    fix_colorbar = get(kwargs,:fix_colorbar,true)
	fig = figure()
	fix_colorbar ? imshow(exp_occ;vmin=0,vmax=maximum(exp_occ)) : imshow(exp_occ)
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

function density_operator(lattice_params::Dict,which_site::Tuple{Int64,Int64})
    full_basis = lattice_params["full_basis"]
    linear_site = linear_index(which_site,lattice_params["Lx"],lattice_params["Ly"])

    nhat = spzeros(Int64,size(full_basis)[2],size(full_basis)[2])

    for i in 1:size(full_basis)[2]
        linear_site in full_basis[:,i] ? nhat[i,i] += 1 : nothing
    end

    return nhat
end

# make two site density correlations
function make_density_correlations(wavefunc::Vector{ComplexF64},lattice_params::Dict; kwargs...)
    if_zero_shift::Bool = get(kwargs,:if_zero_shift,true)
    if_save::Bool = get(kwargs,:if_save,false)
    Lphys::Int64,Lsynth::Int64 = lattice_params["Lx"],lattice_params["Ly"]

    density_correlations::Array{Float64,4} = zeros(Lphys,Lphys,Lsynth,Lsynth)

    for s in 1:Lsynth
        for j in 1:Lphys
            nhat_first = density_operator(lattice_params,(j,s))
            for ss in 1:Lsynth
                for jj in 1:Lphys
                    nhat_prime = density_operator(lattice_params,(jj,ss))
                    corr_val::Float64 = real(conj(transpose(wavefunc)) * (nhat_first * (nhat_prime * wavefunc)))
                    if if_zero_shift
                        corr_val -= real(conj(transpose(wavefunc)) * nhat_first * wavefunc) * real(conj(transpose(wavefunc)) * nhat_prime * wavefunc)
                    end
                    density_correlations[j,jj,s,ss] = corr_val
                end
            end
        end
    end

    if if_save
        filepath = kwargs[:filepath]
        save_dd_correlation(density_correlations,filepath)
    end

    return density_correlations
end

# finds the distance correlation of density-density along physical direction
function get_dd_distance_correlations_physical(dens_corr_mat::Array{Float64}; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)
    Lphys::Int64 = size(dens_corr_mat)[1]
    Lsynth::Int64 = size(dens_corr_mat)[3]

    distances::Vector{Int64} = collect(1:Lphys-1)
    distance_correlations::Vector{Float64} = zeros(Float64,Lphys-1)

    for j in 1:Lphys
        for d in 1:Lphys-1
            next_site = mod1(j+d,Lphys)
            for s in 1:Lsynth
                distance_correlations[d] += dens_corr_mat[j,next_site,s,s] / Lsynth
            end
        end
    end

    if_plot ? plot_distdensdenscorrs(distances,distance_correlations,"Physical"; kwargs...) : nothing

    return distances,distance_correlations
end

# finds the distance correlation of density-density along synthetic direction
function get_dd_distance_correlations_synthetic(dens_corr_mat::Array{Float64}; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)
    Lphys::Int64 = size(dens_corr_mat)[1]
    Lsynth::Int64 = size(dens_corr_mat)[3]

    distances::Vector{Int64} = collect(1:Lsynth-1)
    distance_correlations::Vector{Float64} = zeros(Float64,Lsynth-1)

    for s in 1:Lsynth
        for d in 1:Lsynth-1
            next_site = mod1(s+d,Lsynth)
            for j in 1:Lphys
                distance_correlations[d] += dens_corr_mat[j,j,s,next_site] / Lphys
            end
        end
    end

    if_plot ? plot_distdensdenscorrs(distances,distance_correlations,"Synthetic"; kwargs...) : nothing

    return distances,distance_correlations
end

# finds the distance correlation of density-density for given direction
function get_distDDcorrs(dens_corr_mat::Array{Float64},direction::String; kwargs...)
    if direction == "x"
        return get_dd_distance_correlations_physical(dens_corr_mat; kwargs...)
    elseif direction == "y"
        return get_dd_distance_correlations_synthetic(dens_corr_mat; kwargs...)
    elseif direction == "both"
        dists1,corrs1 = get_dd_distance_correlations_physical(dens_corr_mat; kwargs...)
        dists2,corrs2 = get_dd_distance_correlations_synthetic(dens_corr_mat; kwargs...)
        return dists1,corrs1,dists2,corrs2
    else
        error("Invalid direction")
    end
end

function ft_densitydensity_correlation(momentum_angle::Float64,momentum_radius::Float64,wavefunc::Union{Nothing,Vector{ComplexF64}},lattice_params::Union{Dict,Nothing}; kwargs...)
    denscorrs = get(kwargs,:denscorrs,nothing)
    if isnothing(denscorrs)
        denscorrs = make_density_correlations(wavefunc,lattice_params; kwargs...)
        lx,ly = lattice_params["Lx"],lattice_params["Ly"]
    else
        ly,lx = size(denscorrs)[1],size(denscorrs)[3]
    end
    if_save::Bool = get(kwargs,:if_save,false)

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

    result::ComplexF64 = sum(denscorrs .* exp.(im .* all_distances)) / ((lx * ly)^2)

    if if_save
        filepath = kwargs[:filepath]
        save_ft_dd(result,round(momentum_angle/pi,digits=3),filepath)
    end

    return result
end

function ft_densitydensity_correlation(momentum::Vector{Float64},wavefunc::Union{Nothing,Vector{ComplexF64}},lattice_params::Union{Dict,Nothing}; kwargs...)
    if momentum == [0.0,0.0]
        return ft_densitydensity_correlation(0.0,0.0,wavefunc,lattice_params; kwargs...)
    else
        return ft_densitydensity_correlation(atan(momentum[2]/momentum[1]),sqrt(momentum[1]^2 + momentum[2]^2),wavefunc,lattice_params; kwargs...)
    end
end

function ft_density(momentum::Vector{Float64},occs::Matrix{Float64}; kwargs...)
    ly,lx = size(occs)

    all_positions::Array{Float64,2} = zeros(Float64,ly,lx)
    for j in 1:lx
        for s in 1:ly
            all_positions[s,j] = dot(momentum,[j,s])
        end
    end

    result::ComplexF64 = sum(occs .* exp.(im .* all_positions)) / (lx * ly)

    return result
end

function ft_density(momentum::Vector{Float64},wavefunc::Vector{ComplexF64},lattice_params::Dict; kwargs...)
    occs = get_occupancy(wavefunc,lattice_params; if_plot=false)
    return ft_density(momentum,occs; kwargs...)
end


function findall_ft_dd(lx::Int64,ly::Int64,n::Int64; kwargs...)
    hanis::Float64 = get(kwargs,:hanis,1.0)
    if_plot::Bool = get(kwargs,:if_plot,false)

    ks = range(0.0,2*pi,length=51)

    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis)])
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0)
    
    filter!(x -> !occursin("twist_angle1",x),all_files)
    filter!(x -> !occursin("mk",x),all_files)

    intstrens = Float64[]
    max_results = Float64[]
    location_results = []
    for f in all_files

        filepath = dataloc * "/" * f
        d,m = read_data_jld2(filepath; output_level=0)

        if !haskey(m,"dens_corr_mat")
            continue
        end

        push!(intstrens,m["U"][end])

        println("Working on $(lx)x$(ly) n=$n at Interaction Strength $(m["U"][end])")

        all_ft_vals::Matrix{Float64} = zeros(Float64,length(ks),length(ks))
        for (idx,kx) in enumerate(ks)
            for (idx2,ky) in enumerate(ks)
                all_ft_vals[idx,idx2] = abs(ft_densitydensity_correlation([kx,ky],nothing,nothing; denscorrs=m["dens_corr_mat"]))
            end
        end

        normalization_factor::Float64 = integrate_2d_matrix(all_ft_vals)

        append!(max_results,[maximum(all_ft_vals) / normalization_factor])
        append!(location_results,[argmax(all_ft_vals)])
    end

    if if_plot
        fig = figure()
        scatter(intstrens,max_results)
        xlabel("Interaction Strength")
        ylabel("Fourier Transform Density-Density Correlation")
        title("FT-DD Maximum for $(lx)x$(ly) N=$n")
    end
    
    return intstrens,max_results,location_results
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

# measures the flatness parameter which is max E2 - E1 / E3 - E1
function twist_flatness_ed(lx::Int,ly::Int,n::Int; kwargs...)
    hanis = get(kwargs,:hanis,1.0)
    if_plot = get(kwargs,:if_plot,true)
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    max_intstren::Float64 = get(kwargs,:max_intstren,2.0)
    
    all_flatnesses = Dict()
    tw1s = Dict()
    tw2s = Dict()
    intstrens = Float64[]

    params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis)])
    all_files = find_data_file(params_dict,"ed",dataloc; output_level=0)
    for f in all_files

        filename_dict = get_params_dict_from_filename(f)
        if !haskey(filename_dict,"twist_angle1")
            continue
        end

        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

        intstren = m["U"][end]
        if intstren > max_intstren
            continue
        end

        if !haskey(all_flatnesses,string(intstren))
            all_flatnesses[string(intstren)] = []
            tw1s[string(intstren)] = []
            tw2s[string(intstren)] = []
        end

        #display(f)

        append!(tw1s[string(intstren)],m["twist_angle"][1])
        append!(tw2s[string(intstren)],m["twist_angle"][2])

        flatness = (d["nrg"][2] - d["nrg"][1]) / (d["nrg"][3] - d["nrg"][1])
        append!(all_flatnesses[string(intstren)],[flatness])
    end

    flatnesses = Float64[]
    for (k,v) in all_flatnesses
        filter!(x->x >= 0.0 && x <= 1.0,v)
        append!(flatnesses,[maximum(v)])
        append!(intstrens,parse(Float64,k))
    end

    if_plot ? plot_twistflatness_vs_intstren(intstrens,flatnesses; kwargs...) : nothing

    return intstrens,flatnesses
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

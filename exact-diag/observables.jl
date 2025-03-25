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
            println("Working on site $(j) $(s)")
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
    if_pinning = get(kwargs,:if_pinning,false)
    if_plot = get(kwargs,:if_plot,true)
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    max_intstren::Float64 = get(kwargs,:max_intstren,2.0)
    
    all_flatnesses = Dict()
    tw1s = Dict()
    tw2s = Dict()
    intstrens = Float64[]

    params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis)])
    if if_pinning
        params_dict["if_pinning"] = true
    end
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

    #= this is the old verion, saved data from before 12 Jan 25 uses this
    if s1 != s2
        return conj(transpose(psi)) * (conj(transpose(hopmat)) * hopmat) * psi
    else
        return conj(transpose(psi)) * (conj(transpose(hopmat)) * hopmat - hopmat) * psi
    end=#

    return conj(transpose(psi)) * (conj(transpose(hopmat)) * hopmat) * psi
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

function fourpoint_alberto(psi::Vector{ComplexF64},lattice_params::Dict; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)
    
    # only need single site occs, would be good to make function for this instead of global occs
    occs::Union{Nothing,Matrix{Float64}} = get(kwargs,:occs,nothing)
    if isnothing(occs)
        occs = get_occupancy(psi,lattice_params; if_plot=false)
    end

    center_site::Vector{Int64} = [5,3]#[Int64(ceil(lattice_params["Lx"]/2)),Int64(ceil(lattice_params["Ly"]/2))]
    center_linear::Int = linear_index(center_site,lattice_params["Lx"],lattice_params["Ly"])

    println("Center site is $(center_site) with linear index $(center_linear)")

    rez::Matrix{Float64} = zeros(Float64,lattice_params["Ly"],lattice_params["Lx"])
    for j in 1:lattice_params["Lx"]
        for s in 1:lattice_params["Ly"]
            println("Working on site ($(j), $(s))")
            rez[s,j] = pointfourhopping(psi,center_linear,linear_index((j,s),lattice_params["Lx"],lattice_params["Ly"]),lattice_params)
        end
    end
    
    # shift from the commutation relation
    rez .-= occs[center_site[2],center_site[1]]
    rez .*= -1

    if_plot ? plot_fourpointcorrelator(rez; kwargs...) : nothing

    return rez
end

function pairdistribution(psi::Vector{ComplexF64},lattice_params::Dict; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)
    occs::Union{Nothing,Matrix{Float64}} = get(kwargs,:occs,nothing)
    if isnothing(occs)
        occs = get_occupancy(psi,lattice_params; kwargs...,if_plot=false)
    end

    Lsynth::Int64,Lphys::Int64 = size(occs)

    fourpoint::Matrix{Float64} = fourpoint_alberto(psi,lattice_params; kwargs...)
    centersite::Vector{Int64} = [Int64(ceil(Lphys/2)),Int64(ceil(Lsynth/2))]
    pairdist::Matrix{Float64} = fourpoint ./ (occs[centersite[2],centersite[1]] .* occs)

    if_plot ? plot_pairdistribution(pairdist; kwargs...) : nothing

    return pairdist
end

# these are not exactly correct I think
# calculate canonical momentum -i * d H(j,s) / d flux_density
function canonical_momentum(psi::Vector{ComplexF64},direction::Int64,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)

    Lx::Int64 = lattice_params["Lx"]
    Ly::Int64 = lattice_params["Ly"]
    flux_density::Float64 = hamilt_params["alpha"][direction]

    canonical_momentum::Matrix{Float64} = zeros(Float64,Lx,Ly)
    for j in 1:Lx
        for s in 1:Ly
            starting_site::Tuple = (j,s)
            starting_site_linear::Int64 = linear_index(starting_site,Lx,Ly)

            change_vector = zeros(Int64,2)
            change_vector[direction] = 1
            ending_site = starting_site .+ change_vector
            ending_site[1] = mod1(ending_site[1],Lx)
            ending_site[2] = mod1(ending_site[2],Ly)
            ending_site_linear::Int64 = linear_index(ending_site,Lx,Ly)

            hopmat = buildHopping(lattice_params,starting_site_linear,ending_site_linear)
            exp_val = conj(transpose(psi)) * (hopmat .- conj(transpose(hopmat))) * psi
            canonical_momentum[j,s] = real(2*pi*im*flux_density*exp_val)
        end
    end

    if_plot ? plot_canonical_momentum(canonical_momentum; kwargs...) : nothing

    return canonical_momentum
end

function average_canonical_momentum(psi::Vector{ComplexF64},lattice_params::Dict,hamilt_params::Dict; kwargs...)
    canonical_momentum_x = canonical_momentum(psi,1,lattice_params,hamilt_params; kwargs...)
    canonical_momentum_y = canonical_momentum(psi,2,lattice_params,hamilt_params; kwargs...)

    return sum(canonical_momentum_x) / prod(size(canonical_momentum_x)),sum(canonical_momentum_y) / prod(size(canonical_momentum_y))
end


function find_subsystem_B(subsystem_A::Vector{Int64}, lx::Int64, ly::Int64)
    subsystem_B = []
    for i in 1:lx*ly
        if !(i in subsystem_A)
            push!(subsystem_B,i)
        end
    end
    return subsystem_B
end

# Function to partition a configuration into subsystems A and B
function partition_config(config, subsystem_A, subsystem_B)
    config_A = [pos for pos in config if pos in subsystem_A]
    config_B = [pos for pos in config if pos in subsystem_B]
    return config_A, config_B
end

# Group coefficients by configurations in A and B
function group_by_subsystems(wavefunc, full_basis, subsystem_A, subsystem_B)
    grouped_coeffs = Dict{Tuple{Vector{Int}, Vector{Int}}, ComplexF64}()
    for i in 1:size(full_basis,2)
        config = full_basis[:,i]
        config_A, config_B = partition_config(config, subsystem_A, subsystem_B)
        key = (config_A, config_B)
        grouped_coeffs[key] = get(grouped_coeffs, key, 0.0) + wavefunc[i]
    end
    return grouped_coeffs
end

# Compute the reduced density matrix for subsystem A
function compute_reduced_density_matrix(wavefunc, full_basis, subsystem_A, subsystem_B)
    grouped_coeffs = group_by_subsystems(wavefunc, full_basis, subsystem_A, subsystem_B)
    unique_A_configs = unique(key[1] for key in keys(grouped_coeffs))
    dim_A = length(unique_A_configs)
    rho_A = zeros(ComplexF64, dim_A, dim_A)

    # Map each unique_A_config to an index
    config_A_to_index = Dict(config => i for (i, config) in enumerate(unique_A_configs))

    # Compute rho_A matrix elements
    for ((config_A1, config_B1), coeff1) in grouped_coeffs
        #println(round(100*idx/length(grouped_coeffs),digits=4)," % done")
        for ((config_A2, config_B2), coeff2) in grouped_coeffs
            if config_B1 == config_B2
                i = config_A_to_index[config_A1]
                j = config_A_to_index[config_A2]
                rho_A[i, j] += coeff1 * conj(coeff2)
            end
        end
    end

    return rho_A, unique_A_configs
end

function entanglement_entropy(rho_A::Matrix{ComplexF64}; tolerance=1e-8)
    # Diagonalize rho_A to get eigenvalues
    eigenvalues = real.(eigen(rho_A).values)

    # Filter out very small eigenvalues (numerical precision issues)
    eigenvalues = eigenvalues[eigenvalues .> tolerance]

    entanglement_entropy = -sum(eigenvalues .* log.(eigenvalues))

    return entanglement_entropy
end


function four_point_operator_old(site1::Int64,site2::Int64,site3::Int64,site4::Int64,lattice_params::Dict)
    
    mat1 = buildHopping(lattice_params,site3,site1)
    mat2 = buildHopping(lattice_params,site4,site2)

    result = mat1 * mat2

    site2 == site3 && (result -= buildHopping(lattice_params,site4,site1))

    return result
end

function four_point_operator_singlebasis(site1::Int64,site2::Int64,site3::Int64,site4::Int64,this_basis_state)
    coeff::Float64 = 1.0
    N = length(this_basis_state)
    #println("Working on basis state $this_basis_state")
    #println("Hopping from sites $site3 and $site4 to sites $site1 and $site2 using basis state $this_basis_state")
    
    output_state = zeros(Int64,length(this_basis_state)) + this_basis_state
    
    n4 = length(filter(x->x==site4,output_state))
    deleteat!(output_state,findfirst(x->x==site4,output_state))
    coeff *= sqrt(n4)

    n3 = length(filter(x->x==site3,output_state))
    deleteat!(output_state,findfirst(x->x==site3,output_state))
    coeff *= sqrt(n3)

    n2 = length(filter(x->x==site2,output_state))
    push!(output_state,site2)
    coeff *= sqrt(n2+1)
    
    n1 = length(filter(x->x==site1,output_state))
    push!(output_state,site1)
    coeff *= sqrt(n1+1)

    sort!(output_state,rev=true)
    #println("Output state is $output_state with coeff $coeff")
    
    if length(unique(output_state)) < N
        #println("Error: output state $output_state has double occupancy")
        return nothing,0.0
    end
    if length(output_state) > N
        error("Too long output state: $output_state")
        #continue
    end

    return output_state,coeff
end

function four_point_operator(site1::Int64,site2::Int64,site3::Int64,site4::Int64,lattice_params::Dict)
    full_basis = lattice_params["full_basis"]
    N = lattice_params["N"]

    dubhop = spzeros(Int64,size(full_basis)[2],size(full_basis)[2])

    site3 == site4 && (return dubhop)
    site1 == site2 && (return dubhop)

    for j in 1:size(full_basis)[2]
        this_basis_state = full_basis[:,j]
        #println("Working on basis state $this_basis_state")

        if site3 in this_basis_state && site4 in this_basis_state
            coeff::Float64 = 1.0
            #println("Working on basis state $this_basis_state")
            #println("Hopping from sites $site3 and $site4 to sites $site1 and $site2 using basis state $this_basis_state")
            
            output_state = zeros(Int64,length(this_basis_state)) + this_basis_state
            
            n4 = length(filter(x->x==site4,output_state))
            deleteat!(output_state,findfirst(x->x==site4,output_state))
            coeff *= sqrt(n4)

            n3 = length(filter(x->x==site3,output_state))
            deleteat!(output_state,findfirst(x->x==site3,output_state))
            coeff *= sqrt(n3)

            n2 = length(filter(x->x==site2,output_state))
            push!(output_state,site2)
            coeff *= sqrt(n2+1)
            
            n1 = length(filter(x->x==site1,output_state))
            push!(output_state,site1)
            coeff *= sqrt(n1+1)

            sort!(output_state,rev=true)
            #println("Output state is $output_state with coeff $coeff")
            
            if length(unique(output_state)) < N
                #println("Error: output state $output_state has double occupancy")
                continue
            end
            if length(output_state) > N
                error("Too long output state: $output_state")
                #continue
            end
            bind = find_basis_index(output_state)
            dubhop[j,bind] += coeff

        end
    end

    return dubhop
end

function ft_fourpt(psi::Vector{ComplexF64},momentum1::Vector{Float64},momentum2::Vector{Float64},lattice_params::Dict; kwargs...)
    
    lx::Int64 = lattice_params["Lx"]
    ly::Int64 = lattice_params["Ly"]

    fourpt::ComplexF64 = 0.0
    for s1 in 1:lx*ly
        coord1 = coordinate(s1,lx,ly)
        coeff1::ComplexF64 = ft_coeff(coord1,momentum1,"Adag")
        for s2 in 1:lx*ly
            coord2 = coordinate(s2,lx,ly)
            coeff2::ComplexF64 = ft_coeff(coord2,momentum2,"Adag")
            println("Working on s1=$(s1) s2=$(s2)")
            for s3 in 1:lx*ly
                coord3 = coordinate(s3,lx,ly)
                coeff3::ComplexF64 = ft_coeff(coord3,momentum2,"A")
                for s4 in 1:lx*ly
                    coord4 = coordinate(s4,lx,ly)
                    coeff4::ComplexF64 = ft_coeff(coord4,momentum1,"A")

                    coeff::ComplexF64 = coeff1 * coeff2 * coeff3 * coeff4

                    big_operator = four_point_operator(s1,s2,s3,s4,lattice_params)
                    
                    fourpt += coeff * (conj(transpose(psi)) * (big_operator * psi))
                end
            end
        end
    end

    return fourpt / (lx*ly)^2
end

function ft_fourpt(wavefuncs::Vector{Vector{ComplexF64}},momentum1::Vector{Float64},momentum2::Vector{Float64},lattice_params::Dict; kwargs...)
    lx::Int64 = lattice_params["Lx"]
    ly::Int64 = lattice_params["Ly"]

    #=fourpt::Matrix{ComplexF64} = zeros(ComplexF64,length(wavefuncs),length(wavefuncs))
    for s1 in 1:lx*ly
        coord1 = coordinate(s1,lx,ly)
        coeff1::ComplexF64 = ft_coeff(coord1,momentum1,"Adag")
        for s2 in 1:lx*ly
            coord2 = coordinate(s2,lx,ly)
            coeff2::ComplexF64 = ft_coeff(coord2,momentum2,"Adag")
            println("Working on s1=$(s1) s2=$(s2)")
            for s3 in 1:lx*ly
                coord3 = coordinate(s3,lx,ly)
                coeff3::ComplexF64 = ft_coeff(coord3,momentum2,"A")
                for s4 in 1:lx*ly
                    coord4 = coordinate(s4,lx,ly)
                    coeff4::ComplexF64 = ft_coeff(coord4,momentum1,"A")

                    coeff::ComplexF64 = coeff1 * coeff2 * coeff3 * coeff4

                    big_operator = four_point_operator(s1,s2,s3,s4,lattice_params)
                    
                    for i in 1:length(wavefuncs)
                        for j in 1:length(wavefuncs)
                            fourpt[i,j] += coeff * (conj(transpose(wavefuncs[j])) * (big_operator * wavefuncs[i]))
                        end
                    end

                end
            end
        end
    end=#

    fourpt::Matrix{ComplexF64} = zeros(ComplexF64,length(wavefuncs),length(wavefuncs))
    for s1 in 1:lx*ly
        coord1 = coordinate(s1,lx,ly)
        coeff1::ComplexF64 = ft_coeff(coord1,momentum1,"Adag")
        for s3 in 1:lx*ly
            coord3 = coordinate(s3,lx,ly)
            coeff3::ComplexF64 = ft_coeff(coord3,momentum2,"Adag")
            mat1 = buildHopping(lattice_params,s1,s3)
            println("Working on s1=$(s1) s3=$(s3)")
            for s4 in 1:lx*ly
                coord4 = coordinate(s4,lx,ly)
                coeff4::ComplexF64 = ft_coeff(coord4,momentum2,"A")
                mat3 = buildHopping(lattice_params,s1,s4)
                for s2 in 1:lx*ly
                    coord2 = coordinate(s2,lx,ly)
                    coeff2::ComplexF64 = ft_coeff(coord2,momentum1,"A")
                    mat2 = buildHopping(lattice_params,s2,s3)

                    big_operator = mat1 * mat2
                    s2 == s3 && (big_operator -= mat3)

                    for i in 1:length(wavefuncs)
                        for j in 1:length(wavefuncs)
                            fourpt[i,j] += coeff1 * coeff2 * coeff3 * coeff4 * (conj(transpose(wavefuncs[i])) * (big_operator * wavefuncs[j]))
                        end
                    end
                end
            end
        end
    end

    display(fourpt)

    return eigvals(fourpt ./ (lx*ly)^2)
end

function two_point_operator(site1::Int64,site2::Int64,lattice_params::Dict)
    mat1 = buildHopping(lattice_params,site2,site1)
    return mat1
end

function ft_twopt(psi::Vector{ComplexF64},momentum1::Vector{Float64},momentum2::Vector{Float64},lattice_params::Dict; kwargs...)
    
    lx::Int64 = lattice_params["Lx"]
    ly::Int64 = lattice_params["Ly"]

    twopt::ComplexF64 = 0.0
    for s1 in 1:lx*ly
        coord1 = coordinate(s1,lx,ly)
        coeff1::ComplexF64 = ft_coeff(coord1,momentum1,"Adag")
        for s2 in 1:lx*ly
            coord2 = coordinate(s2,lx,ly)
            coeff2::ComplexF64 = ft_coeff(coord2,momentum2,"A")
            println("Working on s1=$(s1) s2=$(s2)")
            big_operator = two_point_operator(s1,s2,lattice_params)

            twopt += coeff1 * coeff2 * (conj(transpose(psi)) * (big_operator * psi))
        end
    end

    return twopt / (lx*ly)
end

function ft_twopt(wavefuncs::Vector{Vector{ComplexF64}},momentum1::Vector{Float64},momentum2::Vector{Float64},lattice_params::Dict; kwargs...)
    
    lx::Int64 = lattice_params["Lx"]
    ly::Int64 = lattice_params["Ly"]

    twopt::Matrix{ComplexF64} = zeros(ComplexF64,length(wavefuncs),length(wavefuncs))
    for s1 in 1:lx*ly
        coord1 = coordinate(s1,lx,ly)
        coeff1::ComplexF64 = ft_coeff(coord1,momentum1,"Adag")
        for s2 in 1:lx*ly
            coord2 = coordinate(s2,lx,ly)
            coeff2::ComplexF64 = ft_coeff(coord2,momentum2,"A")
            #println("Working on s1=$(s1) s2=$(s2)")
            big_operator = two_point_operator(s1,s2,lattice_params)

            for i in 1:length(wavefuncs)
                for j in 1:length(wavefuncs)
                    twopt[i,j] += coeff1 * coeff2 * (conj(transpose(wavefuncs[j])) * (big_operator * wavefuncs[i]))
                end
            end

        end
    end

    #display(twopt)

    return eigvals(twopt ./ (lx*ly))
end


function ft_fourpt_numop(wavefunc::Vector{ComplexF64},momentum::Vector{Float64},lattice_params::Dict; kwargs...)
    
    lx::Int64 = lattice_params["Lx"]
    ly::Int64 = lattice_params["Ly"]
    n::Int64 = lattice_params["N"]

    term1::ComplexF64 = 0.0
    term2::ComplexF64 = 0.0
    for s1 in 1:lx*ly
        coord1 = coordinate(s1,lx,ly)
        coeff1::ComplexF64 = ft_coeff(coord1,momentum,"A")
        mat1 = density_operator(lattice_params,coord1)
        term2 += coeff1 * (conj(transpose(wavefunc)) * (mat1 * wavefunc))
        for s2 in 1:lx*ly
            coord2 = coordinate(s2,lx,ly)
            coeff2::ComplexF64 = ft_coeff(coord2,momentum,"A")
            mat2 = density_operator(lattice_params,coord2)
            term1 += (coeff1 * coeff2) * (conj(transpose(wavefunc)) * ((mat1 * mat2) * wavefunc))
        end
    end

    return term1 * (1 / (lx*ly)),term2 * (1 / sqrt(lx*ly))
end


function ft_fourpt_alberto(psi::Vector{ComplexF64},momentum1::Vector{Float64},momentum2::Vector{Float64},lattice_params::Dict; kwargs...)
    Lx::Int64 = lattice_params["Lx"]
    Ly::Int64 = lattice_params["Ly"]

    m = Int(momentum1[2] * Ly)
    mp = Int(momentum2[2] * Ly)
    if mp > Lx || m > Lx
        error("Momentum out of bounds: m=$m mp=$mp Lx=$Lx")
    elseif mp == 0 || m == 0
        error("Momentum cannot be zero: m=$m mp=$mp")
    end
    alpha = 1/Ly

    fourpt::ComplexF64 = 0.0
    for y1 in 1:Ly
        coord1 = [m,y1]
        s1 = linear_index(coord1,Lx,Ly)
        coeff1::ComplexF64 = ft_coeff_alberto(coord1,momentum1,"Adag",Lx,Ly,m,alpha)
        for y2 in 1:Ly
            coord2 = [mp,y2]
            s2 = linear_index(coord2,Lx,Ly)
            coeff2::ComplexF64 = ft_coeff_alberto(coord2,momentum2,"Adag",Lx,Ly,mp,alpha)
            #println("Working on y1=$(y1) y2=$(y2)")
            for y3 in 1:Ly
                coord3 = [mp,y3]
                s3 = linear_index(coord3,Lx,Ly)
                coeff3::ComplexF64 = ft_coeff_alberto(coord3,momentum2,"A",Lx,Ly,mp,alpha)
                for y4 in 1:Ly
                    coord4 = [m,y4]
                    s4 = linear_index(coord4,Lx,Ly)
                    coeff4::ComplexF64 = ft_coeff_alberto(coord4,momentum1,"A",Lx,Ly,m,alpha)

                    coeff::ComplexF64 = coeff1 * coeff2 * coeff3 * coeff4

                    #println("Working on s1=$(s1) s2=$(s2) s3=$(s3) s4=$(s4)")

                    big_operator = four_point_operator(s1,s2,s3,s4,lattice_params)
                    
                    expectval = (conj(transpose(psi)) * (big_operator * psi))

                    #println("At $y1 $y2 $y3 $y4 coeff is: $(round(coeff,digits=10))")
                    #println("expectation value = $(round(expectval,digits=10))")

                    fourpt += coeff * expectval
                end
            end
        end
    end

    return fourpt
end

function ft_twopt_alberto(wavefunc::Vector{ComplexF64},momentum1::Vector{Float64},momentum2::Vector{Float64},lattice_params::Dict; kwargs...)
    Lx::Int64 = lattice_params["Lx"]
    Ly::Int64 = lattice_params["Ly"]

    mval = Int(momentum1[2] * Ly)
    mval2 = Int(momentum2[2] * Ly)
    alpha = 1/Ly

    twopt::ComplexF64 = 0.0
            
    for y3 in 1:ly
        coord3 = (mval,y3)
        coeff3::ComplexF64 = ft_coeff_alberto(coord3,momentum1,"Adag",Lx,Ly,mval,alpha)
        for y4 in 1:ly
            coord4 = (mval2,y4)
            println("Working on y1=$(y3) y2=$(y4)")
            coeff4::ComplexF64 = ft_coeff_alberto(coord4,momentum2,"A",Lx,Ly,mval2,alpha)

            coeff::ComplexF64 = coeff3 * coeff4
            
            local_exppart = hopping_probability(wavefunc,Tuple(coord3),Tuple(coord4),lattice_params)
            local_twopt = coeff * local_exppart
            #println("At site $coord3 and $coord4, expectation value is $(round(local_exppart,digits=8))")
            twopt += local_twopt
        end
    end

    return twopt
end

function ft_twopt_alberto(wavefuncs::Vector{Vector{ComplexF64}},momentum1::Vector{Float64},momentum2::Vector{Float64},lattice_params::Dict; kwargs...)
    Lx::Int64 = lattice_params["Lx"]
    Ly::Int64 = lattice_params["Ly"]

    mval = Int(momentum1[2] * Ly)
    mval2 = Int(momentum2[2] * Ly)
    alpha = 1/Ly

    twopt::Matrix{ComplexF64} = zeros(ComplexF64,length(wavefuncs),length(wavefuncs))
            
    for y3 in 1:ly
        coord3 = (mval,y3)
        coeff3::ComplexF64 = ft_coeff_alberto(coord3,momentum1,"Adag",Lx,Ly,mval,alpha)
        lin1 = linear_index(coord3,Lx,Ly)
        for y4 in 1:ly
            coord4 = (mval2,y4)
            lin2 = linear_index(coord4,Lx,Ly)
            println("Working on y1=$(y3) y2=$(y4)")
            coeff4::ComplexF64 = ft_coeff_alberto(coord4,momentum2,"A",Lx,Ly,mval2,alpha)

            coeff::ComplexF64 = coeff3 * coeff4
            big_operator = two_point_operator(lin1,lin2,lattice_params)
            for i in 1:length(wavefuncs)
                for j in 1:length(wavefuncs)
                    local_exppart = adjoint(wavefuncs[j]) * (big_operator * wavefuncs[i])
                    local_twopt = coeff * local_exppart
                    #println("At site $coord3 and $coord4, expectation value is $(round(local_exppart,digits=8))")
                    twopt[i,j] += local_twopt
                end
            end
        end
    end

    return eigvals(twopt)
end




















"fin"

#####################################################
#=

This file contains observable calculations for TTNs

Depends on:
    review-practice-codes/ttn.jl

=#
######################################################

include("../other-funcs/include-other-files.jl")

include_other_files(["review-practice-codes/ttn.jl","other-funcs/basic-2d-observables.jl"])

function get_occupancy(ttn::TTN.TreeTensorNetwork; kwargs...)

	exp_occ = abs.(TTN.expect(ttn,"N"))
	
	if_plot = get(kwargs, :if_plot, true)
	
	if_plot && plot_occupancy(transpose(exp_occ); kwargs...)
		
	return transpose(exp_occ)
end

function get_occupancy(densmat::Matrix; kwargs...)
	if isinteger(sqrt(size(densmat)[1]))
		phys_length = Int(sqrt(size(densmat)[1]))
		virt_length = Int(phys_length)
	else
		phys_length = Int(sqrt(2*size(densmat)[1]))
		virt_length = Int(phys_length/2)
	end

	lat = TTN.SimpleLattice((phys_length,virt_length),TTN.ITensorNode,"Boson")

	exp_occ = zeros(phys_length,virt_length)
	for j in 1:phys_length
		for s in 1:virt_length
			linear_index = TTN.linear_ind(lat,(j,s))
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
	write_data(filename,occs_data_dict,location,metadata)
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

function physical_correlation(wavefunc::TTN.TreeTensorNetwork; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
    densmat = get(kwargs,:densmat,nothing)
	if !isnothing(densmat)
        lat = TTN.physical_lattice(TTN.network(wavefunc))
	    Lx,Ly = size(lat)
        return physical_correlation(densmat,Lx,Ly; kwargs...)
    else
        densmat = density_matrix(wavefunc)
    end

	lat = TTN.physical_lattice(TTN.network(wavefunc))
	Lx,Ly = size(lat)

    phys_corrs = Array{Float64,2}(undef,Lx-1,Ly)
    for j in 2:Lx
        for s in 1:Ly
            phys_corr = densmat[TTN.linear_ind(lat,(j,s)),TTN.linear_ind(lat,(1,s))]
            phys_corr /= sqrt(densmat[TTN.linear_ind(lat,(j,s)),TTN.linear_ind(lat,(j,s))] * densmat[TTN.linear_ind(lat,(1,s)),TTN.linear_ind(lat,(1,s))])
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

function synthetic_correlation(wavefunc::TTN.TreeTensorNetwork; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
	densmat = get(kwargs,:densmat,nothing)
	if isnothing(densmat)
		densmat = density_matrix(wavefunc)
    else
        lat = TTN.physical_lattice(TTN.network(wavefunc))
        Lx,Ly = size(lat)
        return synthetic_correlation(densmat,Lx,Ly; kwargs...)
	end

	lat = TTN.physical_lattice(TTN.network(wavefunc))
	Lx,Ly = size(lat)

    syn_corrs = Array{Float64,2}(undef,Lx,Ly-1)
    for j in 1:Lx
        for s in 2:Ly
			syn_corr = densmat[TTN.linear_ind(lat,(j,s)),TTN.linear_ind(lat,(j,1))]
			syn_corr /= sqrt(densmat[TTN.linear_ind(lat,(j,s)),TTN.linear_ind(lat,(j,s))] * densmat[TTN.linear_ind(lat,(j,1)),TTN.linear_ind(lat,(j,1))])
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

function physical_current(wavefunc::TTN.TreeTensorNetwork; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
	densmat = get(kwargs,:densmat,nothing)
	if isnothing(densmat)
		densmat = density_matrix(wavefunc)
    else
        lat = TTN.physical_lattice(TTN.network(wavefunc))
        Lx,Ly = size(lat)
        return physical_current(densmat,Lx,Ly; kwargs...)
	end

    lat = TTN.physical_lattice(TTN.network(wavefunc))
	Lx,Ly = size(lat)

    currents = Array{Float64,2}(undef,Lx,Ly)
    for s in 1:Ly
        for j in 1:Lx
            site1 = TTN.linear_ind(lat,(j,s))
            site2 = TTN.linear_ind(lat,(mod1(j+1,Lx),s))
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
        for s in 1:Ly
            site1 = linear_index((j,s),Lx,Ly)
            site2 = linear_index((j,mod1(s+1,Ly)),Lx,Ly)
            #println("Working on site $(j),$(s) with linear index $(site1) and $(site2)")
            current_val = imag(densmat[site1,site2] - densmat[site2,site1])
            current_normalization = densmat[site1,site1] + densmat[site2,site2]
            current_val /= current_normalization
            currents[j,s] = real(current_val)
        end
    end

    if_plot ? plot_synthetic_current(currents; kwargs...) : nothing

    return currents
end

function synthetic_current(wavefunc::TTN.TreeTensorNetwork; kwargs...)
    if_plot = get(kwargs,:if_plot,true)
	densmat = get(kwargs,:densmat,nothing)
	if isnothing(densmat)
		densmat = density_matrix(wavefunc)
    else
        lat = TTN.physical_lattice(TTN.network(wavefunc))
        Lx,Ly = size(lat)
        return synthetic_current(densmat,Lx,Ly; kwargs...)
	end

	lat = TTN.physical_lattice(TTN.network(wavefunc))
	Lx,Ly = size(lat)

    currents = Array{Float64,2}(undef,Lx,Ly)
    for j in 1:Lx
        for s in 1:Lx
            site1 = TTN.linear_ind(lat,(j,s))
            site2 = TTN.linear_ind(lat,(j,mod1(s+1,Ly)))
            current_val = imag(densmat[site1,site2] - densmat[site2,site1])
            current_normalization = densmat[site1,site1] + densmat[site2,site2]
            current_val /= current_normalization
            currents[j,s] = real(current_val)
        end
    end

    if_plot ? plot_synthetic_current(currents; kwargs...) : nothing

    return currents
end

function check_nrg_convergence(metadata::Dict,if_perfect::Bool=true; kwargs...)
    nrg_tol::Float64 = get(metadata,"nrgtol",5e-5)
    key_type::String = get(kwargs,:key_type,"observer")

    if_perfect ? nothing : nrg_tol *= 10

    all_nrgs::Dict{String,Array{Float64,1}} = Dict()
    all_nrgs["0"] = key_type == "observer" ? metadata[key_type].nrg : metadata[key_type]

    found_excited_nrg::Bool = haskey(metadata,"$(key_type)_1")
    next_nrg_level::Int = 1
    while found_excited_nrg
        all_nrgs[string(next_nrg_level)] = key_type == "observer" ? metadata["$(key_type)_$(next_nrg_level)"].nrg : metadata["$(key_type)_$(next_nrg_level)"]
        next_nrg_level += 1
        found_excited_nrg = haskey(metadata,"$(key_type)_$(next_nrg_level)")
    end

    if_converged::Dict{String,Bool} = Dict()
    for (k,v) in all_nrgs
        if_converged[k] = abs(v[end] - v[end-1]) < nrg_tol
    end

    return all(values(if_converged)),if_converged
end

# find some way to speed this up
function make_density_correlations(wavefunc::TTN.TreeTensorNetwork; kwargs...)
    Lphys::Int64,Lsynth::Int64 = get_lattice_dims(wavefunc)

    density_correlations::Array{Float64,4} = zeros(Lphys,Lphys,Lsynth,Lsynth)

    for s in 1:Lsynth
        for j in 1:Lphys
            for ss in 1:Lsynth
                for jj in 1:Lphys
                    corr_val::Float64 = real(TTN.correlation(wavefunc,"N","N",(j,s),(jj,ss)))
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

function pointfourhopping(psi::TTN.TreeTensorNetwork,s1::Int64,s2::Int64)
    result = TTN.correlation(psi,"N","N",s1,s2)

    if s1 == s2
        result -= TTN.expect(psi,"N",s1)
    end

    return result
end

# to be used when using smaller lattice than TTN actual size
function fourpoint_alberto(psi::TTN.TreeTensorNetwork,restricted_size::Vector{Int64}; kwargs...)
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

function fourpoint_alberto(psi::TTN.TreeTensorNetwork; kwargs...)
    if_restricted_size::Bool = get(kwargs,:if_restricted_size,false)
    if_restricted_size && fourpoint_alberto(psi,if_restricted_size; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,true)

    lx,ly = get_lattice_dims(psi)

    center_site::Vector{Int64} = [5,3]#[Int64(ceil(lx/2)),Int64(ceil(ly/2))]
    println("Center site is $(center_site)")
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

function pairdistribution(psi::TTN.TreeTensorNetwork; kwargs...)
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

function spatial_entanglement_spectrum(psi::TTN.TreeTensorNetwork,layers_down::Int; kwargs...)
    cap = get(kwargs,:if_cap,nothing)

    numlayers = TTN.number_of_layers(psi)
    TTN.move_ortho!(psi,(numlayers - layers_down,1))
    top_tensor = psi[numlayers - layers_down,1]
    idx_left = TTN.inds(top_tensor; tags = "Link,nl=$(numlayers-1-layers_down),np=1")
    u,s,v,spec = svd(top_tensor,idx_left)

    result = isnothing(cap) ? spec.eigs : spec.eigs[1:cap]

    return result
end

function spatial_entanglement_spectrum(psi::TTN.TreeTensorNetwork; kwargs...)
    if_save::Bool = get(kwargs,:if_save,false)
    filepath::Union{Nothing,String} = get(kwargs,:filepath,nothing)
    output_level::Int64 = get(kwargs,:output_level,1)

    layer_count = TTN.number_of_layers(psi)
    entspec = zeros(Float64,TTN.maxlinkdim(psi),layer_count-2)
    
    for ld in 2:layer_count-1
        output_level > 0 && println("Calculating Layer $(ld) Entanglement Spectrum")
        
        spec_vals = spatial_entanglement_spectrum(psi,ld; kwargs...)
        extra_length = size(entspec,1) - length(spec_vals)
        entspec[:,ld-1] = vcat(spec_vals,zeros(extra_length))
    end

    if_save && save_spatial_entspec(entspec,filepath; kwargs...)

    return entspec
end

function save_spatial_entspec(entspec::Matrix{Float64},filepath::String; kwargs...)
    isnothing(filepath) && error("No filepath provided for saving")

    modify_data(Dict([("entanglement_spectrum",entspec)]),filepath,"metadata"; output_level=0)
end

function calculate_perimeter(which_layer::Int64)
    all_perims = [6,8,12,16,24,32,48,62]
    return all_perims[which_layer]
end

function get_layer_from_linkname(linkname::String)
    return parse(Int64,split(split(linkname,",")[2],"=")[end])
end

function calculate_perimeter(linkname::String)
    return calculate_perimeter(get_layer_from_linkname(linkname))
end

function tee(psi::TTN.TreeTensorNetwork,top_layer::Int64; kwargs...)
    if_save::Bool = get(kwargs,:if_save,false)

    numlayers = TTN.number_of_layers(psi)

    if numlayers < top_layer
        error("Top layer is higher than the number of layers in the TTN")
    end

    ee_data::Dict{String,Float64} = Dict()
    cutlink_data::Dict{String,Tuple{Int64,Int64}} = Dict()

    net = TTN.network(psi)

    # start with S_ABDC at the top node
    top_node = (top_layer,1)

    # make sure the ortho_center is at the cutting tensor
    TTN.move_ortho!(psi,top_node)

    # get the tensor, link (and index) to cut (at top always cut the left link thus np=1)
    tensor_abcd = psi[top_node]
    link_abcd = "Link,nl=$(top_layer-1),np=1"
    index_abcd = inds(tensor_abcd; tags = link_abcd)

    # perform svd on the link and sum schmidt values for entanglement entropy
    u,s,v,spec = svd(tensor_abcd,index_abcd)
    s_abcd = entanglement_entropy(spec.eigs)

    # save data to dictionary
    ee_data["s_abcd"] = s_abcd
    cutlink_data["s_abcd"] = (top_layer-1,1)

    # find the child nodes to loop over to find S_AB and S_CD
    middle_node = (top_node[1]-1,1)
    middle_layer_children = TTN.child_nodes(net,middle_node)
    for (i,middle_child) in enumerate(middle_layer_children)
        
        # move ortho_center to the tensor site to cut
        TTN.move_ortho!(psi,middle_node)

        # get tensor, link, and index
        tensor_middle = psi[middle_node]
        link_middle = "Link,nl=$(middle_child[1]),np=$(middle_child[2])"
        index_middle = inds(tensor_middle; tags = link_middle)

        # perform svd on the link and sum schmidt values for entanglement entropy
        u,s,v,spec = svd(tensor_middle,index_middle)
        s_middle = entanglement_entropy(spec.eigs)

        # save data to dictionary
        dict_key = i == 1 ? "s_ab" : "s_cd"
        ee_data[dict_key] = s_middle
        cutlink_data[dict_key] = middle_child

        # now do loop over children of this middle node
        baby_layer_children = TTN.child_nodes(net,middle_child)
        for (j,baby_child) in enumerate(baby_layer_children)
                
                # move ortho_center to the tensor site to cut
                TTN.move_ortho!(psi,middle_child)
    
                # get tensor, link, and index
                tensor_baby = psi[middle_child]
                link_baby = "Link,nl=$(baby_child[1]),np=$(baby_child[2])"
                index_baby = inds(tensor_baby; tags = link_baby)
    
                # perform svd on the link and sum schmidt values for entanglement entropy
                u,s,v,spec = svd(tensor_baby,index_baby)
                s_baby = entanglement_entropy(spec.eigs)
    
                # save data to dictionary
                if i == 1
                    dict_key = j == 1 ? "s_a" : "s_b"
                else
                    dict_key = j == 1 ? "s_c" : "s_d"
                end
                ee_data[dict_key] = s_baby
                cutlink_data[dict_key] = baby_child
        end
    end

    # now for the mixed middle section S_BC and S_AD
    TTN.move_ortho!(psi,middle_node)
    combined_tensor = psi[middle_node] * psi[middle_layer_children[1]] * psi[middle_layer_children[2]]

    # build indices for S_BC
    link_data_b,link_data_c = (cutlink_data["s_b"],cutlink_data["s_c"])
    link_bc = ["Link,nl=$(link_data_b[1]),np=$(link_data_b[2])","Link,nl=$(link_data_c[1]),np=$(link_data_c[2])"]
    index_bc = [inds(combined_tensor; tags = link_bc[1]),inds(combined_tensor; tags = link_bc[2])]

    # perform svd on the link and sum schmidt values for entanglement entropy
    u,s,v,spec = svd(combined_tensor,index_bc)
    s_bc = entanglement_entropy(spec.eigs)

    # save data to dictionary (don't need to save cutlink data because it comes from the S_B and S_C)
    ee_data["s_bc"] = s_bc

    # build indices for S_AD
    link_data_a,link_data_d = (cutlink_data["s_a"],cutlink_data["s_d"])
    link_ad = ["Link,nl=$(link_data_a[1]),np=$(link_data_a[2])","Link,nl=$(link_data_d[1]),np=$(link_data_d[2])"]
    index_ad = [inds(combined_tensor; tags = link_ad[1]),inds(combined_tensor; tags = link_ad[2])]

    # perform svd on the link and sum schmidt values for entanglement entropy
    u,s,v,spec = svd(combined_tensor,index_ad)
    s_ad = entanglement_entropy(spec.eigs)

    # save data to dictionary (don't need to save cutlink data because it comes from the S_A and S_D)
    ee_data["s_ad"] = s_ad

    gamma = -(ee_data["s_a"] + ee_data["s_b"] + ee_data["s_c"] + ee_data["s_d"]) + (ee_data["s_ab"] + ee_data["s_cd"] + ee_data["s_bc"] + ee_data["s_ad"]) - ee_data["s_abcd"]

    if_save && save_tee(gamma,ee_data,cutlink_data; kwargs...)

    return gamma,ee_data,cutlink_data

end

function save_tee(gamma::Float64,tee_data::Dict{String,Float64},cutlink_data::Dict{String,Tuple{Int64,Int64}}; kwargs...)
    filepath = get(kwargs,:filepath,nothing)

    isnothing(filepath) && error("No filepath provided for saving")

    modify_data(Dict([("tee",gamma),("tee_data",tee_data),("tee_cutlink_data",cutlink_data)]),filepath,"metadata"; output_level=0)
end

function find_particle_type(ttn::TTN.TreeTensorNetwork)
    if join(split(split(string(TTN.tags(TTN.inds(ttn[1,1])[1])),",")[1],"")[2:end]) == "Boson"
        return "Boson"
    else
        return "Fermion"
    end
end

function density_matrix(ttn::TTN.TreeTensorNetwork; kwargs...)
	
    if find_particle_type(ttn) != "Boson"
        return zeros(ComplexF64,2,2)
    end

	creation = "Adag"
	annihilation = "A"
	output_level = get(kwargs, :output_level, false)
	
	lat = TTN.physical_lattice(TTN.network(ttn))
	num_sites = prod(size(lat))
	densmat = zeros(ComplexF64,num_sites,num_sites)
	for i in 1:num_sites
		for j in 1:i
			output_level ? println(i,", ",j) : nothing
			densmat[i,j] = TTN.correlation(ttn,creation,annihilation,i,j)
			densmat[j,i] = conj(densmat[i,j])
		end
	end
	
	return densmat
end

function ft_density_matrix(rho::Matrix{ComplexF64},momentum::Vector{Float64},lx::Int,ly::Int)
    result::ComplexF64 = 0.0

    for s1 in 1:size(rho,1)
        for s2 in 1:size(rho,2)
            s1_coord = coordinate(s1,lx,ly)
            s2_coord = coordinate(s2,lx,ly)
            result += rho[s1,s2] * ft_coeff(s1_coord .- s2_coord,momentum,"A")
        end
    end

    return result
end

function two_point_densmat(rho::Matrix{ComplexF64},lx::Int,ly::Int)
    moms = [n/lx for n in 0:lx]

    result::Vector{Float64} = zeros(Float64,lx+1)

    for (idx,k1) in enumerate(moms)
        result[idx] = real(ft_density_matrix(rho,[k1,0.0],lx,ly))
    end

    return result
end



function local_hilbert_space_dimension(ttn::TTN.TreeTensorNetwork)
    return minimum(TTN.dims(ttn[1,1]))
end

function ft_coeff(phys_site::TTN.Index,momentum::Vector{Float64},op_type::String; kwargs...)
    Lx::Int = kwargs[:Lx]
    Ly::Int = kwargs[:Ly]

    index_tag = string(TTN.tags(phys_site))
    @assert occursin("Site",index_tag)

    lin_ind = parse(Int,match(r"n=(\d+)",index_tag)[1])
    coord_label = coordinate(lin_ind,Lx,Ly)
    return ft_coeff(coord_label,momentum,op_type; kwargs...)
end

function diocane(phys_site::TTN.Index,momentum::Vector{Float64},op_type::String; kwargs...)
    Lx::Int = kwargs[:Lx]
    Ly::Int = kwargs[:Ly]

    index_tag = string(TTN.tags(phys_site))
    @assert occursin("Site",index_tag)

    lin_ind = parse(Int,match(r"n=(\d+)",index_tag)[1])
    coord_label = coordinate(lin_ind,Lx,Ly)

    other_op_type = op_type == "A" ? "Adag" : "A"

    return diocane(coord_label,momentum,other_op_type; kwargs...)
end

function diocane(phys_site::TTN.Index,momentum::Vector{Float64},op_type::String,ttn_size::Tuple{Int,Int}; kwargs...)
    # this is the functional lattice size
    #Lx::Int = kwargs[:Lx]
    #Ly::Int = kwargs[:Ly]

    index_tag = string(TTN.tags(phys_site))
    @assert occursin("Site",index_tag)

    lin_ind = parse(Int,match(r"n=(\d+)",index_tag)[1])
    coord_label = coordinate(lin_ind,ttn_size[1],ttn_size[2])

    other_op_type = op_type == "A" ? "Adag" : "A"

    return diocane(coord_label,momentum,other_op_type; kwargs...)
end

function construct_top_node_environments(ttn1::TTN.TreeTensorNetwork, ttn2::TTN.TreeTensorNetwork, tpo::TTN.MPOWrapper; kwargs...)
    opl::Int = get(kwargs,:output_level,1)

    # need to do some checks at the start
    TTN.move_ortho!(ttn2,(TTN.number_of_layers(ttn2),1))
    TTN.move_ortho!(ttn1,(TTN.number_of_layers(ttn1),1))

    net = ttn1.net

    mapping = tpo.mapping
    ham = tpo.data

    bEnvironment = map(eachindex(net,1)) do pp
        chdnds = TTN.child_nodes(net, (1,pp))
        map(1:TTN.number_of_child_nodes(net, (1,pp))) do nn
            ham[TTN.inverse_mapping(mapping)[chdnds[nn][2]]]
        end
    end

    for ll in Iterators.drop(TTN.eachlayer(net), 1)
        opl > 1 && println("Constructing top node environments for layer $ll")
		bEnvironment_new = Vector{Vector{TTN.ITensor}}(undef, TTN.number_of_tensors(net, ll))
		for pp in eachindex(net, ll)
			n_chds = TTN.number_of_child_nodes(net, (ll,pp))
			bEnvironment_new[pp] = Vector{TTN.ITensor}(undef, n_chds)
		
			for chd in TTN.child_nodes(net, (ll,pp))
                #println("Making environment for child $chd")
				chd_idx = TTN.index_of_child(net, chd)
				Tn1 = ttn1[chd]
                Tn2 = ttn2[chd]
				
				tensorListBottom = map(TTN.child_nodes(net, chd)) do cc
					bEnvironment[chd[2]][TTN.index_of_child(net, cc)]
				end
                #println("At layer $ll and child $chd the tensorListBottom is")
                #display(TTN.inds.(tensorListBottom))
                #display(TTN.tags.(TTN.inds(Tn)))
                #display([TTN.tags.(in) for in in TTN.inds.(tensorListBottom)])
				tlist = vcat(Tn1, tensorListBottom, TTN.prime(TTN.dag(Tn2)))
                #display(prod(prod.(TTN.ITensorMPS.dims.(tlist))))
                #display(TTN.dims.(tlist))
				opt_seq = TTN.optimal_contraction_sequence(tlist)
				bEnvironment_new[pp][chd_idx] = contract(tlist; sequence = opt_seq)
                #println("Now showing after contraction tags \n")
                #display(TTN.dims(bEnvironment_new[pp][chd_idx]))
                #display(TTN.ITensorMPS.dims(bEnvironment_new[pp][chd_idx]))
                #display(TTN.tags.(TTN.inds(bEnvironment_new[pp][chd_idx])))
			end
		end
		bEnvironment = bEnvironment_new
	end
	return only(bEnvironment)
end

function construct_top_node_environments(ttn::TTN.TreeTensorNetwork, tpo::TTN.MPOWrapper; kwargs...)
    opl::Int = get(kwargs,:output_level,1)

	net = ttn.net

    TTN.move_ortho!(ttn,(TTN.number_of_layers(ttn),1))

	mapping = tpo.mapping
	ham = tpo.data
	
	bEnvironment = map(eachindex(net,1)) do pp
        chdnds = TTN.child_nodes(net, (1,pp))
        map(1:TTN.number_of_child_nodes(net, (1,pp))) do nn
          ham[TTN.inverse_mapping(mapping)[chdnds[nn][2]]]
        end
    end
	
	for ll in Iterators.drop(TTN.eachlayer(net), 1)
        opl > 1 && println("Constructing top node environments for layer $ll")
		bEnvironment_new = Vector{Vector{TTN.ITensor}}(undef, TTN.number_of_tensors(net, ll))
		for pp in eachindex(net, ll)
			n_chds = TTN.number_of_child_nodes(net, (ll,pp))
			bEnvironment_new[pp] = Vector{TTN.ITensor}(undef, n_chds)
		
			for chd in TTN.child_nodes(net, (ll,pp))
                #println("Making environment for child $chd")
				chd_idx = TTN.index_of_child(net, chd)
				Tn = ttn[chd]
				
				tensorListBottom = map(TTN.child_nodes(net, chd)) do cc
					bEnvironment[chd[2]][TTN.index_of_child(net, cc)]
				end
                #println("At layer $ll and child $chd the tensorListBottom is")
                #display(TTN.inds.(tensorListBottom))
                #display(TTN.tags.(TTN.inds(Tn)))
                #display([TTN.tags.(in) for in in TTN.inds.(tensorListBottom)])
				tlist = vcat(Tn, tensorListBottom, TTN.prime(TTN.dag(Tn)))
                #tlist = vcat(TTN.prime(Tn), tensorListBottom, TTN.dag(Tn))
                #display(prod(prod.(TTN.ITensorMPS.dims.(tlist))))
                #display(TTN.dims.(tlist))
				opt_seq = TTN.optimal_contraction_sequence(tlist)
				bEnvironment_new[pp][chd_idx] = contract(tlist; sequence = opt_seq)
                #println("Now showing after contraction tags \n")
                #display(TTN.dims(bEnvironment_new[pp][chd_idx]))
                #display(TTN.ITensorMPS.dims(bEnvironment_new[pp][chd_idx]))
                #display(TTN.tags.(TTN.inds(bEnvironment_new[pp][chd_idx])))
			end
		end
		bEnvironment = bEnvironment_new
	end
	return only(bEnvironment)
end

function calculate_mpo_expectation(ttn1::TTN.TreeTensorNetwork, ttn2::TTN.TreeTensorNetwork, tpo::TTN.MPOWrapper; kwargs...)
    opl::Int = get(kwargs,:output_level,1)

	topenvs = construct_top_node_environments(ttn1, ttn2, tpo)
    opl > 1 && println("Finished making environments")
    #display(inds.(topenvs))
	T1 = ttn1[TTN.number_of_layers(ttn1), 1]
    T2 = ttn2[TTN.number_of_layers(ttn2), 1]
	tlist = [T1, topenvs..., TTN.prime(TTN.dag(T2))]
	opt_seq = TTN.optimal_contraction_sequence(tlist)
	return TTN.scalar(contract(tlist; sequence = opt_seq))
end

function calculate_mpo_expectation(ttn::TTN.TreeTensorNetwork, tpo::TTN.MPOWrapper; kwargs...)
    opl::Int = get(kwargs,:output_level,1)
	topenvs = construct_top_node_environments(ttn, tpo; kwargs...)
    opl > 1 && println("Finished making environments")
    #display(inds.(topenvs))
	T = ttn[TTN.number_of_layers(ttn), 1]
	tlist = [T, topenvs..., TTN.prime(TTN.dag(T))]
	opt_seq = TTN.optimal_contraction_sequence(tlist)
	return TTN.scalar(contract(tlist; sequence = opt_seq))
end

function zigzag_curve(lat::TTN.AbstractLattice)
    return zigzag_curve(lat.dims[1],lat.dims[2])
end

function zigzag_curve(lx::Int,Ly::Int)
    curve::Vector{Int} = zeros(Int,Ly*lx)
    num_quadrants::Int = Int(lx*Ly/4)

    starting_point = (1,1)
    for q in 1:num_quadrants
        # add starting point
        curve[4*(q-1) + 1] = linear_index(starting_point,lx,Ly)

        # move right one point
        p2 = starting_point .+ (1,0)
        curve[4*(q-1) + 2] = linear_index(p2,lx,Ly)

        # move diagonally up and left
        p3 = p2 .+ (-1,1)
        curve[4*(q-1) + 3] = linear_index(p3,lx,Ly)

        # move right one point
        p4 = p3 .+ (1,0)
        curve[4*(q-1) + 4] = linear_index(p4,lx,Ly)

        # reset starting point
        isodd(q) && (starting_point = p4 .+ (1,-1))

        haystack = 8 .* [1,2,3,4]
        if iseven(q)
            if q in haystack
                renorm_val = Int(q/8)
                isodd(renorm_val) && (starting_point = p4 .+ (-7,1))
                iseven(renorm_val) && (starting_point = p4 .+ (1,-7))
            else
                if (q-2) % 4 == 0
                    starting_point = p4 .+ (-3,1)
                elseif (q-2) % 4 == 2
                    starting_point = p4 .+ (1,-3)
                end
            end
        end
    end

    return curve
end

function make_mpowrapper(mpo::TTN.MPO, lat::L; mapping::Vector{Int} = collect(eachindex(lat))) where{L}
    @assert TTN.is_physical(lat)
    @assert length(lat) == length(mpo)
    #@assert isone(dimensionality(lat))
    idx_lat = TTN.siteinds(lat)

    mpoc = TTN.deepcopy(mpo)
    idx_mpo = last.(TTN.siteinds(mpoc,plev = 0))
    println("Starting wrapping")
    for (idx,jj) in enumerate(mapping)
        #println("working on wrapping site $jj")
        sj_lat = idx_lat[jj]
        sj_mpo = idx_mpo[jj]
        mpoc[idx] = TTN.replaceinds!(mpoc[jj], sj_mpo => sj_lat, TTN.prime(sj_mpo) => TTN.prime(sj_lat))
    end
    return TTN.MPOWrapper{L, TTN.MPO}(lat, mpoc, mapping)
end

function easy_mpowrapper(mpo::TTN.MPO, lat::L; mapping::Vector{Int} = collect(eachindex(lat))) where{L}
    @assert TTN.is_physical(lat)
    @assert length(lat) == length(mpo)

    return TTN.MPOWrapper{L, TTN.MPO}(lat, mpo, mapping)
end

function adag_matrix(hilbdim::Int)
    mat = zeros(Float64,hilbdim,hilbdim)
    for i in 1:(hilbdim-1)
        mat[i,i+1] = sqrt(i)
    end
    return mat
end

function a_matrix(hilbdim::Int)
    mat = zeros(Float64,hilbdim,hilbdim)
    for i in 2:hilbdim
        mat[i,i-1] = sqrt(i-1)
    end
    return mat
end

function n_matrix(hilbdim::Int)
    mat = zeros(Float64,hilbdim,hilbdim)
    for i in 1:hilbdim
        mat[i,i] = i-1
    end
    return mat
end

function make_qnset(op_string::String)
    if op_string == "A"
        return [TTN.QN("Number",0)=>1,TTN.QN("Number",-1)=>1]
    elseif op_string == "Adag"
        return [TTN.QN("Number",0)=>1,TTN.QN("Number",1)=>1]
    elseif op_string == "N"
        return [TTN.QN("Number",0)=>2]
    elseif op_string == "2pt"
        return [TTN.QN("Number",0)=>2,TTN.QN("Number",-1)=>1,TTN.QN("Number",1)=>1]
    else
        error("Invalid operator string")
    end
end

function build_W_singlepoint(op_type::String,coeff::ComplexF64,hilbdim::Int)
    if op_type == "Adag"
        mat = adag_matrix(hilbdim)
    elseif op_type == "A"
        mat = a_matrix(hilbdim)
    elseif op_type == "N"
        mat = n_matrix(hilbdim)
    else
        error("Invalid operator type")
    end

    A11 = I(hilbdim)
    A21 = zeros(hilbdim,hilbdim)
    A12 = coeff * mat
    A22 = I(hilbdim)
    M = Array{ComplexF64, 4}(undef, 2,hilbdim,2,hilbdim)
    M[1,:,1,:] = A11
    M[2,:,1,:] = A21
    M[1,:,2,:] = A12
    M[2,:,2,:] = A22

    return M
end

function build_links_singlepoint(op_type::String,L::Int; kwargs...)
    mapping = kwargs[:mapping]
    links = Vector{TTN.Index}(undef,L+1)

    qnset_start = make_qnset(op_type)

    # do zeroth link
    links[1] = TTN.Index(qnset_start; tags="Link,Left=Start,Right=$(mapping[1])")

    # do interior links
    for i in 1:L-1
        qnset = make_qnset(op_type)
        #println("At link $i which is physical site $(mapping[i]) which has left tag $(string(mapping[i])) and right tag $(string(mapping[i+1]))")
        left_tag = string(mapping[i])
        right_tag = string(mapping[i+1])
        links[i+1] = TTN.Index(qnset; tags="Link,Left=$left_tag,Right=$right_tag")
    end

    # do last link
    qnset_end = make_qnset(op_type)
    links[L+1] = TTN.Index(qnset_end; tags="Link,Left=$(mapping[L]),Right=End")

    return links
end

function single_point_mpo(wavefunc::TTN.TreeTensorNetwork,op_type::String; kwargs...)
    opl::Int = get(kwargs,:output_level,1)
    mom::Vector{Float64} = get(kwargs,:momentum,[0.0,0.0])
    which_coeff::Function = get(kwargs,:which_coeff,diocane)
    coeff_kwargs::NamedTuple = kwargs[:coeff_kwargs]

    hilbdim = local_hilbert_space_dimension(wavefunc)

    mapping = kwargs[:mapping]

    phys_sites = TTN.sites(wavefunc)[mapping]

    links = build_links_singlepoint(op_type,length(phys_sites); mapping=mapping)

    tensor_train = Vector{TTN.ITensor}(undef,length(phys_sites))
    for (idx,s) in enumerate(phys_sites)
        opl > 1 && println("Working on Physical Site $(TTN.tags(s))")

        coeff::ComplexF64 = which_coeff(s,mom,op_type; coeff_kwargs...)

        mat = build_W_singlepoint(op_type,coeff,hilbdim)

        local_inds = [links[idx],TTN.dag(s),TTN.dag(links[idx+1]),TTN.prime(s)]
        #local_inds = [links[idx],s,TTN.dag(links[idx+1]),TTN.prime(TTN.dag(s))]

        mit = TTN.ITensor(mat,local_inds)

        if idx == 1
			mit = mit * TTN.dag(TTN.onehot(links[1] => 1))
		elseif idx == length(phys_sites)
			mit = mit * (TTN.onehot(links[end] => 2))
		end
        
        tensor_train[idx] = mit
    end

    return TTN.MPO(tensor_train)
end

# same function for restricted size TTNs means set coeffs on dead sites to zero
function single_point_mpo(wavefunc::TTN.TreeTensorNetwork,op_type::String,restricted_size::Vector{Int}; kwargs...)
    opl::Int = get(kwargs,:output_level,1)
    mom::Vector{Float64} = get(kwargs,:momentum,[0.0,0.0])
    which_coeff::Function = get(kwargs,:which_coeff,diocane)
    coeff_kwargs::NamedTuple = kwargs[:coeff_kwargs]

    ttn_size = get_lattice_dims(wavefunc)

    dead_sites_linear = find_dead_linear_sites(restricted_size,ttn_size)

    hilbdim = local_hilbert_space_dimension(wavefunc)

    mapping = kwargs[:mapping]

    phys_sites = TTN.sites(wavefunc)[mapping]

    links = build_links_singlepoint(op_type,length(phys_sites); mapping=mapping)

    tensor_train = Vector{TTN.ITensor}(undef,length(phys_sites))
    for (idx,s) in enumerate(phys_sites)
        opl > 1 && println("Working on Physical Site $(TTN.tags(s))")

        #if idx in dead_sites_linear
        #    coeff::ComplexF64 = 0.0
        #else
        coeff = which_coeff(s,mom,op_type,ttn_size; coeff_kwargs...)
        #end

        mat = build_W_singlepoint(op_type,coeff,hilbdim)

        local_inds = [links[idx],TTN.dag(s),TTN.dag(links[idx+1]),TTN.prime(s)]
        #local_inds = [links[idx],s,TTN.dag(links[idx+1]),TTN.prime(TTN.dag(s))]

        mit = TTN.ITensor(mat,local_inds)

        if idx == 1
			mit = mit * TTN.dag(TTN.onehot(links[1] => 1))
		elseif idx == length(phys_sites)
			mit = mit * (TTN.onehot(links[end] => 2))
		end
        
        tensor_train[idx] = mit
    end

    return TTN.MPO(tensor_train)
end

function two_point_mpo(wavefunc::TTN.TreeTensorNetwork; kwargs...)
    opl::Int = get(kwargs,:output_level,1)

    k1::Vector{Float64} = get(kwargs,:momentum1,[0.0,0.0])
    k2::Vector{Float64} = get(kwargs,:momentum2,[0.0,0.0])
    mapping::Vector{Int} = get(kwargs,:mapping,collect(1:TTN.number_of_sites(wavefunc.net)))

    creat = single_point_mpo(wavefunc,"Adag"; momentum=k1,mapping=mapping,kwargs...)
    opl > 0 && println("Made Creation")
    annih = single_point_mpo(wavefunc,"A"; momentum=k2,mapping=mapping,kwargs...)
    opl > 0 && println("Made Annihilation")

    return apply(creat,annih)
end

function two_point_mpowrapped(wavefunc::TTN.TreeTensorNetwork,momentum1::Vector{Float64},momentum2::Vector{Float64}; kwargs...)
    lat = TTN.physical_lattice(wavefunc.net)
    Lx,Ly = size(lat)
    mapss = ttn_2d_mapping([Lx,Ly])#zigzag_curve(Lx,Ly)
    coeff_kwargs = get(kwargs,:coeff_kwargs,(Lx=Lx,Ly=Ly,))

    twopt = two_point_mpo(wavefunc; momentum1 = momentum1, momentum2 = momentum2, mapping = mapss, coeff_kwargs=coeff_kwargs, kwargs...)
    twopt_wrapped = easy_mpowrapper(twopt, lat; mapping=mapss)
    return twopt_wrapped
end

function two_point(wavefunc::TTN.TreeTensorNetwork,momentum1::Vector{Float64},momentum2::Vector{Float64}; kwargs...)
    lat = TTN.physical_lattice(wavefunc.net)
    Lx,Ly = size(lat)
    mapss = ttn_2d_mapping([Lx,Ly])#zigzag_curve(Lx,Ly)
    coeff_kwargs = get(kwargs,:coeff_kwargs,(Lx=Lx,Ly=Ly,))

    twop = two_point_mpo(wavefunc; momentum1 = momentum1, momentum2 = momentum2, mapping = mapss, coeff_kwargs=coeff_kwargs, kwargs...)
    twop_wrapped = easy_mpowrapper(twop, lat; mapping=mapss)
    return calculate_mpo_expectation(wavefunc, twop_wrapped; kwargs...)
end

function two_point(wavefuncs::Vector{TTN.TreeTensorNetwork},momentum1::Vector{Float64},momentum2::Vector{Float64}; kwargs...)

    lat = TTN.physical_lattice(wavefuncs[1].net)
    Lx,Ly = size(lat)
    mapss = ttn_2d_mapping([Lx,Ly])#zigzag_curve(Lx,Ly)
    coeff_kwargs = get(kwargs,:coeff_kwargs,(Lx=Lx,Ly=Ly,))

    twop = two_point_mpo(wavefuncs[1]; momentum1 = momentum1, momentum2 = momentum2, mapping = mapss, coeff_kwargs=coeff_kwargs, kwargs...)
    twop_wrapped = easy_mpowrapper(twop, lat; mapping=mapss)

    mat::Matrix{ComplexF64} = zeros(Float64,length(wavefuncs),length(wavefuncs))
    for i in 1:length(wavefuncs)
        for j in 1:length(wavefuncs)
            mat[i,j] = calculate_mpo_expectation(wavefuncs[i], wavefuncs[j], twop_wrapped; kwargs...)
        end 
    end

    #display(mat)

    return eigvals(mat)
end

function two_point(wavefunc::TTN.TreeTensorNetwork; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,false)
    opl::Int = get(kwargs,:outpu_level,1)

    Lx,Ly = get_lattice_dims(wavefunc)

    momenta = [n/Ly for n in 0:Lx-1]
    twopt_vals = zeros(Float64,Lx,Lx)
    for (idx1,k1) in enumerate(momenta)
        for (idx2,k2) in enumerate(momenta)
            opl > 0 && println("Working on momenta $(k1) and $(k2)")
            twopt_vals[idx1,idx2] = abs(two_point(wavefunc,[0.0,k1],[0.0,k2]; kwargs...))
        end
    end

    if_plot && plot_four_point(twopt_vals; kwargs...,if_2pt=true) 

    return twopt_vals
end

function four_point_mpo(wavefunc::TTN.TreeTensorNetwork; kwargs...)
    opl::Int = get(kwargs,:output_level,1)

    k1::Vector{Float64} = get(kwargs,:momentum1,[0.0,0.0])
    k2::Vector{Float64} = get(kwargs,:momentum2,[0.0,0.0])
    mapping::Vector{Int} = get(kwargs,:mapping,collect(1:TTN.number_of_sites(wavefunc.net)))

    #println("Momenta are $(k1[2]) and $(k2[2])")

    creat1 = single_point_mpo(wavefunc,"Adag"; momentum=k1,mapping=mapping, kwargs...)
    #println("Made Creation 1")
    creat2 = single_point_mpo(wavefunc,"Adag"; momentum=k2,mapping=mapping, kwargs...)
    #println("Made Creation 2")
    annih1 = single_point_mpo(wavefunc,"A"; momentum=k2,mapping=mapping, kwargs...)
    #println("Made Annihilation 1")
    annih2 = single_point_mpo(wavefunc,"A"; momentum=k1,mapping=mapping, kwargs...)
    #println("Made Annihilation 2")
    opl > 1 && println("Made Suboperators")

    return TTN.ITensors.apply(TTN.ITensors.apply(creat1, creat2), TTN.ITensors.apply(annih1, annih2))
    #bothcreats = TTN.replaceprime(TTN.contract(creat1',creat2; alg="naive", truncate=false), 2 => 1)
    #bothannihs = TTN.replaceprime(TTN.contract(annih1',annih2; alg="naive", truncate=false), 2 => 1)
    #return TTN.replaceprime(TTN.contract(bothcreats',bothannihs; alg="naive", truncate=false), 2 => 1)
end

function four_point_mpo(wavefunc::TTN.TreeTensorNetwork, restricted_size::Vector{Int}; kwargs...)
    opl::Int = get(kwargs,:output_level,1)

    k1::Vector{Float64} = get(kwargs,:momentum1,[0.0,0.0])
    k2::Vector{Float64} = get(kwargs,:momentum2,[0.0,0.0])
    mapping::Vector{Int} = get(kwargs,:mapping,collect(1:TTN.number_of_sites(wavefunc.net)))

    #println("Momenta are $(k1[2]) and $(k2[2])")

    creat1 = single_point_mpo(wavefunc,"Adag",restricted_size; momentum=k1,mapping=mapping, kwargs...)
    #println("Made Creation 1")
    creat2 = single_point_mpo(wavefunc,"Adag",restricted_size; momentum=k2,mapping=mapping, kwargs...)
    #println("Made Creation 2")
    annih1 = single_point_mpo(wavefunc,"A",restricted_size; momentum=k2,mapping=mapping, kwargs...)
    #println("Made Annihilation 1")
    annih2 = single_point_mpo(wavefunc,"A",restricted_size; momentum=k1,mapping=mapping, kwargs...)
    #println("Made Annihilation 2")
    opl > 1 && println("Made Suboperators")

    return TTN.ITensors.apply(TTN.ITensors.apply(creat1, creat2), TTN.ITensors.apply(annih1, annih2))
    #bothcreats = TTN.replaceprime(TTN.contract(creat1',creat2; alg="naive", truncate=false), 2 => 1)
    #bothannihs = TTN.replaceprime(TTN.contract(annih1',annih2; alg="naive", truncate=false), 2 => 1)
    #return TTN.replaceprime(TTN.contract(bothcreats',bothannihs; alg="naive", truncate=false), 2 => 1)
end

function four_point_mpowrapped(wavefunc::TTN.TreeTensorNetwork,momentum1::Vector{Float64},momentum2::Vector{Float64}; kwargs...)
    lat = TTN.physical_lattice(wavefunc.net)
    Lx,Ly = size(lat)
    mapss = ttn_2d_mapping([Lx,Ly])#zigzag_curve(Lx,Ly)
    coeff_kwargs = get(kwargs,:coeff_kwargs,(Lx=Lx,Ly=Ly,))

    fourpt = four_point_mpo(wavefunc; momentum1 = momentum1, momentum2 = momentum2, mapping = mapss, coeff_kwargs=coeff_kwargs, kwargs...)
    fourpt_wrapped = easy_mpowrapper(fourpt, lat; mapping=mapss)
    return fourpt_wrapped
end

function four_point(wavefunc::TTN.TreeTensorNetwork,momentum1::Vector{Float64},momentum2::Vector{Float64}; kwargs...)
    lat = TTN.physical_lattice(wavefunc.net)
    Lx,Ly = size(lat)
    mapss = ttn_2d_mapping([Lx,Ly])#zigzag_curve(Lx,Ly)
    coeff_kwargs = get(kwargs,:coeff_kwargs,(Lx=Lx,Ly=Ly,))

    fourpt = four_point_mpo(wavefunc; momentum1 = momentum1, momentum2 = momentum2, mapping = mapss, coeff_kwargs=coeff_kwargs, kwargs...)
    fourpt_wrapped = easy_mpowrapper(fourpt, lat; mapping=mapss)

    #=matver = focking_matrix(wavefunc,fourpt_wrapped,lattice_params["full_basis"]; kwargs...)
    ed_wavefunc = focking_vector(wavefunc,lattice_params["full_basis"])
    return abs(adjoint(ed_wavefunc) * matver * ed_wavefunc)=#

    return abs(calculate_mpo_expectation(wavefunc, fourpt_wrapped; kwargs...))
end

function four_point(wavefunc::TTN.TreeTensorNetwork,momentum1::Vector{Float64},momentum2::Vector{Float64},restricted_size::Vector{Int}; kwargs...)
    lat = TTN.physical_lattice(wavefunc.net)
    Lx,Ly = size(lat)
    mapss = ttn_2d_mapping([Lx,Ly])#zigzag_curve(Lx,Ly)
    coeff_kwargs = get(kwargs,:coeff_kwargs,(Lx=restricted_size[1],Ly=restricted_size[2],))

    fourpt = four_point_mpo(wavefunc, restricted_size; momentum1 = momentum1, momentum2 = momentum2, mapping = mapss, coeff_kwargs=coeff_kwargs, kwargs...)
    fourpt_wrapped = easy_mpowrapper(fourpt, lat; mapping=mapss)

    #=matver = focking_matrix(wavefunc,fourpt_wrapped,lattice_params["full_basis"]; kwargs...)
    ed_wavefunc = focking_vector(wavefunc,lattice_params["full_basis"])
    return abs(adjoint(ed_wavefunc) * matver * ed_wavefunc)=#

    return abs(calculate_mpo_expectation(wavefunc, fourpt_wrapped; kwargs...))
end

function four_point(wavefuncs::Vector{TTN.TreeTensorNetwork},momentum1::Vector{Float64},momentum2::Vector{Float64}; kwargs...)

    lat = TTN.physical_lattice(wavefuncs[1].net)
    Lx,Ly = size(lat)
    mapss = ttn_2d_mapping([Lx,Ly])#zigzag_curve(Lx,Ly)
    coeff_kwargs = get(kwargs,:coeff_kwargs,(Lx=Lx,Ly=Ly,))

    fourpt = four_point_mpo(wavefuncs[1]; momentum1 = momentum1, momentum2 = momentum2, mapping = mapss, coeff_kwargs=coeff_kwargs, kwargs...)
    fourpt_wrapped = easy_mpowrapper(fourpt, lat; mapping=mapss)

    mat::Matrix{ComplexF64} = zeros(Float64,length(wavefuncs),length(wavefuncs))
    for i in 1:length(wavefuncs)
        for j in 1:length(wavefuncs)
            mat[i,j] = calculate_mpo_expectation(wavefuncs[i], wavefuncs[j], fourpt_wrapped; kwargs...)
        end 
    end

    return abs.(eigvals(mat))
end

function four_point(wavefunc::TTN.TreeTensorNetwork; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,false)
    opl::Int = get(kwargs,:output_level,1)

    Lx,Ly = get_lattice_dims(wavefunc)

    momenta = [n/Ly for n in 0:Lx-1]
    fourpt_vals = zeros(Float64,Lx,Lx)
    for (idx1,k1) in enumerate(momenta)
        for (idx2,k2) in enumerate(momenta)
            opl > 0 && println("Working on momenta $(k1) and $(k2)")
            fourpt_vals[idx1,idx2] = four_point(wavefunc,[0.0,k1],[0.0,k2]; kwargs...)
        end
    end

    if_plot && plot_four_point(fourpt_vals; kwargs...) 

    return fourpt_vals
end

function four_point(wavefunc::TTN.TreeTensorNetwork, restricted_size::Vector{Int}; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,false)
    opl::Int = get(kwargs,:output_level,1)

    Lx,Ly = restricted_size

    momenta = [n/Ly for n in 0:Lx-1]
    fourpt_vals = zeros(Float64,Lx,Lx)
    display(fourpt_vals)
    for (idx1,k1) in enumerate(momenta)
        for (idx2,k2) in enumerate(momenta)
            opl > 0 && println("Working on momenta $(k1) and $(k2)")
            fourpt_vals[idx1,idx2] = four_point(wavefunc,[0.0,k1],[0.0,k2],restricted_size; kwargs...)
        end
    end

    if_plot && plot_four_point(fourpt_vals; kwargs...) 

    return fourpt_vals
end

function four_point(wavefuncs::Vector{TTN.TreeTensorNetwork}; kwargs...)
    if_plot::Bool = get(kwargs,:if_plot,false)
    opl::Int = get(kwargs,:output_level,1)

    Lx,Ly = get_lattice_dims(wavefuncs[1])

    momenta = [n/Ly for n in 0:Lx-1]
    fourpt_vals = [zeros(Float64,Lx,Lx) for i in 1:length(wavefuncs)]
    for (idx1,k1) in enumerate(momenta)
        for (idx2,k2) in enumerate(momenta)
            opl > 0 && println("Working on momenta $(k1) and $(k2)")
            result = four_point(wavefuncs,[0.0,k1],[0.0,k2]; kwargs...)
            for i in 1:length(wavefuncs)
                fourpt_vals[i][idx1,idx2] = result[i]
            end
        end
    end

    if if_plot
        for i in 1:length(wavefuncs)
            plot_four_point(fourpt_vals[i]; kwargs...)
        end
    end

    return fourpt_vals
end

function four_point_mpo_real(wavefunc::TTN.TreeTensorNetwork,s1::Int,s2::Int,s3::Int,s4::Int; kwargs...)

    #mapping::Vector{Int} = get(kwargs,:mapping,collect(1:TTN.number_of_sites(wavefunc.net)))

    phys_sites = TTN.sites(wavefunc)

    os_creat1 = OpSum()
    os_creat1 += "Adag",s1
    creat1 = TTN.MPO(os_creat1,phys_sites)

    os_creat2 = OpSum()
    os_creat2 += "Adag",s2
    creat2 = TTN.MPO(os_creat2,phys_sites)

    os_annih1 = OpSum()
    os_annih1 += "A",s3
    annih1 = TTN.MPO(os_annih1,phys_sites)

    os_annih2 = OpSum()
    os_annih2 += "A",s4
    annih2 = TTN.MPO(os_annih2,phys_sites)

    #println("Made Suboperators")

    term1 = apply(apply(creat1,creat2),apply(annih1,annih2))
    return term1

end

function four_point_real(wavefunc::TTN.TreeTensorNetwork, momentum1::Vector{Float64}, momentum2::Vector{Float64}; kwargs...)

    lat = TTN.physical_lattice(wavefunc.net)
    Lx,Ly = size(lat)
    alpha = 1/Ly
    mapss = zigzag_curve(Lx,Ly)

    mval = Int(momentum1[2] * Ly)
    mval2 = Int(momentum2[2] * Ly)

    if mval > Lx || mval2 > Lx
        error("Momentum out of bounds: m=$mval mp=$mval2 Lx=$Lx")
    elseif mval2 == 0 || mval == 0
        error("Momentum cannot be zero: m=$mval mp=$mval2")
    end

    rez::ComplexF64 = 0.0
    for s1 in 1:Ly
        coeff1::ComplexF64 = ft_coeff_alberto([mval,s1],momentum1,"Adag",Lx,Ly,mval,alpha)
        lin1 = linear_index([mval,s1],Lx,Ly)
        for s2 in 1:Ly
            println("Working on site $(s1) and $(s2)")
            coeff2::ComplexF64 = ft_coeff_alberto([mval2,s2],momentum2,"Adag",Lx,Ly,mval2,alpha)
            lin2 = linear_index([mval2,s2],Lx,Ly)
            for s3 in 1:Ly
                coeff3::ComplexF64 = ft_coeff_alberto([mval2,s3],momentum2,"A",Lx,Ly,mval2,alpha)
                lin3 = linear_index([mval2,s3],Lx,Ly)
                for s4 in 1:Ly
                    lin4 = linear_index([mval,s4],Lx,Ly)
                    coeff4::ComplexF64 = ft_coeff_alberto([mval,s4],momentum1,"A",Lx,Ly,mval,alpha)

                    fourpt = four_point_mpo_real(wavefunc,lin1,lin2,lin3,lin4; mapping = mapss)
                    fourpt_wrapped = easy_mpowrapper(fourpt, lat; mapping=mapss)

                    coeff = coeff1 * coeff2 * coeff3 * coeff4
                    local_exppart = calculate_mpo_expectation(wavefunc, fourpt_wrapped; opl=0)

                    println("At $s1 $s2 $s3 $s4 coeff is: $(round(coeff,digits=10))")
                    println("expectation value = $(round(local_exppart,digits=10))")

                    rez += coeff * local_exppart
                end
            end
        end
    end

    return rez
end

function two_point_mpo_real(wavefunc::TTN.TreeTensorNetwork,s1::Int,s2::Int; kwargs...)

    #mapping::Vector{Int} = get(kwargs,:mapping,collect(1:TTN.number_of_sites(wavefunc.net)))

    phys_sites = TTN.sites(wavefunc)

    os_creat1 = OpSum()
    os_creat1 += "Adag",s1
    creat1 = TTN.MPO(os_creat1,phys_sites)

    os_annih1 = OpSum()
    os_annih1 += "A",s2
    annih1 = TTN.MPO(os_annih1,phys_sites)

    #println("Made Suboperators")

    term1 = apply(creat1,annih1)
    return term1

end

function two_point_real(wavefunc::TTN.TreeTensorNetwork, momentum1::Vector{Float64}, momentum2::Vector{Float64}; kwargs...)

    lat = TTN.physical_lattice(wavefunc.net)
    Lx,Ly = size(lat)
    alpha = 1/Ly
    mapss = zigzag_curve(Lx,Ly)

    mval = Int(momentum1[2] * Ly)
    mval2 = Int(momentum2[2] * Ly)

    if mval > Lx || mval2 > Lx
        error("Momentum out of bounds: m=$mval mp=$mval2 Lx=$Lx")
    elseif mval2 == 0 || mval == 0
        error("Momentum cannot be zero: m=$mval mp=$mval2")
    end

    rez::ComplexF64 = 0.0
    for s1 in 1:Ly
        coord1 = [mval,s1]
        coeff1::ComplexF64 = ft_coeff_alberto(coord1,momentum1,"Adag",Lx,Ly,mval,alpha)
        lin1 = linear_index(coord1,Lx,Ly)
        for s2 in 1:Ly
            println("Working on site $(s1) and $(s2)")
            coord2 = [mval2,s2]
            coeff2::ComplexF64 = ft_coeff_alberto(coord2,momentum2,"A",Lx,Ly,mval2,alpha)
            lin2 = linear_index(coord2,Lx,Ly)

            twopt = two_point_mpo_real(wavefunc,lin1,lin2; mapping = mapss)
            twopt_wrapped = easy_mpowrapper(twopt, lat; mapping=mapss)

            local_exppart = calculate_mpo_expectation(wavefunc, twopt_wrapped; opl=0)
            local_rez = (coeff1 * coeff2) * local_exppart

            #println("At site $coord1 and $coord2: expectation is $(round(local_exppart, digits=8))")

            rez += local_rez
        end
    end

    return rez
end













































"fin"
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
    cap = get(kwargs,:if_cap,nothing)
    layers_down = get(kwargs,:layers_down,0)

    numlayers = TTNKit.number_of_layers(psi)
    TTNKit.move_ortho!(psi,(numlayers - layers_down,1))
    top_tensor = psi[numlayers - layers_down,1]
    idx_left = inds(top_tensor; tags = "Link,nl=$(numlayers-1-layers_down),np=1")
    u,s,v,spec = svd(top_tensor,idx_left)

    result = isnothing(cap) ? spec.eigs : spec.eigs[1:cap]

    if_save && save_spatialentanglementspectrum(result,filepath)

    return result
end

function save_spatialentanglementspectrum(spec::Vector{Float64},filepath::Nothing)
    error("File Path not provided")
end

function save_spatialentanglementspectrum(spec::Vector{Float64},filepath::String)
    modify_data_jld2(Dict([("entanglement_spectrum",spec)]),filepath,"metadata"; output_level=0)
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

function tee(psi::TTNKit.TreeTensorNetwork,top_layer::Int64; kwargs...)
    if_save::Bool = get(kwargs,:if_save,false)

    numlayers = TTNKit.number_of_layers(psi)

    if numlayers < top_layer
        error("Top layer is higher than the number of layers in the TTN")
    end

    ee_data::Dict{String,Float64} = Dict()
    cutlink_data::Dict{String,Tuple{Int64,Int64}} = Dict()

    net = TTNKit.network(psi)

    # start with S_ABDC at the top node
    top_node = (top_layer,1)

    # make sure the ortho_center is at the cutting tensor
    TTNKit.move_ortho!(psi,top_node)

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
    middle_layer_children = TTNKit.child_nodes(net,middle_node)
    for (i,middle_child) in enumerate(middle_layer_children)
        
        # move ortho_center to the tensor site to cut
        TTNKit.move_ortho!(psi,middle_node)

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
        baby_layer_children = TTNKit.child_nodes(net,middle_child)
        for (j,baby_child) in enumerate(baby_layer_children)
                
                # move ortho_center to the tensor site to cut
                TTNKit.move_ortho!(psi,middle_child)
    
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
    TTNKit.move_ortho!(psi,middle_node)
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

    modify_data_jld2(Dict([("tee",gamma),("tee_data",tee_data),("tee_cutlink_data",cutlink_data)]),filepath,"metadata"; output_level=0)
end

function construct_top_node_environments(ttn::TTNKit.TreeTensorNetwork, tpo::TTNKit.MPOWrapper)
	
	net = ttn.net

	n_sites = TTNKit.number_of_sites(net)
	n_tensors = TTNKit.number_of_tensors(net) + n_sites

	mapping = tpo.mapping
	ham = tpo.data
	
	bEnvironment = map(eachindex(net,1)) do pp
        chdnds = TTNKit.child_nodes(net, (1,pp))
        map(1:TTNKit.number_of_child_nodes(net, (1,pp))) do nn
          ham[TTNKit.inverse_mapping(mapping)[chdnds[nn][2]]]
        end
    end
	
	for ll in Iterators.drop(TTNKit.eachlayer(net), 1)
        println("Constructing top node environments for layer $ll")
		bEnvironment_new = Vector{Vector{ITensor}}(undef, TTNKit.number_of_tensors(net, ll))
		for pp in eachindex(net, ll)
			n_chds = TTNKit.number_of_child_nodes(net, (ll,pp))
			bEnvironment_new[pp] = Vector{ITensor}(undef, n_chds)
		
			for chd in TTNKit.child_nodes(net, (ll,pp))
                println("Making environment for child $chd")
				chd_idx = TTNKit.index_of_child(net, chd)
				Tn = ttn[chd]
				
				tensorListBottom = map(TTNKit.child_nodes(net, chd)) do cc
					bEnvironment[chd[2]][TTNKit.index_of_child(net, cc)]
				end
                #println("At layer $ll and child $chd the tensorListBottom is")
                #display(inds.(tensorListBottom))
				tlist = vcat(Tn, tensorListBottom, prime(dag(Tn)))
                #display(prod(prod.(TTNKit.ITensors.dims.(tlist))))
                #display(inds.(tlist))
				opt_seq = ITensors.optimal_contraction_sequence(tlist)
				bEnvironment_new[pp][chd_idx] = contract(tlist; sequence = opt_seq)
                #println("Now showing after contraction tags \n")
                #display(inds(bEnvironment_new[pp][chd_idx]))
                #display(prod(TTNKit.ITensors.dims(bEnvironment_new[pp][chd_idx])))
			end
		end
		bEnvironment = bEnvironment_new
	end
	return only(bEnvironment)
end

function calculate_mpo_expectation(ttn::TTNKit.TreeTensorNetwork, tpo::TTNKit.MPOWrapper)
	topenvs = construct_top_node_environments(ttn, tpo)
    println("Finished making environments")
    #display(inds.(topenvs))
	T = ttn[TTNKit.number_of_layers(ttn), 1]
	tlist = [T, topenvs..., prime(dag(T))]
	opt_seq = ITensors.optimal_contraction_sequence(tlist)
	return scalar(contract(tlist; sequence = opt_seq))
end

function make_mpowrapper(mpo::MPO, lat::L; mapping::Vector{Int} = collect(eachindex(lat))) where{L}
    @assert TTNKit.is_physical(lat)
    @assert length(lat) == length(mpo)
    #@assert isone(dimensionality(lat))
    idx_lat = TTNKit.siteinds(lat)

    mpoc = TTNKit.deepcopy(mpo)
    idx_mpo = last.(TTNKit.siteinds(mpoc,plev = 0))
    println("Starting wrapping")
    for (idx,jj) in enumerate(mapping)
        println("working on wrapping site $jj")
        sj_lat = idx_lat[jj]
        sj_mpo = idx_mpo[jj]
        mpoc[idx] = replaceinds!(mpoc[jj], sj_mpo => sj_lat, prime(sj_mpo) => prime(sj_lat))
    end
    return TTNKit.MPOWrapper{L, MPO, TTNKit.ITensorsBackend}(lat, mpoc, mapping)
end

#=function build_W_singlepoint(which_ladder::Int,coeff::ComplexF64)
    mat::Array{ComplexF64} = zeros(ComplexF64,2,2,2,2)
    mat[1,:,1,:] = I(2)
    mat[2,:,2,:] = I(2)
    if which_ladder == 1
        mat[2,:,1,:] = [0.0 0.0; 1.0 0.0]
        #mat[2,:,1,:] = zeros(2,2)
    elseif which_ladder == -1
        mat[2,:,1,:] = [0.0 1.0; 0.0 0.0]
        #mat[1,:,2,:] = zeros(2,2)
    elseif which_ladder == 0
        mat[1,:,2,:] = [0.0 0.0; 0.0 1.0] # counts the total particle number
    else
        error("Invalid ladder type")
    end

    return mat
end

function build_W_4pt(coeff::ComplexF64)
    error("Not implemented yet")
end

function build_W_2pt(coeff::ComplexF64)
    error("Still not working with quantum numbers")
    
    mat::Array{ComplexF64} = zeros(ComplexF64,4,4,2,2)

    mat[1,1,:,:] = I(2)
    mat[2,2,:,:] = I(2)
    #mat[3,3,:,:] = I(2)
    mat[4,4,:,:] = I(2)

    mat[1,2,:,:] = [0.0 0.0; 0.0 1.0] # N
    mat[1,3,:,:] = [0.0 0.0; 0.0 1.0] # N
    mat[2,4,:,:] = [0.0 0.0; 0.0 1.0] # N
    mat[3,4,:,:] = I(2) # Id

    return mat
end

function build_W(op_string::String,coeff::ComplexF64)
    if op_string == "A"
        return build_W_singlepoint(-1,coeff)
    elseif op_string == "Adag"
        return build_W_singlepoint(1,coeff)
    elseif op_string == "N"
        return build_W_singlepoint(0,coeff)
    elseif op_string == "2pt"
        return build_W_2pt(coeff)
    elseif op_string == "4pt"
        return build_W_4pt(coeff)
    else
        error("Invalid operator string")
    end
end=#

function make_qnset(op_string::String)
    if op_string == "A"
        return [QN("Number",0)=>1,QN("Number",-1)=>1]
    elseif op_string == "Adag"
        return [QN("Number",0)=>1,QN("Number",1)=>1]
    elseif op_string == "N"
        return [QN("Number",0)=>2]
    elseif op_string == "2pt"
        return [QN("Number",0)=>2,QN("Number",-1)=>1,QN("Number",1)=>1]
    else
        error("Invalid operator string")
    end
end

function make_sites_list(what_given::Vector{String},wavefunc::TTNKit.TreeTensorNetwork; kwargs...)
    phys_sites::Vector{Vector{Index}} = [Vector{Index}(undef,TTNKit.number_of_sites(wavefunc.net)) for i in 1:length(what_given)]
    for i in 1:length(what_given)
        if what_given[i] == "pull"
            phys_sites[i] = TTNKit.sites(wavefunc)
        elseif what_given[i] == "new"
            phys_sites[i] = [Index(make_qnset("Adag"); tags="MPO-Internal,Site,n=$i") for i in 1:TTNKit.number_of_sites(wavefunc.net)]
        elseif what_given[i] == "provided"
            phys_sites[i] = kwargs[:provided_sites][i]
        else
            error("Invalid option")
        end
    end
    return phys_sites
end

function projected_op_mpo(wavefunc::TTNKit.TreeTensorNetwork,op_type::String; kwargs...)
    lat = TTNKit.physical_lattice(wavefunc.net)
    mapping = get(kwargs,:mapping,collect(1:TTNKit.number_of_sites(lat)))#TTNKit.hilbert_curve(lat)
    if_wrap::Bool = get(kwargs,:if_wrap,true)
    what_given_inds::Vector{String} = get(kwargs,:what_given_inds,["pull","pull"])

    top_phys_sites,bottom_phys_sites = make_sites_list(what_given_inds,wavefunc; kwargs...)

    tensor_train = Vector{ITensor}(undef,length(top_phys_sites))
    all_indices = Vector{Index}(undef,length(top_phys_sites)+1)

    for s in 0:length(top_phys_sites)
        qnset = make_qnset(op_type)
        if s == 0
            left_tag = "Start"
            right_tag = string(s+1)
        elseif s == length(top_phys_sites)
            left_tag = string(s)
            right_tag = "End"
        else
            left_tag = string(s)
            right_tag = string(s+1)
        end
        all_indices[s+1] = Index(qnset; tags="Link,Left=$left_tag,Right=$right_tag")
    end


    for (idx,s) in enumerate(top_phys_sites)
        println("Working on Physical Site $(TTNKit.tags(s))")

        coeff::ComplexF64 = 1.0 + 0.0*im

        mat = build_W(op_type,coeff)
        local_inds = [all_indices[idx],dag(bottom_phys_sites[idx]),dag(all_indices[idx+1]),prime(top_phys_sites[idx])]
        mit = ITensor(mat,local_inds)

        if idx == 1
            mit = mit * dag(onehot(all_indices[1] => 1))
        elseif idx == length(top_phys_sites)
            mit = mit * onehot(all_indices[end] => 2)
        end

        tensor_train[idx] = mit

    end

    if if_wrap
        println("Mapping only linear")
        return make_mpowrapper(MPO(tensor_train),lat)
    else
        return MPO(tensor_train)
    end

end

function find_physical_inflow_sites(ind_list)
    all_inflowing_sites = Vector{Index}(undef,length(ind_list))
    for (idx,s) in enumerate(ind_list)
        # remove the virtual indices and find the inflowing physical site
        correct_index = filter(x -> occursin("Site",string(tags(x))) && string(x.dir)=="In",s)[1]
        all_inflowing_sites[idx] = correct_index
    end
    return all_inflowing_sites
end

function combine_mpos(m1::MPO,m2::MPO; kwargs...)
    tensor_train = Vector{ITensor}(undef,length(m1))

    for i in 1:length(m1)
        println("Combine MPOs on site $(i)")
        ti = m1[i] * m2[i]

        #=all_inds = inds(ti)
        for i1 in 1:length(all_inds)
            in1 = all_inds[i1]
            for i2 in i1+1:length(all_inds)
                in2 = all_inds[i2]
                #println("Look at indices $(in1) and $(in2)")
                if occursin("Left",string(tags(in1))) && string(tags(in1)) == string(tags(in2)) #&& true # include that the flow is in the same direction as check
                    #println("Found that the indices are the same: $(tags(in1)) and $(tags(in2))")
                    ci = combiner(in1,in2)
                    ti = ti * ci
                end
            end
        end=#

        tensor_train[i] = ti
    end

    return MPO(tensor_train)
end

function two_point_mpo(wavefunc::TTNKit.TreeTensorNetwork; kwargs...)
    if_wrap::Bool = get(kwargs,:if_wrap,true)

    annih = projected_op_mpo(wavefunc,"A"; if_wrap=false,what_given_inds=["pull","new"])

    inflow_inds = dag.(find_physical_inflow_sites(inds.(annih)))

    creat = projected_op_mpo(wavefunc,"Adag"; if_wrap=false,what_given_inds=["provided","pull"],provided_sites=(inflow_inds,nothing))

    rho = combine_mpos(annih,noprime(creat))
    println("Combined MPOs")
    if if_wrap
        lat = TTNKit.physical_lattice(wavefunc.net)
        return make_mpowrapper(rho,lat)
    else
        return rho
    end
end

function two_point_mpo_reverse(wavefunc::TTNKit.TreeTensorNetwork; kwargs...)
    if_wrap::Bool = get(kwargs,:if_wrap,true)

    creat = projected_op_mpo(wavefunc,"Adag"; if_wrap=false,what_given_inds=["pull","new"])

    inflow_inds = dag.(find_physical_inflow_sites(inds.(creat)))

    annih = projected_op_mpo(wavefunc,"A"; if_wrap=false,what_given_inds=["provided","pull"],provided_sites=(inflow_inds,nothing))

    rho = combine_mpos(noprime(annih),creat)
    println("Combined MPOs")
    if if_wrap
        lat = TTNKit.physical_lattice(wavefunc.net)
        return make_mpowrapper(rho,lat)
    else
        return rho
    end
end







function build_W_singlepoint(op_type::String,coeff::ComplexF64)
    if op_type == "A"
        mat = [0.0 0.0; 1.0 0.0]
    elseif op_type == "Adag"
        mat = [0.0 1.0; 0.0 0.0]
    elseif op_type == "N"
        mat = [1.0 0.0; 0.0 0.0]
    else
        error("Invalid operator type")
    end

    A11 = I(2)
    A12 = zeros(2,2)
    A21 = coeff * mat
    A22 = I(2)
    M = Array{Float64, 4}(undef, 2,2,2,2)
    M[1,:,1,:] = A11
    M[1,:,2,:] = A21
    M[2,:,1,:] = A12
    M[2,:,2,:] = A22

    return M
end

function build_links_singlepoint(op_type::String,L::Int)
    links = Vector{Index}(undef,L+1)
    for i in 0:L
        qnset = make_qnset(op_type)
        if i == 0
            left_tag = "Start"
            right_tag = string(i+1)
        elseif i == L
            left_tag = string(i)
            right_tag = "End"
        else
            left_tag = string(i)
            right_tag = string(i+1)
        end
        links[i+1] = Index(qnset; tags="Link,Left=$left_tag,Right=$right_tag")
    end
    return links
end

function single_point_mpo(wavefunc::TTNKit.TreeTensorNetwork,op_type::String; kwargs...)
    if_wrap::Bool = get(kwargs,:if_wrap,true)

    lat = TTNKit.physical_lattice(wavefunc.net)

    phys_sites = TTNKit.sites(wavefunc)

    links = build_links_singlepoint(op_type,length(phys_sites))

    tensor_train = Vector{ITensor}(undef,length(phys_sites))
    for (idx,s) in enumerate(phys_sites)
        println("Working on Physical Site $(TTNKit.tags(s))")

        coeff::ComplexF64 = 1.0 + 0.0*im

        mat = build_W_singlepoint(op_type,coeff)

        local_inds = [links[idx],dag(s),dag(links[idx+1]),prime(s)]

        mit = ITensor(mat,local_inds)

        if idx == 1
			mit = mit * dag(onehot(links[1] => 1))
		elseif idx == length(phys_sites)
			mit = mit * (onehot(links[end] => 2))
		end
        
        tensor_train[idx] = mit
    end

    if if_wrap
        return make_mpowrapper(MPO(tensor_train),lat)
    else
        return MPO(tensor_train)
    end
end

function reshape_mpo_to_matrix(mpo::MPO)
    s1 = [siteinds(mpo)[i][1] for i in 1:4]
    s2 = [siteinds(mpo)[i][2] for i in 1:4]
    cm1 = combiner(s1)
    cm2 = combiner(s2)
    
    rez = (contract(mpo) * cm1) * cm2
    return Matrix(rez,inds(rez)[1],inds(rez)[2])
end











































"fin"
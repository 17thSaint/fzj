#####################################################
#=

This file is for converting the old Backend TTN data type to the main_new version without the Backend variable

    Need to convert:
        TreeTensorNetwork
        Lattice
        Network
        TPO
        MPO

=#
######################################################

include("data-storage-funcs.jl")

# convert ITensorNode to regular Node
function convert_ITensorNode_to_Node(itn::TTNKit.AbstractNode)
    s = itn.s
    desc = itn.desc
    return TTNKit.Node{Index,typeof(s)}(s,desc)
end

# get raw data for SimpleLattice
function rawdata_SimpleLattice(dataname::String,old_lattice::TTNKit.SimpleLattice)
    rawdata_dict = Dict()
    lat_vec = [convert_ITensorNode_to_Node(n) for n in old_lattice.lat]
    rawdata_dict[dataname*"-lat"] = lat_vec
    rawdata_dict[dataname*"-dims"] = old_lattice.dims

    return rawdata_dict
end

# rebuild the new SimpleLattice from raw data
function rebuild_SimpleLattice(data_dict::Dict)
    lat = nothing
    dims = nothing
    for k in collect(keys(data_dict))
        if "lat" == split(k,"-")[end]
            lat = data_dict[k]
        elseif "dims" == split(k,"-")[end]
            dims = data_dict[k]
        end
    end
    return TTNKit.SimpleLattice{length(dims), TTNKit.spacetype(lat[1]), TTNKit.sectortype(lat[1])}(lat, dims)
end

# get raw data for BinaryNetwork
function rawdata_BinaryNetwork(dataname::String,old_network::TTNKit.BinaryNetwork)
    rawdata_dict = Dict([(dataname*"-lattices",[])])
    for (idx,l) in enumerate(old_network.lattices)
        append!(rawdata_dict[dataname*"-lattices"],[rawdata_SimpleLattice(dataname*"-lattice$idx",l)])
    end

    return rawdata_dict
end

# rebuild the new BinaryNetwork from raw data
function rebuild_BinaryNetwork(data_dict::Dict)
    lattices = []
    for v in data_dict["net-lattices"]
        append!(lattices,[rebuild_SimpleLattice(v)])
    end
    return TTNKit.BinaryNetwork{typeof(lattices[1])}(lattices)
end

# get raw data for TreeTensorNetwork
function rawdata_TreeTensorNetwork(dataname::String,old_ttn::TTNKit.TreeTensorNetwork)
    rawdata_dict = Dict()
    rawdata_dict[dataname*"-data"] = old_ttn.data
    rawdata_dict[dataname*"-ortho_direction"] = old_ttn.ortho_direction
    rawdata_dict[dataname*"-ortho_center"] = old_ttn.ortho_center
    rawdata_dict[dataname*"-net"] = rawdata_BinaryNetwork(dataname*"-net",old_ttn.net)

    return rawdata_dict
end

# rebuild the new TreeTensorNetwork from raw data
function rebuild_TreeTensorNetwork(data_dict::Dict)
    data = data_dict["ttn-data"]
    ortho_direction = data_dict["ttn-ortho_direction"]
    ortho_center = data_dict["ttn-ortho_center"]
    net = rebuild_BinaryNetwork(data_dict["ttn-net"])
    return TTNKit.TreeTensorNetwork{typeof(net),typeof(data[1][1])}(data, ortho_direction, ortho_center, net)
end

# save raw data and kill the original data
function save_rawdata(rawdata,dataname::String,filepath::String,which_group::String="metadata")
    data_dict = Dict([(dataname*"-rawdata",rawdata),(dataname,nothing)])
    modify_data_jld2(data_dict,filepath,which_group; output_level=2)
end

# save new data and kill the raw data
function save_newdata(newdata,dataname::String,filepath::String,which_group::String="metadata")
    data_dict = Dict([(dataname,newdata),(dataname*"-rawdata",nothing)])
    modify_data_jld2(data_dict,filepath,which_group; output_level=2)
end

#filename = "ttn-if_periodic_phys-true-onsite_strength-10.0-lr-7-particles-4-alpha-0.125-layers-6-hopping_anisotropy-1.0.jld2"
#d,m = read_data_jld2(filename; output_level=0)
#bf = jldopen(filename,"r+")

# save the raw data to the original file
if false
    old_net = m["net"]
    rawdata_net = rawdata_BinaryNetwork("net",old_net)
    save_rawdata(rawdata_net,"net",filename)

    old_ttn = d["ttn"]
    rawdata_ttn = rawdata_TreeTensorNetwork("ttn",old_ttn)
    save_rawdata(rawdata_ttn,"ttn",filename,"all_data")
end

# extract rawdata from file and then build the new network and save it
if false
    rawdata_net = bf["metadata"]["net-rawdata"]
    close(bf)
    new_net = rebuild_BinaryNetwork(rawdata_net)
    save_newdata(new_net,"net",filename)
end









































"fin"
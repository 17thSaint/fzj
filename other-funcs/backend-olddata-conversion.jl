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
    return TTNKit.Node{TTNKit.Index,typeof(s)}(s,desc)
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
    return TTNKit.SimpleLattice{length(dims),TTNKit.spacetype(lat[1]),TTNKit.sectortype(lat[1])}(lat, dims)
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
    all_keys = collect(keys(data_dict))
    for v in data_dict[all_keys[findfirst(x -> occursin("net-lattices",x), all_keys)]]
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
    return TTNKit.TreeTensorNetwork(data, ortho_direction, ortho_center, net)
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

#=filename = "ttn-if_periodic_phys-true-onsite_strength-10.0-lr-7-particles-4-alpha-0.125-layers-6-hopping_anisotropy-1.0.jld2"
bf = jldopen(filename,"r+")

# save the raw data to the original file
if false
    old_net = bf["metadata"]["net"]
    rawdata_net = rawdata_BinaryNetwork("net",old_net)
    
    old_ttn = bf["all_data"]["ttn"]
    rawdata_ttn = rawdata_TreeTensorNetwork("ttn",old_ttn)

    close(bf)

    save_rawdata(rawdata_net,"net",filename,"metadata")
    save_rawdata(rawdata_ttn,"ttn",filename,"all_data")
end

# extract rawdata from file and then build the new network and save it
if true
    rawdata_net = bf["metadata"]["net-rawdata"]
    new_net = rebuild_BinaryNetwork(rawdata_net)

    rawdata_ttn = bf["all_data"]["ttn-rawdata"]
    new_ttn = rebuild_TreeTensorNetwork(rawdata_ttn)

    close(bf)

    save_newdata(new_net,"net",filename)
    save_newdata(new_ttn,"ttn",filename,"all_data")
end=#

function extract_data(binary_file,filename::String)
    
    # extract all keys from the metadata group of the binary file
    metadata_keys = collect(keys(binary_file["metadata"]))
    all_raw_metadata = Dict()
    for k in metadata_keys

        # look for the string TTN within each key
        if occursin("ttn",k)
            local_old_ttn = binary_file["metadata"][k]
            # check if the data value is empty
            if !isnothing(local_old_ttn)
                # extract rawdata and add to rawdata dictionary
                local_rawdata_ttn = rawdata_TreeTensorNetwork(k,local_old_ttn)
                all_raw_metadata[k] = local_rawdata_ttn
            end
        # look for the string net within each key
        elseif occursin("net",k)
            local_old_net = binary_file["metadata"][k]
            # check if the data value is empty
            if !isnothing(local_old_net)
                # extract rawdata and add to rawdata dictionary
                local_rawdata_net = rawdata_BinaryNetwork(k,local_old_net)
                all_raw_metadata[k] = local_rawdata_net
            end
        end
    end

    # extract keys from the all_data group of binary file, will only include TTN, never just net
    all_data_keys = collect(keys(binary_file["all_data"]))
    all_raw_alldata = Dict()
    for k in all_data_keys
        # look for the string TTN within each key
        if occursin("ttn",k)
            local_old_ttn = binary_file["all_data"][k]
            # check if the data value is empty
            if !isnothing(local_old_ttn)
                # extract rawdata and add to rawdata dictionary
                local_rawdata_ttn = rawdata_TreeTensorNetwork(k,local_old_ttn)
                all_raw_alldata[k] = local_rawdata_ttn
            end
        end
    end

    # close binary file before saving rawdata
    close(binary_file)

    # run through rawdata dictionaries and save them to the file
    for (k,v) in all_raw_metadata
        save_rawdata(v,k,filename,"metadata")
    end

    for (k,v) in all_raw_alldata
        save_rawdata(v,k,filename,"all_data")
    end

end

function rebuild_data(binary_file,filename::String)
    # extract all keys from the metadata group of the binary file
    metadata_keys = collect(keys(binary_file["metadata"]))
    all_new_metadata = Dict()
    for k in metadata_keys

        # look for the string rawdata within each key
        if occursin("rawdata",k)
            local_rawdata = binary_file["metadata"][k]
            # check if the data value is empty
            if !isnothing(local_rawdata)
                if occursin("net",k)
                    # rebuild new NET data and add to new data dictionary
                    local_newdata = rebuild_BinaryNetwork(local_rawdata)
                    all_new_metadata[join(split(k,"-")[1:end-1],"-")] = local_newdata
                elseif occursin("ttn",k)
                    # rebuild new TTN data and add to new data dictionary
                    local_newdata = rebuild_TreeTensorNetwork(local_rawdata)
                    all_new_metadata[join(split(k,"-")[1:end-1],"-")] = local_newdata
                end
            end
        end
        
    end

    # extract keys from the all_data group of binary file, will only include TTN, never just net
    all_data_keys = collect(keys(binary_file["all_data"]))
    all_new_alldata = Dict()
    for k in all_data_keys
        # look for the string rawdata within each key
        if occursin("rawdata",k)
            local_rawdata = binary_file["all_data"][k]
            # check if the data value is empty
            if !isnothing(local_rawdata)
                if occursin("ttn",k)
                    # rebuild new TTN data and add to new data dictionary
                    local_newdata = rebuild_TreeTensorNetwork(local_rawdata)
                    all_new_alldata[join(split(k,"-")[1:end-1],"-")] = local_newdata
                end
            end
        end
    end
   

    # close binary file before saving rawdata
    close(binary_file)

    println("Rebuilding metadata: ")
    display(collect(keys(all_new_metadata)))
    # run through rawdata dictionaries and save them to the file
    for (k,v) in all_new_metadata
        save_newdata(v,k,filename,"metadata")
    end

    println("Rebuilding all_data: ")
    display(collect(keys(all_new_alldata)))
    for (k,v) in all_new_alldata
        save_newdata(v,k,filename,"all_data")
    end

end

# need to get rawdata for DefaultExpander and Observers
if false
    all_files = readdir()
    filter!(x -> occursin(".jld2",x),all_files)
    filter!(x -> occursin("ttn",x),all_files)

    display(all_files)

    for f in all_files
        println("Starting to extract $f")
        bf = jldopen(f,"r+")
        extract_data(bf,f)
    end
end

# not rebuilding ttn from all_data correctly, maybe the TypeError
if true
    all_files = readdir()
    filter!(x -> occursin(".jld2",x),all_files)
    filter!(x -> occursin("ttn",x),all_files)

    display(all_files)

    for f in all_files
        println("Starting to rebuild $f")
        bf = jldopen(f,"r+")
        try
           rebuild_data(bf,f)
        catch err
            if err isa TypeError
                println("Found TypeError but continuing")
            else
                rethrow(err)
            end
        end
    end
end









































"fin"
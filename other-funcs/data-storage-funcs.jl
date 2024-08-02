using JLD2

function get_folder_location(folder_name)
	central_loc=find_center()
	get_to_center = split(pwd(),central_loc)[1]
	full_location = occursin("main-git",pwd()) ? get_to_center * "fzj/main-git/" * folder_name : get_to_center * central_loc * "/" * folder_name
	return full_location
end

function find_center()
	all_folders = split(pwd(),"/")
	if "fzj" in all_folders
		return "fzj"
	elseif "local" in all_folders
		return all_folders[findfirst(x -> all_folders[x] == "local",1:length(all_folders))+1]
	elseif "Local" in all_folders
		return all_folders[findfirst(x -> all_folders[x] == "Local",1:length(all_folders))+1]
	else
		println("Not sure where the center is: $(pwd())")
	end
end

function check_cluster()
	if occursin("local",pwd()) || occursin("Local",pwd())
		return true
	else
		return false
	end
end

function named_tuple_to_dict(namedtuple)
	new_dict = Dict{String,Any}()
	for key in keys(namedtuple)
		new_dict[string(key)] = namedtuple[key]
	end
	return new_dict
end

function check_datafolder_exists(folder_name::String)
	current_location = pwd()
	try
		cd(folder_name)
		cd(current_location)
		correct_location = folder_name
		return correct_location
	catch
		println("The directory $folder_name is not accessible")
		correct_location = get_folder_location("cluster-data"*split(folder_name,"cluster-data")[end])
		return correct_location
	end
end

function dict_to_symbols(given_dict)
	new_dict = Dict()
	for (key,value) in given_dict
		new_dict[Symbol(key)] = value	
	end	
	return new_dict
end	

function turn_string_into_bool(input)
	if typeof(input) != Bool
		if input == "true"
			input = true
		elseif input == "false"
			input = false
		else
			println("Input is not true or false, Input = $input")
		end
	end
	return input
end

function make_args_dict(args,parameter_iteration=0)
	parameters_dict = Dict()
	if length(args) == 1
		args = collect(split(args[1],","))
	end
	if length(args) % 2 != 0
		parameter_shift = parse(Float64,args[3])*parameter_iteration
		deleteat!(args,3)
	else
		parameter_shift = 0.0
	end
	for i in 1:2:length(args)
		key = args[i]
		value = args[i+1]
		if key == "dataloc"
			parameters_dict["dataloc"] = value
			continue
		end
		get_integer = tryparse(Int,value) 
		if isnothing(get_integer)
			get_float = tryparse(Float64,value)
			if isnothing(get_float)
				value = turn_string_into_bool(value)
			else
				value = get_float
			end
		else
			value = get_integer
		end
		parameters_dict[key] = value#i == 1 ? value + parameter_shift : value
	end
	return parameters_dict
end

function make_parameters_filename(param_dict)
	param_filename = ""
	for key in keys(param_dict)
		if key in ["change", "if_change", "open_cores"]
			continue
		else
			if typeof(param_dict[key]) == Float64 && param_dict[key] < 0
				value = "n" * string(abs(param_dict[key]))
			else
				value = string(param_dict[key])
			end
			param_filename *= "$key-$value-"
		end
	end
	param_filename = chop(param_filename,tail=1)
	return param_filename
end

function get_params_dict_from_filename(filename)
	params_dict = Dict()
	
	file_type = split(filename,".")[end]
	if file_type == "png" || file_type == "jld2" || file_type == "hdf5"
		filename = join(split(filename,".")[1:end-1],".")
	end
	if split(filename,"-")[1] in ["virt","phys","Y","X"]
		split_filename = split(filename,"-")[4:end]
	elseif split(filename,"-")[1] in ["mps","ttn","rfa","laugh","basis","ed","wavefuncmps","wavefuncttn"]
		split_filename = split(filename,"-")[2:end]
	else
		split_filename = split(filename,"-")
	end

	for i in 1:Int(length(split_filename)/2)
		key = split_filename[2*i-1]
		value = split_filename[2*i]
		if value[1] == 'n'
			value = "-" * string(value[2:end])
		end
		get_integer = tryparse(Int,value) 
		if isnothing(get_integer)
			get_float = tryparse(Float64,value)
			if isnothing(get_float)
				value = turn_string_into_bool(value)
			else
				value = get_float
			end
		else
			value = get_integer
		end
		params_dict[key] = value
	end
	return params_dict
end

function rewrite_filename(old_name,new_params_dict)
	old_name_dict = get_params_dict_from_filename(old_name)
	params_to_keep = Dict()
	for (k,v) in new_params_dict
		if k in keys(old_name_dict)
			params_to_keep[k] = v
		end
	end
	
	if "phi" in keys(new_params_dict) && "alpha" in keys(old_name_dict)
		params_to_keep["alpha"] = round(new_params_dict["phi"] / (2*pi),digits=4)
	end
	new_name_dict = merge(old_name_dict,params_to_keep)
	new_name = make_parameters_filename(new_name_dict)
	return new_name
end

function separate_filename_location(filename::String)
	if occursin("/",filename)
		fullsplit = split(filename,"/")
		location = join(fullsplit[1:end-1],"/")
		actual_filename = fullsplit[end]
		return location,actual_filename
	else
		println("No directory splits in input: $filename")
		return nothing,filename
	end
end

function change_numparticles_metadata(filename)
	binary_file = jldopen("../cluster-data/$filename","a+")
	try
		metadata = binary_file["metadata"]
		metadata["particles"] = metadata["num_particles"]
		delete!(metadata,"num_particles")
		close(binary_file)
		println("Completed for $filename")
		return "yes"
	catch
		close(binary_file)
		println("Had some error for $filename")
		return "no"
	end
end

function find_data_file(params_dict,calc_type,location="/home/patrick/fzj/main-git/cluster-data"; kwargs...)
	file_type = get(kwargs, :file_type, "jld2")
	og_loc = pwd()
	cd(location)
	
	file_choices = readdir()
	
	remove_indices = []
	for i in 1:length(file_choices)
		split(file_choices[i],".")[end] != file_type ? append!(remove_indices,[i]) : nothing
	end
	deleteat!(file_choices,remove_indices)
	remove_indices = []

	for i in 1:length(file_choices)
		current_calc_type = split(file_choices[i],"-")[1]
		current_calc_type in ["virt","phys","Y","X"] ? current_calc_type = join(split(file_choices[i],"-")[1:3],"-") : nothing
		current_calc_type != calc_type ? append!(remove_indices,[i]) : nothing
	end
	deleteat!(file_choices,remove_indices)

	remove_indices = []
	for i in 1:length(file_choices)
		current_file = file_choices[i]
		current_file_params_dict = get_params_dict_from_filename(current_file)
		for params in keys(params_dict)
			try
				current_file_params_dict[params] in params_dict[params] ? nothing : append!(remove_indices,i)
			catch
				try
					#println("Error in first Attempt, $current_file")
					current_file_metadata_dict = read_data_jld2(current_file; kwargs...)[2]
					current_file_metadata_dict[params] in params_dict[params] ? nothing : append!(remove_indices,i)
				catch
					#println("Parameter $params could not be found in file $current_file, skipping")
					append!(remove_indices,i)
				end
			end
		end
	end
	deleteat!(file_choices,unique(remove_indices))
	
	cd(og_loc)
	return file_choices
end

function make_sure_file_type(file_name,desired_type)
	split_name = split(file_name,".")
	potential_type = split_name[end]
	if potential_type in ["png","jld2","hdf5","txt"]
		if potential_type != desired_type
			println("Wrong File Type: changing $potential_type => $desired_type")
			file_name = join(split_name[1:end-1],".") * ".$desired_type"
		end
	else
		#println("No File Type: adding $desired_type")
		file_name *= "." * desired_type
	end
	return file_name
end

function check_plot_label(file_name,version)
	split_name = split(file_name,"-")
	potential_type = split_name[1]
	if potential_type in ["densdens","occs","ttn","Y-dir-GF","X-dir-GF","Y-dir-current","X-dir-current","mps","rfa","laugh","tevo"]
		if potential_type != version
			println("Wrong Plot Label: changing $potential_type => $version")
			file_name = "$version-" * join(split_name[2:end],"-")
		end
	else
		println("No File Type: adding $version")
		file_name = "$version-" * file_name
	end
	return file_name
end

function check_data_exists(params_dict::Dict,data_type::String; kwargs...)
	location = get(kwargs, :location, "/home/patrick/fzj/main-git/cluster-data")
	possible_files = find_data_file(params_dict,data_type,location; kwargs...)
	if length(possible_files) == 0
		return false,nothing
	else
		if length(possible_files) > 1
			println("Too Many files")
			return false,nothing
		else
			return true,read_data_jld2(possible_files[1],location; kwargs...)
		end
	end
end
	
function check_data_exists(filename::AbstractString,data_type="observer"; kwargs...)
    location = get(kwargs, :location, pwd())
    here = pwd()
    cd(location)
    if check_duplicates(filename)[1]
        cd(here)
        return true,read_data_jld2(filename,location; kwargs...)
    end
    cd(here)
    return false,nothing
end

function check_duplicates(file_name)
	split_period_name = split(file_name,".")
	file_type = split_period_name[end]
	file_string = join(split_period_name[1:end-1],".")
	rename = false
	if file_name in readdir()
		rename = true
	end
	while file_name in readdir()
		string_elems = split(file_string,"-")
		if "mk" in string_elems
			mk_number_loc = findfirst(x -> x == "mk",string_elems) + 1
			string_elems[mk_number_loc] = string(parse(Int,string_elems[mk_number_loc]) + 1)
		else
			append!(string_elems,["mk","2"])
		end
		file_string = join(string_elems,"-")
		file_name = file_string * "." * file_type
	end
	if rename
		println("Found Duplicate File, renaming $file_name")
	end
	return rename,file_name
end

function check_duplicates(filename,location)
	og_loc = pwd()
	cd(location)
	result = check_duplicates(filename)
	cd(og_loc)
	return result
end

function prep_file(file_name,desired_type)
	file_name = make_sure_file_type(file_name,desired_type)
	file_name = check_duplicates(file_name)[2]
	return file_name
end

function save_figure(file_name; kwargs...)
	file_location = get(kwargs, :location, pwd())
	file_type = get(kwargs, :file_type, "png")
	current_location = pwd()
	cd(file_location)
	file_name = prep_file(file_name,file_type)
	PyPlot.savefig(file_name)
	cd(current_location)
	println("Figure Saved, File Closed: $file_name")
	return
end

function check_dict(data)
	if isa(data,Dict)
		return data
	elseif isa(data,NamedTuple)
		println("Data was NamedTuple, converting to Dictionary")
		return named_tuple_to_dict(data)
	else
		println("Data being sent isn't Dictionary or NamedTuple")
		return
	end
end

function remove_wavefunc_from_filename(filename::String)
	split_filename = split(filename,"-")
	if occursin("wavefunc",split_filename[1])
		datatype = string(split(split_filename[1],"wavefunc")[2])
		filename = join(vcat([datatype],split_filename[2:end]),"-")
	end
	return filename
end

function add_wavefunc_to_filepath(filepath::String)
	split_filepath = split(filepath,"/")
	filename = split_filepath[end]
	return join(vcat(split_filepath[1:end-1],["wavefunc$filename"]),"/")
end

function write_data_jld2(file_name::AbstractString,data::Dict,location=pwd(),metadata=nothing; kwargs...)
	og_location = pwd()
	try
		cd(location)
	catch
		println("The directory $location is not accessible")
		location = pwd()
	end

	file_name = prep_file(file_name,"jld2")
	data = check_dict(data)
	
	binary_file = jldopen(file_name,"w")
	
	if !isnothing(metadata)
		metadata = check_dict(metadata)
		metadata_var = JLD2.Group(binary_file,"metadata")
		for metadatum_key in keys(metadata)
			metadata_var[metadatum_key] = metadata[metadatum_key]
		end
	end
	
	alldata = JLD2.Group(binary_file,"all_data")
	for datum_key in keys(data)
		alldata[datum_key] = data[datum_key]
	end
	
	close(binary_file)
	cd(og_location)
	println("Data Added, File Closed: $file_name")
	return file_name
end

function write_data_jld2(location::AbstractString,data::Dict,metadata::Dict; kwargs...)
	split_location = split(location,"/")
	file_name = split_location[end]
	location = join(split_location[1:end-1],"/")
	write_data_jld2(file_name,data,location,metadata; kwargs...)
end

function read_data_jld2(file_name,location=pwd(); kwargs...)
	output_bool = get(kwargs, :output_level, 1)

	if occursin("/",file_name)
		location = join(split(file_name,"/")[1:end-1],"/")
		file_name = split(file_name,"/")[end]
	end

	og_location = pwd()
	try
		cd(location)
	catch
		try
			cd(get_folder_location(location))
		catch
			println("The directory $location is not accessible")
			location = pwd()
		end
	end
	
	file_name = make_sure_file_type(file_name,"jld2")
	
	output = []
	binary_file = jldopen(file_name,"r")
	
	data = Dict()
	alldata = binary_file["all_data"]
	for datum_key in keys(alldata)
		data[datum_key] = alldata[datum_key]
	end
	push!(output,data)
	
	if !isnothing(findfirst(x -> x == "metadata",keys(binary_file)))
		metadata = Dict()
		metadata_vals = binary_file["metadata"]
		for metadatum_key in keys(metadata_vals)
			metadata[metadatum_key] = metadata_vals[metadatum_key]
		end
	push!(output,metadata)
	end
	
	close(binary_file)
	cd(og_location)
	if output_bool > 0
		println("Data Extracted, File Closed: $file_name")
	end
	
	if length(output) < 2
		return data
	else
		return output
	end
end

function modify_data_jld2(key_to_modify::String, new_value, file_path, which_group="all_data"; kwargs...)
	output_level = get(kwargs, :output_level, 0)
	
	file_name = split(file_path,"/")[end]
	file_name = make_sure_file_type(file_name,"jld2")
	file_path = join(split(file_path,"/")[1:end-1],"/") * "/" * file_name

	# Open the JLD2 file in write mode
    jld_file = jldopen(file_path, "a+")
	data_dict = jld_file[which_group]

    # Check if the key exists in the file
    if haskey(data_dict, key_to_modify)
        # Modify the value associated with the given key
		delete!(data_dict, key_to_modify)
        data_dict[key_to_modify] = new_value
        output_level > 1 ? println("Value associated with key $key_to_modify has been modified.") : nothing
    else
        output_level > 1 ? println("Key $key_to_modify does not exist in the file. Making it") : nothing
		data_dict[key_to_modify] = new_value
    end
	#
    
    # Close the JLD2 file
    close(jld_file)
	return split(file_path,"/")[end]
end

function modify_data_jld2(to_modify_dict::Dict,file_path, which_group="all_data"; kwargs...)
	output_level = get(kwargs, :output_level, 0)

	file_name = split(file_path,"/")[end]
	file_name = make_sure_file_type(file_name,"jld2")
	file_path = join(split(file_path,"/")[1:end-1],"/") * "/" * file_name

	# Open the JLD2 file in write mode
	jld_file = jldopen(file_path, "a+")
	data_dict = jld_file[which_group]

	# Check if the key exists in the file
	for (key_to_modify,new_value) in to_modify_dict
		if haskey(data_dict, key_to_modify)
			# Modify the value associated with the given key
			delete!(data_dict, key_to_modify)
			data_dict[key_to_modify] = new_value
			output_level > 1 ? println("Value associated with key $key_to_modify has been modified.") : nothing
		else
			output_level > 1 ? println("Key $key_to_modify does not exist in the file. Making it") : nothing
			data_dict[key_to_modify] = new_value
		end
	end
	#

	# Close the JLD2 file
	close(jld_file)
	return split(file_path,"/")[end]
end

#=function read_data_hdf5(file_name,location=pwd(); kwargs...)
	output_level = get(kwargs, :output_level, true)
	og_location = pwd()
	cd(location)

	file_name = make_sure_file_type(file_name,"hdf5")

	output = []
	binary_file = h5open(file_name,"r")

	data = Dict()
	alldata = binary_file["all_data"]
	for datum_key in keys(alldata)
		data[datum_key] = read(alldata,datum_key)
	end
	push!(output,data)

	if !isnothing(findfirst(x -> x == "metadata",keys(binary_file)))
		metadata = Dict()
		metadata_vals = binary_file["metadata"]
		for metadatum_key in keys(metadata_vals)
			metadata[metadatum_key] = read(metadata_vals,metadatum_key)
		end
	push!(output,metadata)
	end

	close(binary_file)
	cd(og_location)
	output_level ? println("Data Extracted, File Closed: $file_name") : nothing

	if length(output) < 2
		return data
	else
		return output
	end
end=#





















"fin"

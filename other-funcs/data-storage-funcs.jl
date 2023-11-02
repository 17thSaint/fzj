using JLD2,TTNKit

function back2cpu(ttn::TTNKit.TreeTensorNetwork)
	datagpu = deepcopy(ttn.data)
	datac = map(datagpu) do layerdata
		map(T -> TTNKit.cpu(T), layerdata)
	end
	ortho_centerc = deepcopy(ttn.ortho_center)
	netc = deepcopy(ttn.net)
	ortho_directionc = deepcopy(ttn.ortho_direction)
	return TreeTensorNetwork(datac, ortho_directionc, ortho_centerc, netc)
end

function named_tuple_to_dict(namedtuple)
	new_dict = Dict{String,Any}()
	for key in keys(namedtuple)
		new_dict[string(key)] = namedtuple[key]
	end
	return new_dict
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
		parameters_dict[key] = i == 1 ? value + parameter_shift : value
	end
	return parameters_dict
end

function make_parameters_filename(param_dict)
	param_filename = ""
	for key in keys(param_dict)
		if key in ["change", "if_change", "open_cores"]
			continue
		else
			value = string(param_dict[key])
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
	elseif split(filename,"-")[1] in ["mps","ttn","rfa"]
		split_filename = split(filename,"-")[2:end]
	else
		split_filename = split(filename,"-")
	end
	for i in 1:Int(length(split_filename)/2)
		key = split_filename[2*i-1]
		value = split_filename[2*i]
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

function get_folder_location(folder_name)
	get_to_fzj = split(pwd(),"fzj")[1]
	full_location = occursin("main-git",pwd()) ? get_to_fzj * "fzj/main-git/" * folder_name : get_to_fzj * "fzj/" * folder_name
	return full_location
end

function find_data_file(params_dict,calc_type,file_type="jld2",location="/home/patrick/fzj/main-git/cluster-data")
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
					current_file_metadata_dict = read_data_jld2(current_file,; output_bool=false)[2]
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
		println("No File Type: adding $desired_type")
		file_name *= "." * desired_type
	end
	return file_name
end

function check_plot_label(file_name,version)
	split_name = split(file_name,"-")
	potential_type = split_name[1]
	if potential_type in ["densdens","occs","ttn","Y-dir-GF","X-dir-GF","Y-dir-current","X-dir-current","mps","rfa"]
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
	return file_name
end

function prep_file(file_name,desired_type)
	file_name = make_sure_file_type(file_name,desired_type)
	file_name = check_duplicates(file_name)
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

function write_data_jld2(file_name,data,location=pwd(),metadata=nothing; kwargs...)
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
	return
end

function read_data_jld2(file_name,location=pwd(); kwargs...)
	output_bool = get(kwargs, :output_bool, true)
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
	if output_bool
		println("Data Extracted, File Closed: $file_name")
	end
	
	if length(output) < 2
		return data
	else
		return output
	end
end

function read_data_hdf5(file_name,location=pwd(); kwargs...)
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
	println("Data Extracted, File Closed: $file_name")

	if length(output) < 2
		return data
	else
		return output
	end
end





















"fin"

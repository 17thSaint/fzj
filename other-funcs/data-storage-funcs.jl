using JLD2

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

function find_data_file(params_dict,calc_type,file_type="jld2")
	og_loc = pwd()
	cd("/home/patrick/fzj/main-git/cluster-data")
	
	file_choices = readdir()
	
	remove_indices = []
	for i in 1:length(file_choices)
		split(file_choices[i],".")[end] != file_type ? append!(remove_indices,[i]) : nothing
	end
	deleteat!(file_choices,remove_indices)
	remove_indices = []

	for i in 1:length(file_choices)
		split(file_choices[i],"-")[1] != calc_type ? append!(remove_indices,[i]) : nothing
	end
	deleteat!(file_choices,remove_indices)

	remove_indices = []
	for i in 1:length(file_choices)
		current_file = file_choices[i]
		current_file_params_dict = get_params_dict_from_filename(current_file)
		for params in keys(params_dict)
			try
				params_dict[params] in current_file_params_dict[params] ? nothing : append!(remove_indices,i)
			catch
				try
					current_file_metadata_dict = read_data_jld2(current_file,; output_bool=false)[2]
					params_dict[params] in current_file_metadata_dict[params] ? nothing : append!(remove_indices,i)
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
	if potential_type in ["densdens","occs","ttn","Y-dir-GF","X-dir-GF","Y-dir-current","X-dir-current"]
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
	savefig(file_name)
	cd(current_location)
	println("Figure Saved, File Closed: $file_name")
	return
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
	
	binary_file = jldopen(file_name,"w")
	
	if !isnothing(metadata)
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
		println("The directory $location is not accessible")
		location = pwd()
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

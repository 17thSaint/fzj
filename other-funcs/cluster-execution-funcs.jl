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

function make_args_dict(args)
	parameters_dict = Dict()
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
		parameters_dict[key] = value
	end
	return parameters_dict
end

function make_parameters_filename(param_dict)
	param_filename = ""
	for key in keys(param_dict)
		value = string(param_dict[key])
		param_filename *= "$key-$value-"
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
	
	split_filename = split(filename,"-")[2:end]
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








"fin"

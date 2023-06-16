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









"fin"

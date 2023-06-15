

function make_args_dict(args)
	parameters_dict = Dict()
	for i in 1:2:length(args)
		key = args[i]
		value = args[i+1]
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

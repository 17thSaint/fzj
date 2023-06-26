include("../other-funcs/cluster-execution-funcs.jl")

args_dict = make_args_dict(ARGS)

params_string = ""
for key in keys(args_dict)
	global params_string *= "and " * "$key = $(args_dict[key]), "
end
println(params_string)


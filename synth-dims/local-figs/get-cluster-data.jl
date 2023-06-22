include("../long-range-ttn.jl")
using PyPlot



function plot_data(version,params_desired; kwargs...)
	data_files = find_data_file(params_desired,version)
	data_dicts = [Dict() for i in 1:length(data_files)]
	for i in 1:length(data_files)
		filename = data_files[i]
		data_dict = read_data_jld2(filename,"/home/patrick/fzj/main-git/cluster-data")[1]
		data_dicts[i] = data_dict
		num_parts = get_params_dict_from_filename(filename)["particles"]
		plot_title = "LR = $lr, Mag Off = $mag_off, Np = $num_parts"
		if version == "densdens"
			plot_allAVG_densdenscorr(nothing,nothing,nothing; data_dict=data_dict,plot_title=plot_title,name=filename, kwargs...)
		elseif version == "occs"
			plot_occupancy(nothing; data_dict=data_dict,plot_title=plot_title,name=filename, kwargs...)
		end
	end
	return data_dicts,data_files
end

layers = 5
mdim = 100
mag_off = false

lr = 1
params = Dict([("lr",lr),("layers",layers),("mdim",mdim),("mag_off",mag_off)])
rez = plot_data("occs",params; if_save_fig=true,location="/home/patrick/fzj/main-git/local-plots/")




























"fin"

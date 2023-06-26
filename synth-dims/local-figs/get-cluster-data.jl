include("../long-range-ttn.jl")
using PyPlot


function get_ttns(params_desired)
	data_files = find_data_file(params_desired,"ttn")
	data_dicts = [Dict() for i in 1:length(data_files)]
	for i in 1:length(data_files)
		filename = data_files[i]
		data_dict = read_data_jld2(filename,"/home/patrick/fzj/main-git/cluster-data")[1]
		data_dicts[i] = data_dict
	end
	return data_dicts,data_files
end

function plot_data(version,params_desired; kwargs...)
	data_files = find_data_file(params_desired,version)
	data_dicts = [Dict() for i in 1:length(data_files)]
	for i in 1:length(data_files)
		filename = data_files[i]
		data_dict = read_data_jld2(filename,"/home/patrick/fzj/main-git/cluster-data")[1]
		data_dicts[i] = data_dict
		
		#=
		num_parts = get_params_dict_from_filename(filename)["particles"]
		plot_title = "LR = $lr, Mag Off = $mag_off, Np = $num_parts"
		if version == "densdens"
			plot_allAVG_densdenscorr(nothing,nothing,nothing; data_dict=data_dict,plot_title=plot_title,name=filename, kwargs...)
		elseif version == "occs"
			plot_occupancy(nothing; data_dict=data_dict,plot_title=plot_title,name=filename, kwargs...)
		end
		=#
	end
	return data_dicts,data_files
end


layers = 6
mdim = 150
lr = 0
#
if true
params = Dict([("lr",lr),("layers",layers),("mdim",mdim),("mag_off",false)])
all_ttns,all_files = get_ttns(params)
end


#
if false
for i in 1:Int(length(all_ttns)/3)
	println(i)
	ttn_neg = all_ttns[3*(i-1)+1]["ttn"]
	ttn_mid = all_ttns[3*(i-1)+2]["ttn"]
	ttn_pos = all_ttns[3*(i-1)+3]["ttn"]
	filename = all_files_magOFF[i]
	forward_deriv = deriv_bulk_dens(ttn_pos,ttn_mid,0.0001,1)
	back_deriv = deriv_bulk_dens(ttn_mid,ttn_neg,-0.0001,1)
	bds_on[j,i] = bulk_density(ttn,j)
end
end





















"fin"

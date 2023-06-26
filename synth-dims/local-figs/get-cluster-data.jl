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
lr = 1
#
if true
params_magOFF = Dict([("lr",lr),("layers",layers),("mdim",mdim),("mag_off",true)])
all_ttns_magOFF,all_files_magOFF = get_ttns(params_magOFF)
lr0params_magOFF = Dict([("lr",0),("layers",layers),("mdim",mdim),("mag_off",true),("nn_strength",0.0)])
lr0_ttn_magOFF = get_ttns(lr0params_magOFF)
append!(all_ttns_magOFF,[lr0_ttn_magOFF[1][1]])
append!(all_files_magOFF,[lr0_ttn_magOFF[2][1]])
end
if true
params_magON = Dict([("lr",lr),("layers",layers),("mdim",mdim),("mag_off",false)])
all_ttns_magON,all_files_magON = get_ttns(params_magON)
lr0params_magON = Dict([("lr",0),("layers",layers),("mdim",mdim),("mag_off",false),("nn_strength",0.0)])
lr0_ttn_magON = get_ttns(lr0params_magON)
append!(all_ttns_magON,[lr0_ttn_magON[1][1]])
append!(all_files_magON,[lr0_ttn_magON[2][1]])
end

#=
dist1_corrs_off = []
strengths_off = []
if true
for i in 1:length(all_ttns_magOFF)
	println(i)
	ttn = all_ttns_magOFF[i]["ttn"]
	filename = all_files_magOFF[i]
	stren = get_params_dict_from_filename(filename)["nn_strength"]
	title_string = "NN Str = $(stren), Mdim=$mdim, Mag OFF, Np=4"
	append!(strengths_off,[stren])
	rez = get_allAVG_densdenscorr(ttn,[i for i in 1:2*Int(sqrt(2^layers))]; plot_title=title_string,if_save_fig=false,if_plot=false,location=pwd(),name=filename)
	append!(dist1_corrs_off,[rez[1][1]])
	#get_occupancy(ttn; plot_title=title_string,if_save_fig=false,location=pwd(),name=filename)
end
scatter(strengths_off,dist1_corrs_off,label="Mag OFF")
end


dist1_corrs_on = []
strengths_on = []
if true
for i in 1:length(all_ttns_magON)
	println(i)
	ttn = all_ttns_magON[i]["ttn"]
	filename = all_files_magON[i]
	stren = get_params_dict_from_filename(filename)["nn_strength"]
	title_string = "NN Str = $(stren), Mdim=$mdim, Mag ON, Np=4"
	append!(strengths_on,[stren])
	rez = get_allAVG_densdenscorr(ttn,[i for i in 1:2*Int(sqrt(2^layers))]; plot_title=title_string,if_save_fig=false,if_plot=false,location=pwd(),name=filename)
	append!(dist1_corrs_on,[rez[1][1]])
	#get_occupancy(ttn; plot_title=title_string,if_save_fig=false,location=pwd(),name=filename)
end
scatter(strengths_on,dist1_corrs_on,label="Mag On")

legend()
xlabel("NN Strength")
ylabel("NN Correlation")
yscale("log")
end

#
on_poses = []
for i in 1:length(all_ttns_magON)
	ttn = all_ttns_magON[i]["ttn"]
	append!(on_poses,[bulk_density(ttn)])
end
scatter(strengths,on_poses)

off_poses = []
for i in 1:length(all_ttns_magOFF)
	ttn = all_ttns_magOFF[i]["ttn"]
	append!(off_poses,[bulk_density(ttn)])
end
scatter(strengths,off_poses)
=#






















"fin"

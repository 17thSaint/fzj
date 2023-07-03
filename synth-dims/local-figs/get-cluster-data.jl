include("../long-range-ttn.jl")
using PyPlot


function get_ttns(params_desired)
	data_files = find_data_file(params_desired,"ttn")
	if length(data_files) == 0
		println("No Matching Files Found")
		return nothing,nothing
	else
		data_dicts = [Dict() for i in 1:length(data_files)]
		for i in 1:length(data_files)
			filename = data_files[i]
			data_dict = read_data_jld2(filename,"/home/patrick/fzj/main-git/cluster-data")[1]
			data_dicts[i] = data_dict
		end
		return data_dicts,data_files
	end
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

#
layers = 8
mdim = 250
lr = 0
central_alpha = 0.18
change = 0.0001
alpha_vals = [central_alpha - change, central_alpha, round(central_alpha + change,digits=4)]
#nnstren = 0.01
#

if false
params = Dict([("lr",lr),("layers",layers),("mag_off",false),("alpha",alpha_vals)])
all_ttns,all_files = get_ttns(params)
end
#
params_naming = Dict([("lr",lr),("layers",layers),("mag_off",false),("alpha",central_alpha)])
if false
	filename = "derivbulkdens-" * make_parameters_filename(params_naming)
	println(filename)
	metadata_dict = Dict([("all_files",all_files)])
end

#
if false
sizehere = Int(length(all_ttns)/3)
bds = zeros(3,sizehere)
avgs = zeros(3,sizehere)
errs = zeros(3,sizehere)
alphas = []#zeros(1,Int(length(all_ttns)/3))
for j in [1,2,3]
for i in 1:sizehere
	println(i)
	#
	ttn_neg = all_ttns[3*(i-1)+1]["ttn"]
	filename_neg = all_files[3*(i-1)+1]
	alpha_neg = get_params_dict_from_filename(filename_neg)["alpha"]
	#append!(alphas,[alpha_neg])
	#append!(bds,[bulk_density(ttn_neg,3)])
	#
	ttn_mid = all_ttns[3*(i-1)+2]["ttn"]
	filename_mid = all_files[3*(i-1)+2]
	alpha_mid = get_params_dict_from_filename(filename_mid)["alpha"]
	
	#append!(bds,[bulk_density(ttn_mid,3)])
	#rads,vals = radial_box_dist(ttn_mid)
	#plot(rads,vals,"-p",label="$alpha_mid")
	#legend()
	#
	ttn_pos = all_ttns[3*(i-1)+3]["ttn"]
	filename_pos = all_files[3*(i-1)+3]
	alpha_pos = get_params_dict_from_filename(filename_pos)["alpha"]
	#append!(alphas,[alpha_pos])
	#append!(bds,[bulk_density(ttn_pos,3)])
	#
	if j == 1
		get_occupancy(ttn_neg; if_plot=true,plot_title="Alpha=$alpha_neg, LR=$lr, MagON",if_save_fig=true,name=filename_neg)
		get_occupancy(ttn_pos; if_plot=true,plot_title="Alpha=$alpha_pos, LR=$lr, MagON",if_save_fig=true,name=filename_pos)
		get_occupancy(ttn_mid; if_plot=true,plot_title="Alpha=$alpha_mid, LR=$lr, MagON",if_save_fig=true,name=filename_mid)
		append!(alphas,[alpha_mid])
	end
	
	deriv12 = deriv_bulk_dens(ttn_mid,ttn_neg,alpha_mid-alpha_neg,j)
	deriv13 = deriv_bulk_dens(ttn_pos,ttn_neg,alpha_pos-alpha_neg,j)
	deriv23 = deriv_bulk_dens(ttn_pos,ttn_mid,alpha_pos-alpha_mid,j)
	bds[1,i] = deriv12
	bds[2,i] = deriv13
	bds[3,i] = deriv23
end
#scatter(alphas,bds[1,:],label="12, $j")
#scatter(alphas,bds[2,:],label="13, $j")
#scatter(alphas,bds[3,:],label="23, $j")
avg = [mean(bds[:,k]) for k in 1:sizehere]
avgs[j,:] = avg
errors = [std(bds[:,k]) for k in 1:sizehere]
errs[j,:] = errors
println("Value = $avg with error = $errors for Edge Width = $j")
#errorbar(alphas,avgs[j,:],yerr=[errs[j,:],errs[j,:]],marker="o",label="AVG $j")
#plot(alphas,avgs[j,:],"-p",label="AVG $j")
end
#plot(alphas,[0.5 for pp in 1:length(alphas)],label="1/2")
#plot(alphas,[2/3 for pp in 1:length(alphas)],label="2/3")
#xlabel("Bulk Density, Alpha")
#title("Derivative Bulk Density wrt Flux Density = Hall Conductance")
#legend()
end
#

if true
	data_dict = Dict([("avgs",avgs),("errs",errs),("alphas",alphas)])
	write_data_jld2(filename,data_dict,pwd(),metadata_dict)
end
#

















"fin"

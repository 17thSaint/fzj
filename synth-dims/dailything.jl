include("long-range-ttn.jl")
include("fqh_effective.jl")
using PyPlot

#
alpha_start = 0.215
alpha_end = 0.295
alpha_count = 9
new_new_new_alpha_vals = [alpha_start + (i-1)*(alpha_end-alpha_start)/(alpha_count-1) for i in 1:alpha_count]
change = 0.0001
#

#params_dict = Dict([("nbosons",60),("L",96)])
loc = "/home/patrick/fzj/main-git/cluster-data/orsay-sept23"
#all_files = find_data_file(params_dict,"mps","jld2",loc)
#display(all_files)

wavefuncs = Dict()
bulkdenses = [0.0 for i in 1:3*alpha_count]
#alpha_vals = [0.0 for i in 1:length(all_files)]
new_new_new_derivbds = [0.0 for i in 1:alpha_count]
#=
if false
for (idx,file) in enumerate(all_files)
	alpha = get_params_dict_from_filename(file)["alpha"]
	alpha_vals[idx] = alpha
	wavefunc = read_data_jld2(file,loc)[1]["mps"]
	bd = bulk_density(wavefunc,10)
	bulkdenses[idx] = bd
end
plot(alpha_vals,bulkdenses,"-p")
end

if false
deriv_alphas = []
for i in 1:Int(length(alpha_vals)/3)
	append!(derivbds,[(bulkdenses[3*(i-1)+3]-bulkdenses[3*(i-1)+1])/(alpha_vals[3*(i-1)+3]-alpha_vals[3*(i-1)+1])])
	append!(deriv_alphas,[alpha_vals[3*(i-1)+2]])
end
end
=#

seed_data = read_data_jld2("mps-if_periodic_phys-false-nflavors-5-alpha-0.2509-if_periodic_synth-true-nbosons-60-L-96-mk-3.jld2",loc)

#
for i in 1:alpha_count
	println(i/alpha_count)
	alph = new_new_new_alpha_vals[i]
	psi_mid = run_mps_new_variable(seed_data[1]["mps"],seed_data[2],Dict([("phi",2*pi*alph),("location",loc),("mdim",100)]))
	psi_m,psi_p = varied_alpha_wavefuncs(psi_mid,merge(seed_data[2],Dict([("phi",2*pi*alph),("location",loc),("mdim",100)])),change)
	wavefuncs[string(round(alph,digits=4))] = [psi_m,psi_mid,psi_p]
	bd_m = bulk_density(psi_m,10)
	bulkdenses[3*(i-1)+1] = bd_m
	bd_mid = bulk_density(psi_mid,10)
	bulkdenses[3*(i-1)+2] = bd_mid
	bd_p = bulk_density(psi_p,10)
	bulkdenses[3*i] = bd_p
	new_new_new_derivbds[i] = (bd_p - bd_m)/(2*change)
	println("DBD = ",derivbds[i])
end
#

if true
#fig = figure()
scatter(new_new_new_alpha_vals,new_new_new_derivbds)
#title("Deriv Bulk Density vs Flux Density")
end

#fig2 = figure()
#plot(Iterators.flatten([[alpha_vals[i]-change,alpha_vals[i],alpha_vals[i]+change] for i in 1:alpha_count]),bulkdenses,"-p")
#title("Bulk Density vs Flux Density")
































"fin"

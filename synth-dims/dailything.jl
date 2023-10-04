include("long-range-ttn.jl")
include("fqh_effective.jl")
using PyPlot

alpha_start = 0.251
alpha_end = 0.501
alpha_count = 5
alpha_vals = [alpha_start + (i-1)*(alpha_end-alpha_start)/(alpha_count-1) for i in 1:alpha_count]
change = 0.0001

wavefuncs = Dict()
bulkdenses = [0.0 for i in 1:3*alpha_count]
derivbds = [0.0 for i in 1:alpha_count]

filename1 = "mps-if_periodic_phys-false-nflavors-5-alpha-0.125-if_periodic_synth-true-nbosons-60-L-96.jld2"
#filename2 = "mps-if_periodic_phys-false-nflavors-5-alpha-0.25-if_periodic_synth-true-nbosons-60-L-96.jld2"
#seed_data1 = read_data_jld2(filename1,"local-figs/orsay-sept23")
#seed_data2 = read_data_jld2(filename2,"local-figs/orsay-sept23")

cd("local-figs/orsay-sept23")

if false
all_files = readdir()
deleteat!(all_files,findall(i->split(all_files[i],"-")[1] != "mps",1:length(all_files)))
display(all_files)
end

if true
for (idx,f) in enumerate(all_files)
	file = jldopen(f,"r+")
	#try
		file["metadata"]["location"] =  "/home/patrick/fzj/main-git/cluster-data/orsay-sept23/"
		close(file)
	#catch
	#	println("Didnt work on file ",f)
	#	close(file)
	#end
	
end
end



cd("../..")




#=
seed_data = seed_data1
seed_data[1] = seed_data2

for i in 1:alpha_count
	println(i/alpha_count)
	alph = alpha_vals[i]
	psi_mid = run_mps_new_variable(seed_data[1]["mps"],seed_data[2],Dict([("phi",2*pi*alph),("typechange","nothing"),("mdim",200)]))
	psi_m,psi_p = varied_alpha_wavefuncs(psi_mid,merge(seed_data[2],Dict([("phi",2*pi*alph),("mdim",200)])),change)
	wavefuncs[string(round(alph,digits=4))] = [psi_m,psi_mid,psi_p]
	bd_m = bulk_density(psi_m,10)
	bulkdenses[3*(i-1)+1] = bd_m
	bd_mid = bulk_density(psi_mid,10)
	bulkdenses[3*(i-1)+2] = bd_mid
	bd_p = bulk_density(psi_p,10)
	bulkdenses[3*i] = bd_p
	derivbds[i] = (bd_p - bd_m)/(2*change)
	println("DBD = ",derivbds[i])
end

#fig = figure()
plot(alpha_vals,derivbds,"-p")
#title("Deriv Bulk Density vs Flux Density")

#fig2 = figure()
#plot(Iterators.flatten([[alpha_vals[i]-change,alpha_vals[i],alpha_vals[i]+change] for i in 1:alpha_count]),bulkdenses,"-p")
#title("Bulk Density vs Flux Density")

=#
# Work looking at more flux values, 0.25 has an oddly high bond dimension....
































"fin"

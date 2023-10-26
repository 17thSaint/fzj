using PyPlot
include("../review-practice-codes/ttn.jl")

function reconstruct_jld2_network(wavefunc)
	lattices_jld2 = (wavefunc.net).lattices
	allnodes = [lattices_jld2[i].lat for i in 1:length(lattices_jld2)]
	alldims = [lattices_jld2[i].dims for i in 1:length(lattices_jld2)]
	all_simplelattices = [TTNKit.SimpleLattice(allnodes[i],alldims[i]) for i in 1:length(lattices_jld2)]
	return BinaryNetwork(all_simplelattices)
end

function reconstruct_jld2_ttn(wavefunc)
	remade_network = reconstruct_jld2_network(wavefunc)
	return TTNKit.TreeTensorNetwork(wavefunc.data,wavefunc.ortho_direction,wavefunc.ortho_center,remade_network)
end

#=
dataloc = "../cluster-data/chemical-potential"
layers = 6
tstren = 0.5
params = Dict([("layers",layers),("t",tstren)])
filenames = find_data_file(params,"ttn","jld2",dataloc)

metadatas = []
allpsis = []
for f in filenames
	alldata = read_data_jld2(f,dataloc)
	append!(metadatas,[alldata[2]])
	fileversion_psi = alldata[1]["ttn"]
	rebuilt_psi = reconstruct_jld2_ttn(fileversion_psi)
	append!(allpsis,[rebuilt_psi])
end
=#

num_sites = 2^layers
alldens = []
chem_vals = []
#ydir_greens = []
for (i,psi) in enumerate(allpsis)
	println(i/length(allpsis))
	chemstren = metadatas[i]["chem_strength"]
	dens = sum(get_occupancy(psi; if_plot=true, plot_title="$(round(chemstren,digits=4))"))/num_sites
	#ydir_green = get_ydir_greenfunc(psi; if_plot=true, plot_title="$(round(chemstren,digits=4))")
	append!(alldens,[dens])
	append!(chem_vals,[chemstren])
end
fig = figure()
plot(chem_vals,alldens,"-p")
xlabel("Chemical Potential Strength")
ylabel("Density")
#


































"fin"

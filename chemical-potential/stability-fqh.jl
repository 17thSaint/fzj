println("Starting Now")

function include_other_files(all_files)
	get_to_fzj = split(pwd(),"fzj")[1]
	if typeof(all_files) == String
		all_files = [all_files]
	end
	for file in all_files
		occursin("main-git",pwd()) ? include(get_to_fzj * "fzj/main-git/" * file) : include(get_to_fzj * "fzj/" * file)
		println("Included $file")
	end
end
include_other_files("review-practice-codes/ttn.jl")
println("Added All files and packages")


#parameter_iteration = parse(Int,ENV["ENCAP_PROCID"])
#params_dict = make_args_dict(ARGS,parameter_iteration)
params_dict = Dict([("layers",4),("mdim",18),("if_save_data",false)])
display(params_dict)

layers = get(params_dict, "layers", 4)#Done
num_sites = 2^layers
filling = 1/2
mdim = get(params_dict, "mdim", 100)#Done
nsweeps = get(params_dict, "nsweeps", 50)#Done
if_save_data = get(params_dict, "if_save_data", false)#Done
noise = get(params_dict, "noise", [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0.0])#Done
expander = TTNKit.DefaultExpander(get(params_dict, "expander_coeff", 0.2))#Done
u_strength = get(params_dict, "u_strength", 1.0)#Done
syms = get(params_dict, "syms", false)#Done
max_occ = get(params_dict, "maxocc", 2)
if_gpu = get(params_dict, "if_gpu", false)#Done
seed_ttn = get(params_dict, "seed_ttn", nothing)#Done
net = TTNKit.BinaryRectangularNetwork(layers, TTNKit.ITensorNode, "Boson";conserve_qns=syms,dim=max_occ+1)#Done
dataloc = get(params_dict, "dataloc", get_folder_location("cluster-data/chemical-potential"))#Done
open_cores = get(params_dict, "open_cores", 1)#Done
particle_type = "Boson"
if typeof(open_cores) != String
	BLAS.set_num_threads(Int(open_cores))
end

if_periodic = false #Done
if_chem = true#Done
no_magF = false

#t_strength = get(params_dict, :t_strength, 0.02)
#chem_strength = get(params_dict, :chem_strength, 0.0)
alpha = get(params_dict, "alpha", 1/4)#Done
t_strength = round(get(params_dict, "t", 0.5),digits=4)#Done
chem_strength = round(get(params_dict, "chem", 1.0),digits=4)#Done
num_particles = Int(round(filling*alpha*num_sites,digits=0))#Done

println("Made All Variables")

naming_dict = merge(Dict([("layers",layers),("alpha",alpha),("maxocc",max_occ),("chem",chem_strength),("t",t_strength)]),params_dict)
if "dataloc" in collect(keys(naming_dict))
	delete!(naming_dict,"dataloc")
end
if "if_save_data" in collect(keys(naming_dict))
	delete!(naming_dict,"if_save_data")
end
model_paras = (nrgtol=1E-6,if_periodic_synth=false,if_periodic_phys = false, if_chem = if_chem, chem_strength = chem_strength, u_strength = u_strength, mdim = mdim, num_sweeps = nsweeps, noise = noise, if_save_data = if_save_data, sweep_type = "dmrg", syms = syms, phi = alpha, ttn_net = net, seed_ttn = seed_ttn, if_gpu = if_gpu, layers = layers, t_strength = t_strength, filling = filling, location = dataloc, particles = num_particles, open_cores = open_cores, part_type = particle_type, expander = expander)
metadata_dict = named_tuple_to_dict(model_paras)
filename = "ttn-" * make_parameters_filename(naming_dict)
display(dataloc * "/" * filename)

ham_op = get_hofstadter_interacting_hamilt(net,t_strength,alpha; model_paras...)
display(ham_op)
og_ttn, og_ham, dmsp = find_ground_state(layers,num_particles; model_paras...,ham_op = ham_op,metadata=metadata_dict,name=filename)

get_occupancy(dmsp.ttn; if_plot=false)

#=
ts_count = 1
mu_count = 1
mus = [(i-1)*2/(mu_count-1) for i in 1:mu_count]
ts = [0.001 + (i-1)*(0.5-0.001)/(ts_count-1) for i in 1:ts_count]
denses = zeros(mu_count,ts_count)

for (idt,t_strength) in enumerate(ts)
	global naming_dict["t"] = round(t_strength,digits=4)
	for (idc,chem_strength) in enumerate(mus)
		model_paras = (if_periodic = if_periodic, if_chem = if_chem, chem_strength = chem_strength, u_strength = u_strength, max_dim = mdim, num_sweeps = nsweeps, noise = noise, if_save_data = if_save_data, sweep_type = "dmrg", syms = syms, phi = alpha, ttn_net = net, seed_ttn = seed_ttn, if_gpu = if_gpu, layers = layers, t_strength = t_strength, filling = filling, location = dataloc)
		metadata_dict = named_tuple_to_dict(model_paras)
		global naming_dict["chem"] = round(chem_strength,digits=4)
		filename = "ttn-" * make_parameters_filename(naming_dict)
		display(filename)

		ham_op = get_hofstadter_interacting_hamilt(net,t_strength,alpha; model_paras...)

		og_ttn, og_ham, dmsp = build_full_harperhofstadter(layers,num_particles,t_strength,filling; model_paras...,ham_op = ham_op,metadata=metadata_dict,name=filename)
		psi = dmsp.ttn
		occs = get_occupancy(psi; if_plot=false)
		denses[idc,idt] = sum(occs)/num_sites
		println("Density Value = ",denses[idc,idt])
		
	end
end

#imshow(denses)
#colorbar()
#
denses = [[]]
all_files = find_data_file(Dict([("layers",4),("alpha",0.25)]),"ttn","jld2",dataloc)
for f in all_files
	psi = read_data_jld2(f,dataloc)[1]["ttn"]
	dens = sum(get_occupancy(psi; if_plot=false))/num_sites
	params = get_params_dict_from_filename(f)
	scatter3D([params["t"]],[params["chem"]],[dens],c="b")
end
xlabel("Hopping Strength, t")
ylabel("Chemical Potential, mu")
zlabel("Density")
=#






































"fin"

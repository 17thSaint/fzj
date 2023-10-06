include("../review-practice-codes/ttn.jl")


params_dict = Dict([("layers",4)])

layers = get(params_dict, :layers, 4)
num_sites = 2^layers
filling = 1/2
mdim = get(params_dict, :mdim, 100)
nsweeps = 5
if_save_data = true
noise = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0.0]
u_strength = 1.0
syms = false
max_occ = get(params_dict, :max_occ, 1)
if_gpu = true
seed_ttn = get(params_dict, :seed_ttn, nothing)
net = TTNKit.BinaryRectangularNetwork(layers, TTNKit.ITensorNode, "Boson";conserve_qns=syms,dim=max_occ+1)
dataloc = "../cluster-data/chemical-potential"

if_periodic = false
if_chem = true
no_magF = false

#t_strength = get(params_dict, :t_strength, 0.02)
#chem_strength = get(params_dict, :chem_strength, 0.0)
alpha = get(params_dict, :alpha, 1/4)
num_particles = Int(round(filling*alpha*num_sites,digits=0))

naming_dict = Dict([("layers",layers),("alpha",alpha),("maxocc",max_occ)])
#
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
#=
all_files = find_data_file(Dict([("alpha",0.25)]),"ttn","jld2",dataloc)
for f in all_files
	psi = read_data_jld2(f,dataloc)[1]["ttn"]
	dens = sum(get_occupancy(psi; if_plot=false))/num_sites
	params = get_params_dict_from_filename(f)
	scatter3D([params["t"]],[params["chem"]],[dens*0.25*sqrt(num_sites)],c="b")
end
xlabel("Hopping Strength, t")
ylabel("Chemical Potential, mu")
zlabel("LLL Filling")
=#






































"fin"

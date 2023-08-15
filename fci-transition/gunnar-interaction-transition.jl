include("../synth-dims/long-range-ttn.jl")
using PyPlot,CurveFit

function nn_HH_ham(net,t_strength,phi; kwargs...)
	resulting_ham = []
	phys_edge_length,virt_edge_length = get_lattice_dims(net)
	
	onsite_strength = get(kwargs, :onsite_strength, 1.0)	
	if_onsite_int = get(kwargs, :if_onsite_int, true)
	if_periodic = get(kwargs, :if_periodic, true)
	if_hopping = get(kwargs, :if_hopping, true)
	if_chem = get(kwargs, :if_chem, false)
	chem_strength = get(kwargs, :chem_strength, 0.0)
	if_nn_int = get(kwargs, :if_nn_int, false)
	nn_int_strength = get(kwargs, :nn_int_strength, 0.0)
	no_magF = get(kwargs, :no_magF, false)
	
	println(join((onsite_strength,if_onsite_int,if_chem,chem_strength,if_nn_int,nn_int_strength),", "))
	
	lat = TTNKit.physical_lattice(net)
	
	if if_hopping
		hopping = TTNKit.OpSum()
		for (s1,s2) in TTNKit.nearest_neighbours(lat,collect(1:TTNKit.number_of_sites(lat)); periodic=if_periodic)
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)
						
			coeff = get_inter_coeff(s1_coord,s2_coord,t_strength,phi,phys_edge_length,virt_edge_length; kwargs...)
			
			hopping += (coeff,"Adag",s1_coord,"A",s2_coord)
			hopping += (conj(coeff),"Adag",s2_coord,"A",s1_coord)
		end
		append!(resulting_ham,[hopping])
	end
	
	if if_onsite_int
		interaction = TTNKit.OpSum()
		for j in TTNKit.eachindex(lat)
			interaction += (onsite_strength,"Adag * Adag * A * A",TTNKit.coordinate(lat,j))
		end
		append!(resulting_ham,[interaction])
	end
	
	if if_nn_int
		nn_int = TTNKit.OpSum()
		nns = TTNKit.nearest_neighbours(lat,collect(TTNKit.eachindex(lat));periodic=if_periodic)
		for (n1,n2) in nns
			nn_int += (nn_int_strength,"Adag * A",TTNKit.coordinate(lat,n1),"Adag * A",TTNKit.coordinate(lat,n2))
		end
		append!(resulting_ham,[nn_int])
	end
	
	if if_chem
		chem = TTNKit.OpSum()
		for i in TTNKit.eachindex(lat)
			chem -= (chem_strength,"N",TTNKit.coordinate(lat,i))
		end
		append!(resulting_ham,[chem])
	end
	
	if length(resulting_ham) > 1
		return sum(resulting_ham)
	else
		return resulting_ham[1]
	end
end

function correlation_length(ttn,direction="Y")
	greensfunc = direction == "Y" ? get_ydir_greenfunc(ttn;if_plot=false) : get_xdir_greenfunc(ttn;if_plot=false)
	dists = greensfunc[1] .- 1.0
	edge_length = size(dists)[1]
	xi_left = [0.0 for i in 1:edge_length]
	xi_right = [0.0 for i in 1:edge_length]
	for i in 1:edge_length
		center = findfirst(x->greensfunc[2][i,x]==minimum(greensfunc[2][i,:]),[i for i in 1:edge_length])
		front_right,inv_xi_right = exp_fit(dists[i,center:end],greensfunc[2][i,center:end])
		xi_right[i] = 1/inv_xi_right
	end
	return xi_left,xi_right
end

#
layer_count = 6

if layer_count % 2 == 0
	edge_sites = Int(sqrt(2^layer_count))
	num_particles = Int(edge_sites/2)
else
	edge_sites = Int(sqrt(2^(layer_count+1)))
	num_particles = Int(sqrt(2^(layer_count+1))/2)
end
nu = 1/2
tot_sites = 2^layer_count
ts = 0.001
mag_off = false
alpha = num_particles/(nu * (tot_sites))
mdim = 300

sweep_type = "dmrg"
max_occ = 1
if_per = true
evolve = true
chemical = false
mu = 0.0
if_NN = true
if_onsite = true
onsite_stren = 10.0
nswps = 6
expan = TTNKit.DefaultExpander(0.5)

if_change = false
change = if_change ? 0.001 : 0.0

if_gpu = false
plotting = false
save_plot = false
save_data = true
loc = "../cluster-data/"

nns_start = 1.51
nns_end = 1.59
nns_count = 5
nn_strens = [nns_start + (i-1)*(nns_end-nns_start)/(nns_count-1) for i in 1:nns_count]
wavefuncs = []
#append!(nn_strens,[1.4,1.5,1.6,1.7,1.8])

metadata = metadata_dict = Dict([("if_per",if_per),("mag_off",mag_off),("chemical",chemical),("mu",mu),("ts",ts),("nu",nu),("layers",layer_count),("particles",num_particles),("alpha",alpha),("mdim",mdim),("nswps",nswps),("max_occ",max_occ),("sweep_type",sweep_type),("if_change",if_change),("change",change),("if_nn_int",if_NN)])

net = build_HH_net(layer_count; syms=true)
states = fill("0", tot_sites)
old_ttn = TTNKit.ProductTreeTensorNetwork(net,states)		
seed_ttn = initialize_ttn(old_ttn,mdim,num_particles)

model_paras = (ttn_net=net,seed_ttn=seed_ttn,if_nn_int=if_NN,if_onsite_int=if_onsite,onsite_strength=onsite_stren,if_save_data=save_data,location=loc,metadata=metadata_dict,max_dim=mdim, num_sweeps=nswps,phi=alpha, if_periodic=if_per,max_occ=max_occ,if_sweep=evolve,sweep_type=sweep_type,expander=expan,if_chem=chemical,chem_strength=mu,no_magF=mag_off,if_gpu=if_gpu,output_level=0)
fig = figure()
xlabel("Nearest-Neighbor Interaction Strength")
ylabel("Correlation Length")

for nnst in nn_strens
	datafile_name = "layers-$layer_count-particles-$num_particles-mdim-$mdim-mag-$(!mag_off)-nn_strength-$nnst"

	#prev_ttn = read_data_jld2("ttn-"*datafile_name*".jld2",loc)[1]["ttn"]
	#prev_net = TTNKit.network(prev_ttn)
	
	starting = time()
	ham = nn_HH_ham(net,ts,alpha; model_paras...,nn_int_strength=nnst)
	og_ttn, hamilt, dm_sp = build_full_harperhofstadter(layer_count,num_particles,ts,nu; ham_op=ham,model_paras...,nn_int_strength=nnst,name="ttn-"*datafile_name)
	total_time = time() - starting
	println("Running time = $total_time")
	append!(wavefuncs,[dm_sp.ttn])
	cl = correlation_length(dm_sp.ttn)[2]
	scatter([nnst for j in 1:length(cl)],cl)

	#get_ydir_greenfunc(dm_sp.ttn;plot_title="$nnst")
end
#







































"fin"

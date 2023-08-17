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
	if_fermion = get(kwargs, :if_fermion, false)
	
	creation = if_fermion ? "Cdag" : "Adag"
	annihilation = if_fermion ? "C" : "A"
	println(join((onsite_strength,if_onsite_int,if_chem,chem_strength,if_nn_int,nn_int_strength),", "))
	
	lat = TTNKit.physical_lattice(net)
	
	if if_hopping
		hopping = TTNKit.OpSum()
		for (s1,s2) in TTNKit.nearest_neighbours(lat,collect(1:TTNKit.number_of_sites(lat)); periodic=if_periodic)
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)
						
			coeff = get_inter_coeff(s1_coord,s2_coord,t_strength,phi,phys_edge_length,virt_edge_length; kwargs...)
			
			hopping += (coeff,creation,s1_coord,annihilation,s2_coord)
			hopping += (conj(coeff),creation,s2_coord,annihilation,s1_coord)
		end
		append!(resulting_ham,[hopping])
	end
	
	if !if_fermion && if_onsite_int
		interaction = TTNKit.OpSum()
		for j in TTNKit.eachindex(lat)
			interaction += (onsite_strength,join([creation,creation,annihilation,annihilation]," * "),TTNKit.coordinate(lat,j))
		end
		append!(resulting_ham,[interaction])
	end
	
	if if_nn_int
		nn_int = TTNKit.OpSum()
		nns = TTNKit.nearest_neighbours(lat,collect(TTNKit.eachindex(lat));periodic=if_periodic)
		for (n1,n2) in nns
			nn_int += (nn_int_strength,join([creation,annihilation]," * "),TTNKit.coordinate(lat,n1),join([creation,annihilation]," * "),TTNKit.coordinate(lat,n2))
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

function spectra(matrices)
    size_matrices = size(matrices[1])[1]
    num_matrices = length(matrices)
    eigenvalue_vectors = [zeros(size_matrices,size_matrices).*im for i in 1:num_matrices]
    eigenvalues = [[] for i in 1:num_matrices]

    for i in 1:length(matrices)
        matrix = matrices[i]
        eigen_vals, eigen_vecs = eigen(matrix)
        eigenvalue_vectors[i] = eigen_vecs
        eigenvalues[i] = eigen_vals
    end

    num_eigenvalues = size(eigenvalues[1], 1)
    spectra = Matrix(undef, num_eigenvalues, num_matrices)

    for i in 1:num_matrices
        spectra[:, i] = eigenvalues[i]
    end

    return spectra
end

#
if true
layer_count = 4
if_fermion = true

if layer_count % 2 == 0
	edge_sites = Int(sqrt(2^layer_count))
	num_particles = Int(edge_sites/2)
else
	edge_sites = Int(sqrt(2^(layer_count+1)))
	num_particles = Int(sqrt(2^(layer_count+1))/2)
end
nu = 1/3
tot_sites = 2^layer_count
ts = 1.0
mag_off = false
alpha = num_particles/(nu * (tot_sites))
mdim = 100

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
noise = 0.0#[1E-2,0]

if_change = false
change = if_change ? 0.001 : 0.0

if_gpu = false
plotting = false
save_plot = false
save_data = false
loc = "../cluster-data/"
end

if false
all_files = find_data_file(Dict([("ts",1.0)]),"ttn")
wavefuncs = [read_data_jld2(f,loc)[1]["ttn"] for f in all_files]
end

#
nns_start = 1.0
nns_end = 2.0
nns_count = 5
nn_strens = [nns_start + (i-1)*(nns_end-nns_start)/(nns_count-1) for i in 1:nns_count]
wavefuncs = []
rhos = []
#append!(nn_strens,[1.4,1.5,1.6,1.7,1.8])

metadata_dict = Dict([("if_per",if_per),("mag_off",mag_off),("chemical",chemical),("mu",mu),("ts",ts),("nu",nu),("layers",layer_count),("particles",num_particles),("alpha",alpha),("mdim",mdim),("nswps",nswps),("max_occ",max_occ),("sweep_type",sweep_type),("if_change",if_change),("change",change),("if_nn_int",if_NN),("if_fermion",if_fermion)])
#
net = build_HH_net(layer_count; syms=true,if_fermion=if_fermion)
#=states = fill("0", tot_sites)
old_ttn = TTNKit.ProductTreeTensorNetwork(net,states)		
seed_ttn = initialize_ttn(old_ttn,mdim,num_particles; if_fermion=if_fermion)
=#
model_paras = (noise=noise,if_fermion=if_fermion,ttn_net=net,if_nn_int=if_NN,if_onsite_int=if_onsite,onsite_strength=onsite_stren,if_save_data=save_data,location=loc,metadata=metadata_dict,max_dim=mdim, num_sweeps=nswps,phi=alpha, if_periodic=if_per,max_occ=max_occ,if_sweep=evolve,sweep_type=sweep_type,expander=expan,if_chem=chemical,chem_strength=mu,no_magF=mag_off,if_gpu=if_gpu,output_level=0)
#fig = figure()
#xlabel("Nearest-Neighbor Interaction Strength")
#ylabel("Correlation Length")

for nnst in nn_strens
	datafile_name = "layers-$layer_count-particles-$num_particles-mdim-$mdim-mag-$(!mag_off)-nn_strength-$nnst-ts-$ts-if_fermion-$if_fermion"

	#prev_ttn = read_data_jld2("ttn-"*datafile_name*".jld2",loc)[1]["ttn"]
	#prev_net = TTNKit.network(prev_ttn)
	
	starting = time()
	ham = nn_HH_ham(net,ts,alpha; model_paras...,nn_int_strength=nnst)
	og_ttn, hamilt, dm_sp = build_full_harperhofstadter(layer_count,num_particles,ts,nu; ham_op=ham,model_paras...,nn_int_strength=nnst,name="ttn-"*datafile_name)
	total_time = time() - starting
	println("Running time = $total_time")
	append!(wavefuncs,[dm_sp.ttn])
	rho = density_matrix(dm_sp.ttn; model_paras...,name="densmat-"*datafile_name)
	append!(rhos,[rho])
	fig = figure()
	imshow(abs.(rho))
	title("$nnst")

	#get_ydir_greenfunc(dm_sp.ttn;plot_title="$nnst")
end
#
#=
wavefuncs = []
rhos = []
nnindexes = [2,7,12,18,24]
for i in nnindexes
	nnst = nn_strens[i]
	datafile_name = "layers-$layer_count-particles-$num_particles-mdim-$mdim-mag-$(!mag_off)-nn_strength-$nnst-ts-$ts"
	rho = read_data_jld2("densmat-"*datafile_name*".jld2",loc)[1]["densmat"]
	append!(rhos,[rho])
	if i == 2
		imshow(abs.(rho))
	end
end
=#
#nns = [nn_strens[i] for i in nnindexes]
fig = figure()
ss = spectra(rhos)
for i in 1:18
	plot(nn_strens,-1 .* log.(10,ss[i,:]),"-p")
end
#





































"fin"

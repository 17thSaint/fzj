include("fqh_effective.jl")
include("long-range-ttn.jl")
using PyPlot

function fix_filling(L,nflavors,nu)
	prod = L * nu * nflavors
	for nbosons in L/2-1:L-2
		inv_alpha = round(prod/nbosons,digits=5)
		if isinteger(inv_alpha)
			println("Found Alpha = 1/$inv_alpha")
			return Int(nbosons),1/inv_alpha
		end
	end
	println("Not Found")
	return nothing,nothing
end
if true
save_nothing = false
params_dict = Dict()
L = 54#get(params_dict, "L", 4)
#nbosons = 12#get(params_dict, "nbosons", nflavors)
#nflavors = 2#get(params_dict, "nflavors", Int(L/2))
t1 = get(params_dict, "t1", 1.0)
t2 = get(params_dict, "t2", 1.0)
U = get(params_dict, "U", 100)
U1 = 4*t1^2/U
U2 = U1/2
conserve_qns = true
if_nn_int = false#get(params_dict, "if_nn_int", false)
if_2ord_pert = false#get(params_dict, "if_2ord_pert", false)
nsweeps = 20
mdim = get(params_dict, "mdim", 100)
noises = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0.0]
if_save_data = save_nothing ? false : true
data_loc = "/home/patrick/fzj/main-git/synth-dims/local-figs/orsay-sept23"
if_periodic = false
#nflavors = 9
alpha = 1/3
phi = 2*pi*alpha

dmrg_obs = TTNKit.DMRGObserver(;energy_tol=10^-3,minsweeps=3)

other_params_dict = Dict([("U",U),("conserve_qns",conserve_qns),("nsweeps",nsweeps),("mdim",mdim),("noise",noises)])
savefig_data = save_nothing ? false : true
savefig = save_nothing ? false : true
if_lines = false

#nbosons,alpha = fix_filling(L,nflavors,1/2)
#alpha = 0.0

wavefuncs = []
rhos = []
nbosons = Int(L/2)
fillings = ["1/2","2/3","1/3"]
for nflavors in [i for i in 2:8]
filename_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",nbosons),("alpha",round(alpha,digits=4)),("if_nn_int",if_nn_int),("if_2ord_pert",if_2ord_pert),("if_periodic",if_periodic)])
#filename_dict_highdens = Dict([("L",L),("nflavors",nflavors),("nbosons",nbosons_highdens),("alpha",round(alpha,digits=4)),("if_nn_int",if_nn_int),("if_2ord_pert",if_2ord_pert),("if_periodic",if_periodic)])

datafile_name = make_parameters_filename(filename_dict)
#datafile_name_highdens = make_parameters_filename(filename_dict_highdens)

model_paras = (t1 = t1, t2 = t2, phi = phi, U1 = U1, U2 = U2, L = L, nflavors = nflavors, nbosons = nbosons, if_nn_int = if_nn_int, if_2ord_pert = if_2ord_pert, mdim = mdim, nsweeps = nsweeps, noise = noises, if_save_data = if_save_data, location = data_loc, if_periodic=if_periodic, name=datafile_name, observer = dmrg_obs)

#model_paras_highdens = (t1 = t1, t2 = t2, phi = phi, U1 = U1, U2 = U2, L = L, nflavors = nflavors, nbosons = nbosons_highdens, if_nn_int = if_nn_int, if_2ord_pert = if_2ord_pert, mdim = mdim, nsweeps = nsweeps, noise = noise, if_save_data = if_save_data, location = data_loc, if_periodic=if_periodic, name=datafile_name_highdens)

#metadata_dict_highdens = merge(named_tuple_to_dict(model_paras_highdens),other_params_dict)
metadata_dict = merge(named_tuple_to_dict(model_paras),other_params_dict)

psi = execute_mps(U1,U2,phi,L,nflavors,nbosons; model_paras...,metadata=metadata_dict)
append!(wavefuncs,[psi])
densmat = correlation_matrix(psi,"FullDag","FullHat") #./ 2.0
append!(rhos,[densmat])
if false
fig = figure()
imshow(real.(densmat))
colorbar()
title("Virt Dim = $(round(alpha,digits=3))")
end
#append!(all_wavefuncs,[psi])
#psi_highdens = execute_mps(U1,U2,phi,L,nflavors,nbosons_highdens; model_paras_highdens...,metadata=metadata_dict_highdens)
#mrez = momentum_occupation(psi,nbosons,100,10.0; model_paras...,plot_title="Virt Dim = $nflavors",if_log=false,if_plot=true)
#append!(mom_occs,[mrez[2]])
end
#
end

for i in 1:7
	plot(expect(wavefuncs[i],"N"),label="$(i+1)")
end
legend()
#title("Nflavors = $nflavors")

#=end

fig = figure()
momenta = [0.0 + (i-1)*(10.0 - 0.0)/(200-1) for i in 1:200]
for i in 1:length(mom_occs)
	plot(momenta./pi,mom_occs[i]./nbosons,label="$(Ls[i])")
end
title("Momentum Distribution range PhysDim NF = $nflavors, w/ 2nd Order")
xlabel("p/pi")
ylabel("Occupation / nbosons")
legend()

append!(allall_wavefuncs,[all_wavefuncs])
append!(all_mom_occs,[mom_occs])

end
=#
#=
for i in 1:length(fillings)
	filling = fillings[i]
	nbosons,alpha = fix_filling(L,nflavors,filling)
	isnothing(nbosons) ? continue : nothing
	#alpha = alphas[i]

	phi  = 2π*alpha
	#phi_right  = 2π*(alpha+change)
	#phi_left  = 2π*(alpha-change)
	#filling = nbosons / (alpha * L * nflavors)
	println("Nu = $filling with Alpha = $alpha and NBosons = $nbosons")
	filename_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",nbosons),("alpha",round(alpha,digits=4)),("if_nn_int",if_nn_int),("if_2ord_pert",if_2ord_pert)])
	datafile_name = make_parameters_filename(filename_dict)
	model_paras = (t1 = t1, t2 = t2, phi = phi, U1 = U1, U2 = U2, L = L, nflavors = nflavors, nbosons = nbosons, if_nn_int = if_nn_int, if_2ord_pert = if_2ord_pert, mdim = mdim, nsweeps = nsweeps, noise = noise, if_save_data = if_save_data, location = data_loc, name = datafile_name)
	metadata_dict = merge(named_tuple_to_dict(model_paras),other_params_dict)
	println(datafile_name)
	if true
	psi = execute_mps(U1,U2,phi,L,nflavors,nbosons; model_paras...,metadata=metadata_dict)
	append!(all_wavefuncs,[psi])
	#println("Done Center")
	
	
	#model_paras_left = copy(model_paras)
	#model_paras_left[:phi] = phi_left
	#psi_left = execute_mps(U1,U2,phi_left,L,nflavors; model_paras...,psi_guess = psi, if_save_data=false)
	#psi_right = execute_mps(U1,U2,phi_right,L,nflavors; model_paras...,psi_guess = psi, if_save_data=false)
	plottitle = "Filling = $filling"
	rez51 = get_occupancy(psi;if_plot=true, if_save_data=savefig_data,name=datafile_name,location=data_loc,metadata=metadata_dict,if_save_fig=savefig,plot_title=plottitle)
	rez11 = get_greenfunc(psi,"virt"; if_plot=true,if_lines=if_lines, if_save_data=savefig_data,name=datafile_name,location=data_loc,metadata=metadata_dict,if_save_fig=savefig,plot_title=plottitle)
	rez21 = get_greenfunc(psi,"phys"; if_plot=true,if_lines=if_lines, if_save_data=savefig_data,name=datafile_name,location=data_loc,metadata=metadata_dict,if_save_fig=savefig,plot_title=plottitle)
	end
	
	#=
	bd_l = bulk_density(psi_left)
	append!(bds,[bd_l])
	bd_c = bulk_density(psi)
	append!(bds,[bd_c])
	bd_r = bulk_density(psi_right)
	append!(bds,[bd_r])
	
	dbd = (bd_r - bd_l) / (2*change)
	println(alpha,", ",dbd)
	dbds[i] = dbd
	=#
	
end
#=
rez5 = get_occupancy(all_wavefuncs[1];)
rez = get_greenfunc(all_wavefuncs[1],"virt"; if_plot=true,if_lines=if_lines)#, if_save_data=savefig_data,name=datafile_name,location=data_loc,metadata=metadata_dict,if_save_fig=savefig)
rez2 = get_greenfunc(all_wavefuncs[1],"phys"; if_plot=true,if_lines=if_lines)#, if_save_data=savefig_data,name=datafile_name,location=data_loc,metadata=metadata_dict,if_save_fig=savefig)
#re3 = get_greenfunc(all_wavefuncs[1],"virt"; if_plot=true,rev=true,if_lines=if_lines)
#rez4 = get_greenfunc(all_wavefuncs[1],"phys"; if_plot=true,rev=true,if_lines=if_lines)
=#
#println(derivs,", ",errs)
#for i in 1:2
#errorbar(alphas,derivs[i,:],yerr=[errs[i,:],errs[i,:]],label="$i")
#end
#legend()



=#
























"fin"

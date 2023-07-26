include("fqh_effective.jl")
include("long-range-ttn.jl")

if false
nflavors = 2
file_params = Dict([("nflavors",nflavors),("alpha",0.0),("nbosons",10),("if_nn_int",false)])
all_files = find_data_file(file_params,"mps","jld2","/home/patrick/fzj/main-git/cluster-data")
end
if true
Ls = []
wavefuncs = []
moms = []
mom_occs = []
for i in 1:length(all_files)
	dats = read_data_jld2(all_files[i],"/home/patrick/fzj/main-git/cluster-data")[1]["mps"]
	append!(wavefuncs,[dats])
	append!(Ls,get_params_dict_from_filename(all_files[i])["L"])
	mrez = momentum_occupation(dats,200,10.0;if_log=false,if_plot=false)
	append!(moms,[mrez[1]])
	append!(mom_occs,[mrez[2]])
end
#end

for i in 1:length(wavefuncs)
	plot(moms[i]./pi,mom_occs[i]./10,label="$(Ls[i])")
end
title("Momentum Distribution range PhysDim NF = 2, 1st Order")
xlabel("p/pi")
ylabel("Occupation / nbosons")
legend()
end
#=
for i in 1:length(all_files)
	fig = figure()
	mm,oc = momentum_occupation(wavefuncs[i],50,15.0,-15.0)
	params = get_params_dict_from_filename(all_files[i])
	filling = params["nbosons"]/(params["alpha"]*L*nflavors)
	println(filling)
	plot(mm,real.(oc),"-p",label="real")
	#plot(mm,imag.(oc),"-p",label="imag")
	yscale("log")
	#xscale("log")
	title("Mom Dist Filling = $filling")
	#
	legend()
end
=#


















#=
s,j = 2,6
#=
model_paras = (t1 = 1.0, t2 = 1.0, phi = phi, U1 = U1, U2 = U2, L = L, nflavors = nflavors, nbosons = nbosons, if_nn_int = false, if_2ord_pert = false, mdim = 200, nsweeps = 20, noise = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0], if_save_data = false)
sidx = siteinds("ExtendedHardcore", L; conserve_qns = true, nflavors = nflavors)
H = MPO(hamiltonian(1.0,1.0,phi,U1,U2,L,nflavors; if_nn_int = false, if_2ord_pert = false), sidx)
#
states = ["0" for i in 1:L]
states[j] = "$(s)"
psi_start = MPS(sidx,states)
#
hop_phys = OpSum()
hop_phys += (1.0, "Cr$(s)", j+1, "Anh$(s)", j)
mpo_hop_phys = MPO(hop_phys,sidx)
psi_hop_phys = apply(mpo_hop_phys,psi_start)
#
hop_virt = OpSum()
hop_virt += (1.0, "Cr$(s+1) * Anh$(s)", j)
mpo_hop_virt = MPO(hop_virt,sidx)
psi_hop_virt = apply(mpo_hop_virt,psi_start)
=#

=#
















"fin"

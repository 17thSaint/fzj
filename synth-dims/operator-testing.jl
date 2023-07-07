include("fqh_effective.jl")
include("long-range-ttn.jl")

L = 6
nflavors = 20
file_params = Dict([("L",L),("nflavors",nflavors),("if_nn_int",true)])
all_files = find_data_file(file_params,"virt-dir-GF","jld2","/home/patrick/fzj/main-git/synth-dims/local-figs")



























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

include("fqh_effective.jl")
include("long-range-ttn.jl")

if false
L = 6
nflavors = 20
file_params = Dict([("L",L),("nflavors",nflavors),("if_nn_int",false)])
all_files = find_data_file(file_params,"mps","jld2","/home/patrick/fzj/main-git/cluster-data")
alphas = []
wavefuncs = []
for i in 1:length(all_files)
	loc_mps = read_data_jld2(all_files[i],"/home/patrick/fzj/main-git/cluster-data")[1]["mps"]
	append!(wavefuncs,[loc_mps])
	append!(alphas,get_params_dict_from_filename(all_files[i])["alpha"])
end
end

function momentum_occupation(wavefunc::MPS,p_count::Int64,p_end=4.0,p_start=0.0)
	num_sites = length(siteinds(wavefunc))
	momenta = [p_start + (i-1)*(p_end - p_start)/(p_count-1) for i in 1:p_count]
	mom_occ = zeros(p_count)*im
	for i in 1:p_count
		momentum = momenta[i]
		exp_vect =  [exp(momentum*j*im) for j in 1:num_sites]
		pos_occ = expect(wavefunc,"N")
		mom_occ[i] = sum(exp_vect .* pos_occ) / sqrt(num_sites)
	end
	return momenta,mom_occ
end

function new_mom_occ(wavefunc::MPS,p_count::Int64,p_end=4.0,p_start=0.0)
	num_sites = length(siteinds(wavefunc))
	dimension = dim(siteinds(wavefunc)[1])
	momenta = [p_start + (i-1)*(p_end - p_start)/(p_count-1) for i in 1:p_count]#[pi*i/(num_sites+1) for i in 1:p_count]
	mom_occ = zeros(p_count)*im
	for i in 1:p_count
		momentum = momenta[i]
		exp_vect = zeros(num_sites,num_sites)
		pos_occ = zeros(num_sites,num_sites) .* im
		for j in 1:num_sites
			exp_vect[:,j] = [-1 * cos(momentum*(l+j)) + cos(momentum*(j-l)) for l in 1:num_sites]
			#pos_occ[:,j] = [correlation_matrix(wavefunc,"Cr$l","Anh$j")[l,j] for l in 1:num_sites]
		end
		pos_occ = sum([correlation_matrix(wavefunc,"Cr$i1","Anh$j1") for i1 in 1:dimension-1 for j1 in 1:dimension-1])
		#return exp_vect,pos_occ
		#=fig = figure()
		imshow(exp_vect)
		colorbar()
		title("$momentum")
		=#
		mom_occ[i] = sum(exp_vect .* pos_occ) / (num_sites*(dimension-1))
	end
	return momenta,mom_occ
end
#
L = 6
nflavors = 4
org = [[2,2]]
psi = make_wavefunc(L,nflavors,org)
mm,oc = new_mom_occ(psi,100)
plot(mm,real.(oc),"-p")
yscale("log")
xscale("log")
#
#legend()
#


















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

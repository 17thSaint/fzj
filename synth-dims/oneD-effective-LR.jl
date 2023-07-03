include("fqh_effective.jl")
include("long-range-ttn.jl")


L = 10
nflavors = Int(L/2)
nbosons = Int(nflavors / 1.0)
t1 = 1
t2 = 0.5
U = 10
U1 = 4*t1^2/U
U2 = U1/2
conserve_qns = true

change = 0.0001
count = 5
#alpha  = 1/8
derivs = zeros(3,count)
errs = zeros(3,count)
alphas = [0.1 + (i-1)*(0.25-0.1)/(count-1) for i in 1:count]
all_wavefuncs = []

for i in 1:count
	alpha = alphas[i]
	if true
	Φ  = 2π*alpha
	Φ_right  = 2π*(alpha+change)
	Φ_left  = 2π*(alpha-change)

	model_paras = (t1 = t1, t2 = t2, Φ = Φ, U1 = U1, U2 = U2, L = L,
					   nflavors = nflavors)
	model_paras_left = (t1 = t1, t2 = t2, Φ = Φ_left, U1 = U1, U2 = U2, L = L,
					   nflavors = nflavors)
	model_paras_right = (t1 = t1, t2 = t2, Φ = Φ_right, U1 = U1, U2 = U2, L = L,
					   nflavors = nflavors)					   					   
	sidx = siteinds("ExtendedHardcore", L; conserve_qns = conserve_qns, nflavors = nflavors)
	H = MPO(hamiltonian(; model_paras...), sidx)
	H_left = MPO(hamiltonian(; model_paras_left...), sidx)
	H_right = MPO(hamiltonian(; model_paras_right...), sidx)
	println("Built Hams")
	states = make_states(L,nbosons,nflavors)
	psi0 = randomMPS(sidx, states)
	
	E, psi = dmrg(H, psi0; maxdim = [25,50, 100], nsweeps = 5, noise = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0])
	append!(all_wavefuncs,[psi])
	println("Done Center")
	
	E_left,psi_left = dmrg(H_left, psi; maxdim = 100, nsweeps = 3, noise = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0])
	
	E_right,psi_right = dmrg(H_right, psi; maxdim = 100, nsweeps = 3, noise = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0])
	
	end
	occ_mat_c = mps_get_occupancy(psi,L,nflavors)
	occ_mat_left = mps_get_occupancy(psi_left,L,nflavors)
	occ_mat_right = mps_get_occupancy(psi_right,L,nflavors)
	
	#bd_left = bulk_density(nothing; occ_mat=occ_mat_left)
	#bd_c = bulk_density(nothing; occ_mat=occ_mat_c)
	#bd_right = bulk_density(nothing; occ_mat=occ_mat_right)
	#println("Bulk Densities = $bd_left, $bd_c, $bd_right")
	for w in 1:2
		deriv12 = deriv_bulk_dens(nothing,nothing,-change,w; occ_mat1=occ_mat_left,occ_mat2=occ_mat_c)
		deriv23 = deriv_bulk_dens(nothing,nothing,-change,w; occ_mat1=occ_mat_c,occ_mat2=occ_mat_right)
		deriv13 = deriv_bulk_dens(nothing,nothing,-2*change,w; occ_mat1=occ_mat_left,occ_mat2=occ_mat_right)
		deriv = mean([deriv12,deriv23,deriv13])
		err = std([deriv12,deriv23,deriv13])
		derivs[w,i] = deriv
		errs[w,i] = err
		println("Streda Value = $deriv +/- $err")
	end
end

for i in 1:2
errorbar(alphas,derivs[i,:],yerr=[errs[i,:],errs[i,:]],label="$i")
end
legend()




























"fin"

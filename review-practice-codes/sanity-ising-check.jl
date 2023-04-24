using TTNKit

function ising_sweep(edge_length,num_sweeps,max_dim,j_strength=1.0,g_strength=0.1)
	L = (edge_length,edge_length)
	ind = TTNKit.siteinds("S=1/2",prod(L))
	net = TTNKit.BinaryNetwork(L, ind)
	lat = TTNKit.physical_lattice(net)

	states = fill("Up", TTNKit.number_of_sites(net))
	ttn = TTNKit.ProductTreeTensorNetwork(net, states)
	ttn = TTNKit.increase_dim_tree_tensor_network_zeros(ttn, maxdim = max_dim)
	ising = TTNKit.TransverseFieldIsing(J = j_strength, g = g_strength);
	tpo = TTNKit.Hamiltonian(ising, lat, mapping = TTNKit.hilbert_curve(lat));            
	proj_tpo = TTNKit.ProjectedTensorProductOperator(ttn, tpo)
	
	noise = 0.0
	eigsolve_tol = TTNKit.DEFAULT_TOL_DMRG
	eigsolve_krylovdim = TTNKit.DEFAULT_KRYLOVDIM_DMRG
	eigsolve_maxiter = TTNKit.DEFAULT_MAXITER_DMRG
	ishermitian = TTNKit.DEFAULT_ISHERMITIAN_DMRG
	eigsolve_which_eigenvalue = TTNKit.DEFAULT_WHICH_EIGENVALUE_DMRG
	func = (action,T) -> TTNKit.eigsolve(action, T, 1,
			            eigsolve_which_eigenvalue;
			            ishermitian=ishermitian,
			            tol=eigsolve_tol,
			            krylovdim=eigsolve_krylovdim,
			            maxiter=eigsolve_maxiter)

	sp = TTNKit.SimpleSweepHandler(ttn,proj_tpo,func,num_sweeps,[max_dim],[noise],TTNKit.NoExpander())
	TTNKit.sweep(ttn,sp);
	return ttn
end
#
edges = 2
swps = 3
dim = 2
all_ttns = []
for i in 1:10
	gsttn = ising_sweep(edges,swps,dim)
	append!(all_ttns,[gsttn])
end
#
for i in 1:length(all_ttns)
	for j in 1:length(all_ttns)
		if i == j
			continue
		else
			overlap = TTNKit.inner(all_ttns[i],all_ttns[j])
			if real(overlap) != 0.0
				println("$i, $j: ",overlap)
			end
		end
	end
end





















"fin"

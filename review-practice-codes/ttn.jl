using TTNKit,PyPlot

#=
Need to figure out how sweeps works
=#

function get_flattened_index(b_list)
	return sum(b_list .* [2^(length(b_list) - i) for i in 1:length(b_list)]) + 1
end

function get_xy(site_number,side_length)
	site_number 
end

edge_sites = 4
maxdim = 2
num_sweeps = 1
noise = 0.0

if false

square = TTNKit.BinaryNetwork((edge_sites,edge_sites), TTNKit.ITensorNode, "SpinHalf")
lat = TTNKit.physical_lattice(square)
num_sites = length(lat)
#sh = TTNKit.SimpleSweepHandler(chain,num_sweeps)

states = fill("Up",num_sites)
#states[Int(ceil(num_sites/2))] = "Dn"
ttn = TTNKit.ProductTreeTensorNetwork(square,states;orthogonalize=true)
ttn = TTNKit.increase_dim_tree_tensor_network_zeros(ttn, maxdim = maxdim)

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



iters = 10
nrgs = [0.0 for i in 1:iters]
jss = []
js, gs = -1.0, -2.0
ising = TTNKit.TransverseFieldIsing(J = js, g = gs)
tpo = TTNKit.Hamiltonian(ising,lat)
proj_tpo = TTNKit.ProjectedTensorProductOperator(ttn,tpo)

sp = TTNKit.SimpleSweepHandler(ttn,proj_tpo,func,num_sweeps,[maxdim],[noise],TTNKit.NoExpander())
TTNKit.sweep(ttn,sp);

end

adag = "S-"
reg_sites = (1)
prime_sites = (2)
ahat = "S+"

greenfunc = TTNKit.correlations(ttn,adag,ahat,reg_sites)


















"fin"

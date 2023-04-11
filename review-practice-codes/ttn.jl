using TTNKit

#=
Need to figure out how sweeps works
=#

num_sites = 4
maxdim = 4
num_sweeps = 2

chain = TTNKit.BinaryNetwork((num_sites,), TTNKit.ITensorNode, "SpinHalf")
lat = TTNKit.physical_lattice(chain)
#sh = TTNKit.SimpleSweepHandler(chain,num_sweeps)

states = fill("Up",num_sites)
#states[Int(ceil(num_sites/2))] = "Dn"
ttn = TTNKit.ProductTreeTensorNetwork(chain,states;orthogonalize=true)
ttn = TTNKit.increase_dim_tree_tensor_network_zeros(ttn, maxdim = maxdim)

js, gs = -1.0, -2.0
ising = TTNKit.TransverseFieldIsing(J = js, g = gs)
tpo = TTNKit.Hamiltonian(ising,lat)
proj_tpo = TTNKit.ProjectedTensorProductOperator(ttn,tpo)




















"fin"

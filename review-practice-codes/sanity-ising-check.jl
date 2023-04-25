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

function rewrite_inds(tensor,ref_tensor)
	num_layers = TTNKit.number_of_layers(TTNKit.network(tensor))
	for i in 1:num_layers
		num_tensors = length(tensor.data[i])
		for j in 1:num_tensors
			old_inds = TTNKit.inds(tensor.data[i][j])
			new_inds = TTNKit.inds(ref_tensor.data[i][j])
			if all([TTNKit.tags(old_inds[i1])==TTNKit.tags(new_inds[i1]) for i1 in 1:length(old_inds)])
				tensor.data[i][j] = TTNKit.replaceinds(tensor.data[i][j],old_inds,new_inds)
			else
				println("Index Tags Don't Match")
			end
		end
	end
	return tensor
end

function localinner(ttn1::TTNKit.TreeTensorNetwork{N, T}, old_ttn2::TTNKit.TreeTensorNetwork{N, T},nexttt=false) where{N<:TTNKit.BinaryNetwork,T}

    net = TTNKit.network(ttn1)
    
    ttn2 = rewrite_inds(old_ttn2,ttn1)

    elT = promote_type(eltype(ttn1), eltype(ttn2))
    # check in case if symmetric the Top node for qn correspondence
    if !(TTNKit.sectortype(net) == Int64)
        fl1 = flux(ttn1[TTNKit.number_of_layers(net), 1])
        fl2 = flux(ttn2[TTNKit.number_of_layers(net), 2])
        fl1 == fl2 || return zero(elT)
    end

    # contruct the network starting from the first layer upwards
    #ns = number_of_sites(net)
    
    phys_lat = TTNKit.physical_lattice(net)
    if nexttt
	    println(size(res))
    end
    res = map(phys_lat) do nd
    	TTNKit.delta(TTNKit.hilbertspace(nd), TTNKit.prime(TTNKit.hilbertspace(nd)))
    end


    for ll in TTNKit.eachlayer(net)
        nt = TTNKit.number_of_tensors(net,ll)
        res_new = Vector{T}(undef, nt)
        for pp in TTNKit.eachindex(net, ll)
            childs_idx = TTNKit.getindex.(TTNKit.child_nodes(net, (ll,pp)),2)
            tn1 = ttn1[ll,pp]
            tn2 = ttn2[ll,pp]
            rpre1 = res[childs_idx[1]]
            rpre2 = res[childs_idx[2]]
            if prod(size(rpre1)) > 2^5
            	println("Stop here: ",typeof(res))
            	return res
            end
            res_new[pp] = TTNKit._dot_inner(tn1, tn2, rpre1, rpre2)
        end
        res = res_new
    end
    # better exception
    length(res) == 1 || error("Tree Tensor Contraction don't leed to a single resulting tensor.")
    res = res[1]
    
    sres = TTNKit.ITensors.scalar(res)

    return sres
end

function make_ttn_chain_givenorg(org)
	edge_length = length(org)
	net = TTNKit.BinaryNetwork((edge_length,), TTNKit.ITensorNode, "SpinHalf")

	states = fill("Up",edge_length)
	for i in 1:edge_length
		if org[i] == 0
			states[i] = "Dn"
		end
	end
	ttn = TTNKit.ProductTreeTensorNetwork(net,states;orthogonalize=true)
	return ttn
end

#=
org1 = [0,1,0,1]
org2 = [1,0,1,0]
ttn1 = make_ttn_chain_givenorg(org1)
ttn2 = make_ttn_chain_givenorg(org2)
rez = localinner(ttn1,ttn2)
=#

#
edges = 2
swps = 3
dim = 2
howmany = 5
#
#
next = false
for i in 1:howmany
	ttn_one = ising_sweep(edges,swps,dim)	
	ttn_two = ising_sweep(edges,swps,dim)
	overlap = localinner(ttn_one,ttn_two)
	if typeof(overlap) != ComplexF64
		println("Probably broken")
	elseif real(overlap) != 0.0
		println(overlap)
	end
end
#

#this = localinner(all_ttns[2],all_ttns[2])


















"fin"

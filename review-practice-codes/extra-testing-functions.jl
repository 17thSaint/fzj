#include("../review-practice-codes/ttn.jl")


function simple_boson_network(num_layers::Int, conserve_qns::Bool, max_occ::Int)
    return TTNKit.BinaryRectangularNetwork(num_layers, TTNKit.ITensorNode, "Boson";conserve_qns=conserve_qns,dim=max_occ+1)
end

function make_randomconfig_ttn(net::TTNKit.AbstractNetwork, particle_count::Int, max_occ::Int; kwargs...)
    layers = TTNKit.number_of_layers(net)
    site_count = 2^layers
    states::Vector{String} = fill_states(particle_count,site_count,max_occ)
	ttn::TTNKit.TreeTensorNetwork = TTNKit.ProductTreeTensorNetwork(net,states)
    return ttn
end

function make_parton_ttn(net::TTNKit.AbstractNetwork, particle_count::Int, max_dim::Int, max_occ::Int)
    num_sites = TTNKit.number_of_sites(net)
    states::Vector{String} = fill("0", num_sites)
	old_ttn::TTNKit.TreeTensorNetwork = TTNKit.ProductTreeTensorNetwork(net,states)
	ttn = initialize_ttn(old_ttn,max_dim,particle_count)

    return ttn
end

#= This makes a 2x2 hardcore lattice with 2 particles at precise positions and a roughly delocalized mixed state
layers = 2
max_occ = 1
particles = 2

net = simple_boson_network(layers, true, max_occ)
psi1 = make_randomconfig_ttn(net,particles,max_occ)
psi2 = make_randomconfig_ttn(net,particles,max_occ)

occs1 = get_occupancy(psi1; plot_title="Psi 1")
occs2 = get_occupancy(psi2; plot_title="Psi 2")

mixedstate = make_parton_ttn(net,particles,10,max_occ)
occs_mixed = get_occupancy(mixedstate; plot_title="Mixed state")
=#


























"fin"
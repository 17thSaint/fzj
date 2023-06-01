using NBInclude,TTNKit,Test

@nbinclude("parton-model-syms.ipynb")

do_all = true

if do_all | false
@testset "ttn mapping" begin

	two_edge_site_exp = [1,2,3,4]
	two_edge_site_calc = ttn_2d_mapping(2)
	@test two_edge_site_exp == two_edge_site_calc
	
	four_edge_site_exp = [1,2,5,6,3,4,7,8,9,10,13,14,11,12,15,16]
	four_edge_site_calc = ttn_2d_mapping(4)
	@test four_edge_site_exp == four_edge_site_calc

end;
end

if do_all | false
@testset "Localize one boson at given site" begin

	layers = 4
	edge_length = Int(sqrt(2^layers))
	net = TTNKit.BinaryRectangularNetwork(layers, TTNKit.ITensorNode, "Boson"; conserve_qns = true)
	lat = TTNKit.physical_lattice(net)
	
	states = fill("0",Int(2^layers))
	ttn = TTNKit.ProductTreeTensorNetwork(net,states)
	
	wf_coefs = zeros(size(lat))
	
	localized_site = (1,1)
	wf_coefs[localized_site[1],localized_site[2]] = 1.0

	ttn = patron_application!(ttn,wf_coefs,"Adag";maxdim = 4)

	@test isapprox(real.(TTNKit.expect(ttn,"N"))[localized_site[1],localized_site[2]],1.0,atol=10^-5)
	
end;
end

function TransverseFieldIsing2d(J,g,lat)
	ampo = OpSum()
	for p in TTNKit.coordinates(lat)
		ampo += g, "Z", p
	end
	for bond in TTNKit.nearest_neighbours(lat, collect(eachindex(lat)))
        	b1 = TTNKit.coordinate(lat, first(bond))
        	b2 = TTNKit.coordinate(lat, last(bond))
        	ampo += J, "X", b1, "X", b2
       	end
	return ampo
end

function get_energy_pol(tdvp_rez)
	topPos = (TTNKit.number_of_layers(TTNKit.network(tdvp_rez.ttn)), 1)
	n_sites = TTNKit.number_of_sites(TTNKit.network(tdvp_rez.ttn))
	action = TTNKit.∂A(tdvp_rez.pTPO, topPos)
	T_rez = tdvp_rez.ttn[topPos]
	actionT = action(T_rez)
	energy = real(TTNKit.ITensors.scalar(dag(T_rez)*actionT))
	pol = real(sum(TTNKit.expect(tdvp_rez.ttn, "Z"))/n_sites)
	return energy,pol
end

if do_all | false
@testset "Time evolution Transverse Field Ising 2x2 with Parton Initialization" begin

        polarization_exact = 0.9501373043948331
	(J, g) = (-1, -2)
	(tmax,dt) = (1, 1e-1)
	
	L = (2,2)
        ind = TTNKit.ITensors.siteinds("S=1/2", prod(L))
        net = TTNKit.BinaryNetwork(L, ind)
        lat = TTNKit.physical_lattice(net)

        states = fill("0", TTNKit.number_of_sites(net))
        ttn = TTNKit.ProductTreeTensorNetwork(net, states)
        
	wf_coefs = zeros(L) .+ 1/sqrt(prod(L))
        ttn = patron_application!(ttn,wf_coefs,"projUp";maxdim = 16)

        ising = TransverseFieldIsing2d(J, g, lat);
        tpo = TTNKit.TPO(ising, lat)

	rez = TTNKit.tdvp(ttn, tpo, finaltime = tmax, timestep = dt)
	energy, pol = get_energy_pol(rez)
	@test isapprox(pol,polarization_exact,atol=10^-5)
	@test energy ≈ -8.
        
end;
end

























"fin"

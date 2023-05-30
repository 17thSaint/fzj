using NBInclude,TTNKit,Test

@nbinclude("parton-model-syms.ipynb")

do_all = false

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

#if do_all | true
#@testset "Make one boson at random site" begin

	layers = 4
	edge_length = Int(sqrt(2^layers))
	net = TTNKit.BinaryRectangularNetwork(layers, TTNKit.ITensorNode, "Boson"; conserve_qns = true)
	lat = TTNKit.physical_lattice(net)
	
	states = fill("0",Int(2^layers))
	ttn = TTNKit.ProductTreeTensorNetwork(net,states)
	
	wf_coefs = zeros(size(lat))
	rand_site = (rand((1:edge_length)),rand((1:edge_length)))
	wf_coefs[rand_site[1],rand_site[2]] = 1.0

	ttn = patron_application!(ttn,wf_coefs,"Adag";maxdim = 4)
	println(rand_site)
	display(real.(TTNKit.expect(ttn,"N")))
#	@test isapprox(real(TTNKit.expect(ttn,"N",rand_site)),1.0,atol=10^-5)
	
#end;
#end



























"fin"

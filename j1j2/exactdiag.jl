include("j1j2.jl")
using ITensors


function get_exact_j1j2ham(j2,layers=4; kwargs...)
	if_periodic = get(kwargs, :if_periodic, true)
	j1 = get(kwargs, :j1, 1.0)
	
	net = TTNKit.BinaryRectangularNetwork(layers, TTNKit.ITensorNode, "Qubit")
	lat = TTNKit.physical_lattice(net)
	resulting_ham = []
	
	if j1 != 0.0
		j1_term = OpSum()
		for (s1,s2) in TTNKit.nearest_neighbours(lat,collect(1:TTNKit.number_of_sites(lat)); periodic=if_periodic)
			
			j1_term += (j1,"Sx",s1,"Sx",s2)
			j1_term += (j1,"Sy",s1,"Sy",s2)
			j1_term += (j1,"Sz",s1,"Sz",s2)
		end
		append!(resulting_ham,[j1_term])
	end
	
	if j2 != 0.0
		j2_term = OpSum()
		for (s1,s2) in find_next_nearest_neighbors(lat,if_periodic)
			
			j2_term += (j2,"Sx",s1,"Sx",s2)
			j2_term += (j2,"Sy",s1,"Sy",s2)
			j2_term += (j2,"Sz",s1,"Sz",s2)
		end
		append!(resulting_ham,[j2_term])
	end
	
	if length(resulting_ham) > 1
		return sum(resulting_ham),lat
	else
		return resulting_ham[1],lat
	end
end

j2 = 0.0
layers = 2
if layers > 2
	ITensors.disable_warn_order()
end
sumham = get_exact_j1j2ham(j2,layers)
ham_mpo = MPO(sumham[1],siteinds(sumham[2]))
ham_mat = round.(reshape(Array(prod(ham_mpo),siteinds(ham_mpo)),2^(2^layers),2^(2^layers)),digits=6)




































"fin"

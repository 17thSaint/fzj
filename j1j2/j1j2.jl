include("../review-practice-codes/ttn.jl")

function find_next_nearest_neighbors(lat,if_periodic=false)
	nn_neighbors = []
	xlen,ylen = size(lat)
	shift = if_periodic ? 0 : xlen
	for s in collect(1:Int(xlen*ylen - shift))
		if (s-1) % xlen == 0
			append!(nn_neighbors,[(s,mod(s+xlen+1,Int(xlen*ylen)))])
			if if_periodic
				append!(nn_neighbors,[(s,mod(s+2*xlen-1,Int(xlen*ylen)))])
			end
		elseif s % xlen == 0
			append!(nn_neighbors,[(s,mod(s+xlen-1,Int(xlen*ylen)))])
			if if_periodic
				append!(nn_neighbors,[(s,mod(s+1,Int(xlen*ylen)))])
			end
		else
			append!(nn_neighbors,[(s,mod(s+xlen+1,Int(xlen*ylen)))])
			append!(nn_neighbors,[(s,mod(s+xlen-1,Int(xlen*ylen)))])
		end
	end
	return nn_neighbors
end

function get_j1j2_hamilt(net,j2,j1=1; kwargs...)
	if_periodic = get(kwargs, :if_periodic, true)
	
	lat = TTNKit.physical_lattice(net)
	resulting_ham = []
	
	if j1 != 0.0
		j1_term = TTNKit.OpSum()
		for (s1,s2) in TTNKit.nearest_neighbours(lat,collect(1:TTNKit.number_of_sites(lat)); periodic=if_periodic)
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)
			
			j1_term += (j1,"Sx",s1_coord,"Sx",s2_coord)
			j1_term += (j1,"Sy",s1_coord,"Sy",s2_coord)
			j1_term += (j1,"Sz",s1_coord,"Sz",s2_coord)
		end
		append!(resulting_ham,[j1_term])
	end
	
	if j2 != 0.0
		j2_term = TTNKit.OpSum()
		for (s1,s2) in find_next_nearest_neighbors(lat,if_periodic)
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)
			
			j2_term += (j2,"Sx",s1_coord,"Sx",s2_coord)
			j2_term += (j2,"Sy",s1_coord,"Sy",s2_coord)
			j2_term += (j2,"Sz",s1_coord,"Sz",s2_coord)
		end
		append!(resulting_ham,[j2_term])
	end
	
	if length(resulting_ham) > 1
		return sum(resulting_ham)
	else
		return resulting_ham[1]
	end
end

function get_spin_texture(ttn; kwargs...)
	exp_occ = real.(TTNKit.expect(ttn,"Sz"))
	
	if_plot = get(kwargs, :if_plot, true)
	
	if_plot ? plot_spin_texture(exp_occ; kwargs...) : nothing
		
	return exp_occ
end

using PyPlot
function plot_spin_texture(exp_occ; kwargs...)
	data_dict = get(kwargs, :data_dict, nothing)
	if !isnothing(data_dict)
		exp_occ = data_dict["vals"]
	end
	fig = figure()
	imshow(exp_occ)
	colorbar()
	plot_title = get(kwargs, :plot_title, "")
	title_string = "Spin Texture, " * plot_title
	title(title_string)
	
	if get(kwargs, :if_save_fig, false)
		location = get(kwargs, :location, pwd())
		fig_name = get(kwargs, :name, "occs")
		fig_name = check_plot_label(fig_name,"occs")
		save_figure(fig_name; kwargs...)
	end
	return
end



layers = 4
num_parts = 4
syms = true
max_occ = 1
mdim = 100
nsweeps = 5
if_periodic = true

net = TTNKit.BinaryRectangularNetwork(layers, TTNKit.ITensorNode, "Qubit")
j2 = 0.4

hamj1j2 = get_j1j2_hamilt(net,j2; if_periodic = if_periodic)

model_paras = (max_dim = mdim, if_save_data = false, num_sweeps = nsweeps, sweep_type = "dmrg", syms = syms, ttn_net = net, ham_op = hamj1j2, if_fermion = true, part_type = "Qubit")

og_ttn, og_ham, dmsp = find_ground_state(layers,num_parts,nothing; model_paras...)

get_spin_texture(dmsp.ttn; plot_title = "J2 = $j2")








































"fin"

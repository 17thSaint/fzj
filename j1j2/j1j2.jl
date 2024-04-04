using Pkg
Pkg.activate(".")
include("../review-practice-codes/ttn.jl")
using Statistics,MKL
using PyPlot

function find_next_nearest_neighbors(layers,if_periodic=false)
	nn_neighbors = []

	if layers % 2 == 0
		lat_size = (Int(sqrt(2^layers)),Int(sqrt(2^layers)))
	else
		lat_size = (Int(sqrt(2^(layers-1))),Int(sqrt(2^(layers+1))))
	end
	lat = TTNKit.SimpleLattice(lat_size,TTNKit.ITensorNode,"Qubit")

	xlen,ylen = lat_size
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

function get_j1j2_hamilt(layers,j2,j1=1; kwargs...)
	if_periodic::Bool = get(kwargs, :if_periodic, false)
	if_mag_orientation::Bool = get(kwargs, :if_mag_orientation, false)
	mag_direction::Float64 = get(kwargs, :mag_direction, 1.0)
	
	if layers % 2 == 0
		lat_size = (Int(sqrt(2^layers)),Int(sqrt(2^layers)))
	else
		lat_size = (Int(sqrt(2^(layers+1))),Int(sqrt(2^(layers-1))))
	end
	lat = TTNKit.SimpleLattice(lat_size,TTNKit.ITensorNode,"Qubit")
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
		for (s1,s2) in find_next_nearest_neighbors(layers,if_periodic)
			s1_coord = TTNKit.coordinate(lat,s1)
			s2_coord = TTNKit.coordinate(lat,s2)
			
			j2_term += (j2,"Sx",s1_coord,"Sx",s2_coord)
			j2_term += (j2,"Sy",s1_coord,"Sy",s2_coord)
			j2_term += (j2,"Sz",s1_coord,"Sz",s2_coord)
		end
		append!(resulting_ham,[j2_term])
	end

	if if_mag_orientation
		mag_term = TTNKit.OpSum()
		mag_term -= (mag_direction*1000.0,"Sz",(1,1))
		append!(resulting_ham,[mag_term])
	end

	if length(resulting_ham) > 1
		return sum(resulting_ham)
	else
		return resulting_ham[1]
	end
end

function get_spin_texture(ttn::TTNKit.TreeTensorNetwork; kwargs...)
	if_sign = get(kwargs, :if_sign, false)
	
	exp_occ = real.(TTNKit.expect(ttn,"Sz"))
	if_sign ? exp_occ = sign.(exp_occ) : nothing
	
	if_plot = get(kwargs, :if_plot, true)
	
	if_plot ? plot_spin_texture(exp_occ; kwargs...) : nothing
		
	return exp_occ
end

function plot_spin_texture(exp_occ; kwargs...)
	data_dict = get(kwargs, :data_dict, nothing)
	if !isnothing(data_dict)
		exp_occ = data_dict["vals"]
	end
	if_phase = get(kwargs, :if_phase, false)
	if if_phase
		phase = pitch_vec = get_pitch_vector(nothing; exp_occ = exp_occ)[2]
	else
		phase = nothing
	end
	fig = figure()
	imshow(exp_occ)
	colorbar()
	plot_title = get(kwargs, :plot_title, "")
	title_string = if_phase ? "Spin Texture, " * phase * ", " * plot_title : "Spin Texture, " * plot_title
	title(title_string)
	
	if get(kwargs, :if_save_fig, false)
		location = get(kwargs, :location, pwd())
		fig_name = get(kwargs, :name, "occs")
		fig_name = check_plot_label(fig_name,"occs")
		save_figure(fig_name; kwargs...)
	end
	return
end

function get_pitch_vector(wavefunc; kwargs...)
	exp_occ = get(kwargs, :exp_occ, nothing)
	if isnothing(exp_occ)
		texture = get_spin_texture(wavefunc; if_sign = true, if_plot = false)
	else
		if abs(exp_occ[1,1]) != 1.0
			texture = sign.(exp_occ)	
		else
			texture = exp_occ
		end
	end
	if isapprox(mean(texture[1,:]),0.0;atol=0.1) && isapprox(mean(texture[:,1]),0.0;atol=0.1)
		#println("Neel")
		return (pi,pi),"Neel"
	else
		if isapprox(sum(texture[1,:] .> 0.0) / sum(texture[1,:] .< 0.0),1.0;atol=0.1)
			#println("Vertical")
			return (0,pi),"Vertical"
		elseif isapprox(sum(texture[:,1] .> 0.0) / sum(texture[:,1] .< 0.0),1.0;atol=0.1)
			#println("Horizontal")
			return (pi,0),"Horizontal"
		else
			println("Unknown Phase")
			return nothing,nothing
		end
	end
end

function allsitescorrelation(wavefunc,op::String)
	num_sites = TTNKit.number_of_sites(TTNKit.network(wavefunc))
	corrmat = zeros(num_sites,num_sites).*im
	for i in 1:num_sites
		for j in 1:i
			corrmat[i,j] = TTNKit.correlation(wavefunc,op,op,i,j)
			corrmat[j,i] = corrmat[i,j]
		end
	end
	return corrmat
end

function distance_btw(point1,point2,lattice_size; kwargs...)
	if_periodic = get(kwargs, :if_periodic, true)
	
	dx = abs(point1[1] - point2[1])
	dy = abs(point1[2] - point2[2])
	if if_periodic
		dx = min(dx, lattice_size[1] - dx)
		dy = min(dy, lattice_size[2] - dy)
	end
	
	return sqrt(dx^2 + dy^2)
end

function distance_correlation(wavefunc,op::String; kwargs...)
	if_plot = get(kwargs, :if_plot, true)
	starting_site = get(kwargs, :starting_site, 1)
	
	lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))
	num_sites = TTNKit.number_of_sites(TTNKit.network(wavefunc))
	corrmat = []
	distances = []
	for i in 1:num_sites
		dist_btw = distance_btw(TTNKit.coordinate(lat,starting_site),TTNKit.coordinate(lat,i),size(lat); kwargs...)
		if dist_btw != 0.0
			append!(corrmat,[abs(TTNKit.correlation(wavefunc,op,op,starting_site,i))])
			append!(distances,[dist_btw])
		end
	end

	if_plot ? plot_distance_correlation(distances,corrmat,op; kwargs...) : nothing
	
	return distances,corrmat
end

function plot_distance_correlation(distances,corrmat,op::String; kwargs...)
	label_string = get(kwargs, :label_string, nothing)
	
	if isnothing(label_string)
		fig = figure()
		plot(distances,corrmat,"p")
	else
		plot(distances,corrmat,"p",label=label_string)
		legend()
	end
	xlabel("Distances")
	ylabel("$(op)_0 $(op)_r")
	xscale("log")
	yscale("log")
	plot_title = get(kwargs, :plot_title, "")
	title_string = "$op Distance Correlation, " * plot_title
	title(title_string)
end

function spinspin_structure_factor(wavefunc::TTNKit.TreeTensorNetwork; kwargs...)
	lat = TTNKit.physical_lattice(TTNKit.network(wavefunc))
	num_sites = prod(size(lat))
	corr_mat = allsitescorrelation(wavefunc,"Sx") .+ allsitescorrelation(wavefunc,"Sy") .+ allsitescorrelation(wavefunc,"Sz")
	pitch_vector = get(kwargs, :pitch_vector, get_pitch_vector(wavefunc)[1])
	exp_mat = zeros(size(corr_mat)).*im
	for i in 1:size(exp_mat)[1]
		for j in 1:i
			phase_part = dot(pitch_vector,TTNKit.coordinate(lat,i) .- TTNKit.coordinate(lat,j))
			exp_mat[i,j] = exp(im * phase_part)
			exp_mat[j,i] = exp(-im * phase_part)
		end
	end
	return real(sum(exp_mat .* corr_mat) / (num_sites*(num_sites+2)))#, exp_mat
end

if true

#ll = 6
#j_two = 0.0
#bonddims = [100,150]
#for bonddim in bonddims
#mag_dir = idx == 1 ? 1.0 : -1.0
#params_dict = Dict([("layers",ll),("j2",j_two),("mdim",bonddim),("if_save_data",false),("if_mag_orientation",false),("if_periodic",false)])
params_dict = make_args_dict(ARGS)
open_cores = get(params_dict, "open_cores", "all")
if typeof(open_cores) != String
	BLAS.set_num_threads(open_cores)	
	display(BLAS.get_config())
end

layers = get(params_dict, "layers", 4)
num_parts = Int((2^layers)/2)
syms = get(params_dict,"syms",true)
if_periodic = get(params_dict,"if_periodic",false)
j2 = get(params_dict,"j2",0.0)
if_mag_orientation = get(params_dict,"if_mag_orientation",false)
mag_direction = get(params_dict,"mag_direction",1.0)

mdim = get(params_dict,"mdim",100)
num_sweeps = get(params_dict,"num_sweeps",100)

dataloc = get(params_dict, "dataloc", get_folder_location("cluster-data/j1j2"))
if_save_data = get(params_dict,"if_save_data",true)
if_cluster = any([occursin("local",pwd()),occursin("Local",pwd()),occursin("geraghty",pwd())])
if_continuous_saving = get(params_dict,"if_continuous_saving",if_cluster)
if_save_data ? nothing : if_continuous_saving = false
if_densmat = get(params_dict,"if_densmat",false)
if_gpu = get(params_dict,"if_gpu",false)
nrgtol = get(params_dict,"nrgtol",1e-4)
cutoff = get(params_dict,"cutoff",1e-10)

if_redo = get(params_dict,"if_redo",false)

noise = get(params_dict,"noise",0.0)
expander_val = get(params_dict,"expander_val",1.0)
expander = TTNKit.DefaultExpander(expander_val)


filename_dict = Dict([("layers",layers),("j2",j2),("mdim",mdim),("if_mag_orientation",if_mag_orientation),("if_periodic",if_periodic)])
datafile_name = make_parameters_filename(filename_dict)

model_paras = (mdim=mdim,
				j2=j2,
				nrgtol=nrgtol,
				cutoff=cutoff,
				if_periodic=if_periodic,
				if_mag_orientation=if_mag_orientation,
				mag_direction=mag_direction,
				if_densmat=if_densmat,
				location=dataloc,
				if_save_data=if_save_data,
				if_continuous_saving=if_continuous_saving,
				name = "ttn-"*datafile_name,
				if_gpu=if_gpu,
				noise=noise,
				num_sweeps=num_sweeps,
				sweep_type="dmrg",
				syms=syms,
				part_type = "Qubit",
				particles=num_parts)

metadata_dict = merge(named_tuple_to_dict(model_paras),filename_dict)

if_exists,found_data = check_data_exists(filename_dict,"ttn"; location=dataloc,output_level=false)

if if_exists
	println("Found Data")
	data,metadata = found_data
	wavefunc = data["ttn"]
	rezobs = metadata["observer"]
	
	if if_redo 
		starting = time()
		net = TTNKit.BinaryRectangularNetwork(layers, TTNKit.ITensorNode, "Qubit", conserve_qns=syms)
		hamj1j2 = get_j1j2_hamilt(layers,j2; model_paras...)

		og_ttn, hamilt, dm_sp, rezobs, runtime = find_ground_state(layers,num_parts; ttn_net=net,ham_op=hamj1j2,model_paras...,if_redo=true,metadata=merge(metadata_dict,Dict([("ham",hamj1j2),("net",net)])))
		total_time = time() - starting
		wavefunc = dm_sp.ttn
		println("Total running time: $total_time")
	end

else		
	starting = time()
	net = TTNKit.BinaryRectangularNetwork(layers, TTNKit.ITensorNode, "Qubit", conserve_qns=syms)
	hamj1j2 = get_j1j2_hamilt(layers,j2; model_paras...)
	display(hamj1j2)

	og_ttn, hamilt, dm_sp, rezobs, runtime = find_ground_state(layers,num_parts; ttn_net=net,ham_op=hamj1j2,model_paras...,metadata=merge(metadata_dict,Dict([("ham",hamj1j2),("net",net)])))
	total_time = time() - starting
	println("Total running time: $total_time")

	wavefunc = dm_sp.ttn
end

#final_energy = if_mag_orientation ? (rezobs.nrg[end] + 1000*mag_direction*TTNKit.expect(wavefunc,"Sz",(1,1))) / 2^layers : rezobs.nrg[end] / 2^layers

#text = get_spin_texture(wavefunc; if_sign = false,plot_title = "J2 = $j2, NRG = $(round(real(final_energy),digits=5))")

#append!(all_wavefuncs,[wavefunc])

#println("Energy per site = ",round(final_energy,digits=5))
#scatter(1/TTNKit.maxlinkdim(wavefunc),final_energy,c="b")
#xscale("log")
#end

end

#=w1_overlap = TTNKit.inner(all_wavefuncs[1],all_wavefuncs[3])
w2_overlap = TTNKit.inner(all_wavefuncs[2],all_wavefuncs[3])
println("Overlap with Up Corner = ",w1_overlap," Overlap with Down Corner = ",w2_overlap)#

which_wavefunc = 3
s1 = 1
s2 = 2
net = TTNKit.network(all_wavefuncs[which_wavefunc])
ttnc1 = TTNKit.copy(all_wavefuncs[which_wavefunc])
ttnc2 = TTNKit.copy(all_wavefuncs[which_wavefunc])

idx_1 = TTNKit.siteinds(net)[s1]
O1 = TTNKit.convert_cu(TTNKit.op("Sz",idx_1),all_wavefuncs[which_wavefunc][(1,1)])
ch_pos1 = (0,s1)
parent_pos1 = TTNKit.parent_node(net,ch_pos1)
TTNKit.move_ortho!(ttnc1,parent_pos1)
T1 = ttnc1[parent_pos1]
first_application = TTNKit.noprime(O1*T1)
ttnc1[parent_pos1] = first_application

idx_2 = TTNKit.siteinds(net)[s2]
O2 = TTNKit.convert_cu(TTNKit.op("Sz",idx_2),all_wavefuncs[which_wavefunc][(1,1)])
ch_pos2 = (0,s2)
parent_pos2 = TTNKit.parent_node(net,ch_pos2)
TTNKit.move_ortho!(ttnc1,parent_pos2)
T2 = ttnc1[parent_pos2]
second_application = TTNKit.noprime(O2*first_application)

TTNKit.move_ortho!(ttnc2,parent_pos2)

T = ttnc2[parent_pos2]
res = dot(T, second_application)

println("Result = ",res)
=#

































"fin"

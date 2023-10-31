include("../review-practice-codes/ttn.jl")
using Statistics

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
	if_sign = get(kwargs, :if_sign, false)
	
	exp_occ = real.(TTNKit.expect(ttn,"Sz"))
	if_sign ? exp_occ = sign.(exp_occ) : nothing
	
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


#
layers = 4
num_parts = Int((2^layers)/2)
syms = true
max_occ = 1
mdim = 20
nsweeps = 3
if_periodic = true
dataloc = "../cluster-data/j1j2/"

net = TTNKit.BinaryRectangularNetwork(layers, TTNKit.ITensorNode, "Qubit")
#j2 = 0.0

j2_count = 10
j2s = [(i-1)*2/(j2_count-1) for i in 1:j2_count]
s2_vert = [0.0 for i in 1:j2_count]
s2_neel = [0.0 for i in 1:j2_count]
s2_horiz = [0.0 for i in 1:j2_count]
for i in  1:j2_count
j2 = j2s[i]

hamj1j2 = get_j1j2_hamilt(net,j2; if_periodic = if_periodic)

model_paras = (max_dim = mdim, if_save_data = false, num_sweeps = nsweeps, sweep_type = "dmrg", syms = syms, ttn_net = net, ham_op = hamj1j2, part_type = "Qubit")

og_ttn, og_ham, dmsp = find_ground_state(layers,num_parts,nothing; model_paras...)

text = get_spin_texture(dmsp.ttn; if_sign = true,plot_title = "J2 = $j2")
s2_vert[i] = spinspin_structure_factor(dmsp.ttn; pitch_vector = (0,pi))
s2_horiz[i] = spinspin_structure_factor(dmsp.ttn; pitch_vector = (pi,0))
s2_neel[i] = spinspin_structure_factor(dmsp.ttn; pitch_vector = (pi,pi))
println("V = ",s2_vert[i],", H = ",s2_horiz[i],", N = ",s2_neel[i])
end
fig = figure()
plot(j2s,s2_vert,"-p",label="Vert")
plot(j2s,s2_horiz,"-p",label="Horiz")
plot(j2s,s2_neel,"-p",label="Neel")
legend()
xlabel("J2")
ylabel("S^2")
#







































"fin"

#####################################################
#=

This file contains the model parameters (get_fpeps_model_paras), Hamiltonian (long_range_HH_ham) and the
column-exclusion initial state (column_exclusion_peps) for fPEPS

Depends on (all loaded by fpeps/load-fpeps.jl):
	other-funcs/basic-2d-stuff.jl (coordinate, find_dist, check_fluxes, long_range_scaling, get_hopping_strengths)
	other-funcs/data-storage-funcs.jl (get_folder_location, default dataloc)
	ITensors (OpSum), QuantumNaturalfPEPS (PEPS)

All parameters are read from a String-keyed dict, axes are x (1st index, 1:Lx) and
y (2nd index, 1:Ly), matching hilbert = siteinds("S=1/2", Lx, Ly). get_fpeps_model_paras keeps the
inputs (lr, ts, filling, scaling inputs) next to what is derived from them, so it is also the saved
metadata and running it again on that metadata gives the same model. Open boundaries only

Hardcore bosons on S=1/2 sites: Up = occupied, Dn = empty, so b† = S+, b = S-, n = ProjUp

=#
######################################################


# Peierls hopping from site along axis (1 = x, 2 = y), Landau gauge with the phase on the
# flux_direction hoppings: x hops pick up exp(2πiα(y-1)), y hops exp(2πiα(x-1))
function get_hopping_coeff(site,axis::Int,t_x::Float64,t_y::Float64,alpha::Float64,flux_direction::String)
	if axis == 1
		coeff = -t_x
		flux_direction == "x" ? coeff *= exp(im*2*pi*alpha*(site[2]-1)) : nothing
	else
		coeff = -t_y
		flux_direction == "y" ? coeff *= exp(im*2*pi*alpha*(site[1]-1)) : nothing
	end
	return coeff
end

# [#Up = occupied, #Dn = empty] counts of the N particle sector, the max_counts of get_sample / generate_Oks_and_Eks
# (sz_sector_counts(peps, Q) would fix total Sz = Q, i.e. N - Lx*Ly/2)
particle_sector_counts(Lx::Int, Ly::Int, N::Int) = [N, Lx*Ly - N]

# inputs of each scaling with the defaults of long_range_scaling (nothing = required), as (key, default) pairs.
# The one place fpeps lists them: get_fpeps_model_paras keeps them in the metadata, make_fpeps_filename in the filename
function fpeps_scaling_inputs(scaling::String, interaction_axis_length::Int)
	all_inputs = Dict([("flat",[]),("gaussian",[("sigma",1.0)]),("exp",[("corr_length",interaction_axis_length)]),
		("rydberg",[("blockade_radius",1.0)]),("dd",[("magnetic_spacing",1.0)]),("powerlaw",[("magnetic_spacing",nothing),("exponent",nothing)])])
	haskey(all_inputs, scaling) || error("ULR Scaling Type Not Recognized: $scaling")
	return all_inputs[scaling]
end

function get_fpeps_model_paras(params_dict::Dict{String,Any})
	
	Lx,Ly = params_dict["Lx"],params_dict["Ly"]
	L = (Lx,Ly)
	num_sites = Lx*Ly
	outputlevel = get(params_dict, "outputlevel", 1)
	outputlevel > 0 ? println("Lx = ",Lx,", Ly = ",Ly) : nothing

	bond_dim = get(params_dict, "bond_dim", 5)

	# lattice
	if_periodic_x = get(params_dict, "if_periodic_x", false)
	if_periodic_y = get(params_dict, "if_periodic_y", false)
	if_per = (if_periodic_x,if_periodic_y)

	# hopping and flux
	t_strength = get(params_dict, "ts", 1.0)
	hopping_anisotropy = get(params_dict, "hopping_anisotropy", 1.0)
	t_x,t_y = get_hopping_strengths(t_strength, hopping_anisotropy)
	twist_angle = get(params_dict, "twist_angle", [0.0,0.0])
	particles = get(params_dict, "particles", Int(ceil(Lx/2)))
	filling = get(params_dict, "filling", 0.5)
	if haskey(params_dict, "alpha")
		alpha = Float64(params_dict["alpha"])
	else
		alpha = filling == 0.0 ? 0.0 : particles/(filling*(Lx - !if_periodic_x)*(Ly - !if_periodic_y))
	end
	flux_direction = get(params_dict, "flux_direction", "y")
	if if_periodic_y && !if_periodic_x
		flux_direction = "y"
	elseif if_periodic_x && !if_periodic_y
		flux_direction = "x"
	end
	if_check_fluxes = get(params_dict, "if_check_fluxes", true)
	if if_check_fluxes && alpha != 0.0
		flux_direction = check_fluxes(alpha,Lx,Ly,if_periodic_x,if_periodic_y,flux_direction; output_level=outputlevel)
	end

	# interactions, which_dir is the axis the long-range interaction acts along ("x", "y" or "both")
	onsite_strength::Float64 = get(params_dict, "onsite_strength", 0.0)
	scaling = get(params_dict, "scaling", "flat")
	which_dir = scaling == "rydberg" ? "both" : get(params_dict, "which_dir", "y")
	interaction_axis_length = which_dir == "x" ? Lx : which_dir == "y" ? Ly : max(Lx,Ly)
	scaling_distance = get(params_dict, "lr", 0)
	scaling_distance == "all" ? scaling_distance = interaction_axis_length - 1 : nothing
	scaling_inputs = fpeps_scaling_inputs(scaling, interaction_axis_length)
	missing_inputs = [k for (k,d) in scaling_inputs if isnothing(d) && !haskey(params_dict, k)]
	isempty(missing_inputs) || error("scaling $scaling needs $(join(missing_inputs, ", ")) in params_dict")
	scaling_paras = Dict{String,Any}([(k, isnothing(d) ? params_dict[k] : get(params_dict, k, d)) for (k,d) in scaling_inputs])
	scaling_kwargs = Dict(Symbol(k) => v for (k,v) in params_dict)
	long_range_strengths = long_range_scaling(scaling_distance,interaction_axis_length,onsite_strength; scaling_kwargs...)
	outputlevel > 0 ? display(long_range_strengths) : nothing

	# on-site potentials
	chem_strength = get(params_dict, "chem_strength", 0.0)
	if_pinning = get(params_dict, "if_pinning", false)
	pinning_strength::Float64 = get(params_dict, "pinning_strength", 1e-4)

	# saving stuff
	dataloc = get(params_dict, "dataloc", get_folder_location("cluster-data/synth-dims/fpeps"))
	if_save_data = get(params_dict, "if_save_data", true)
	save_every::Int = get(params_dict, "save_every", 2)
	nrgtol = get(params_dict, "nrgtol", nothing)

	# inputs (ts, filling, lr, if_check_fluxes, scaling_paras) are kept next to what is derived from them
	model_paras::Dict{String,Any} = Dict([("Lx",Lx),("Ly",Ly),("particles",particles),("filling",filling),("num_sites",num_sites),("if_periodic_x",if_periodic_x),("if_periodic_y",if_periodic_y),
		("ts",t_strength),("t_x",t_x),("t_y",t_y),("hopping_anisotropy",hopping_anisotropy),("twist_angle",twist_angle),("alpha",alpha),("flux_direction",flux_direction),("if_check_fluxes",if_check_fluxes),
		("onsite_strength",onsite_strength),("scaling",scaling),("which_dir",which_dir),("lr",get(params_dict, "lr", 0)),
		("interaction_axis_length",interaction_axis_length),("long_range_strengths",long_range_strengths),
		("chem_strength",chem_strength),("if_pinning",if_pinning),("pinning_strength",pinning_strength),("outputlevel",outputlevel),
		("dataloc",dataloc),("if_save_data",if_save_data),("save_every",save_every),("nrgtol",nrgtol),("bond_dim",bond_dim)])

	return merge(model_paras, scaling_paras)
end

function long_range_HH_ham(model_paras::Dict{String,Any})

	Lx,Ly = model_paras["Lx"],model_paras["Ly"]
	L = (Lx,Ly)
	num_sites = Lx*Ly
	outputlevel = model_paras["outputlevel"]

	# lattice
	if_periodic_x = false
	if_periodic_y = false
	if_per = (if_periodic_x,if_periodic_y)

	# hopping and flux
	t_x,t_y = model_paras["t_x"],model_paras["t_y"]
	twist_angle = model_paras["twist_angle"]
	alpha = Float64(model_paras["alpha"])
	flux_direction = model_paras["flux_direction"]

	# interactions, which_dir is the axis the long-range interaction acts along ("x", "y" or "both")
	long_range_strengths = model_paras["long_range_strengths"]
	which_dir = model_paras["which_dir"]

	# on-site potentials
	chem_strength = get(model_paras, "chem_strength", 0.0)
	if_pinning = get(model_paras, "if_pinning", false)
	pinning_strength = get(model_paras, "pinning_strength", 1e-4)

	ham = ITensors.OpSum()

	for i in 1:num_sites
		starting_site = coordinate(i,Lx,Ly)
		for axis in [1,2]
			ending_site = collect(starting_site)
			ending_site[axis] += 1
			twist = 1.0

			# enforce boundary conditions
			if ending_site[axis] > L[axis]
				if_per[axis] ? nothing : continue
				ending_site[axis] = 1
				twist = exp(im*twist_angle[axis]*2*pi)
			end

			coeff = get_hopping_coeff(starting_site,axis,t_x,t_y,alpha,flux_direction) * twist
			ham += (coeff,"S+",starting_site,"S-",Tuple(ending_site))
			ham += (conj(coeff),"S+",Tuple(ending_site),"S-",starting_site)
		end
	end

	# long_range_strengths[d+1] couples each pair a distance d apart along the interaction axis, open boundaries,
	# the on-site d = 0 term drops out for hardcore bosons
	for i in 1:num_sites, j in i+1:num_sites
		s1,s2 = coordinate(i,Lx,Ly),coordinate(j,Lx,Ly)
		_,(dx,dy),_ = find_dist(s1,s2,(Lx,Ly))
		on_axis = which_dir == "x" ? dy == 0 : which_dir == "y" ? dx == 0 : dx == 0 || dy == 0
		on_axis || continue
		stren = long_range_strengths[Int(dx + dy) + 1]
		stren == 0.0 && continue
		ham += (stren,"ProjUp",s1,"ProjUp",s2)
	end

	if chem_strength != 0.0
		for i in 1:num_sites
			ham -= (chem_strength,"ProjUp",coordinate(i,Lx,Ly))
		end
	end

	if if_pinning
		ham += (pinning_strength,"ProjUp",(1,1))
	end

	return ham
end

# Initial PEPS with at most one particle per column (fixed x, the interaction axis is y): each column is the
# y-MPS of a|empty> + b|W> (W = one particle spread evenly along y), so the y-links carry a counter
# (1 = no particle yet, 2 = placed) and the x-links only their first state, the rest of bond_dim is zeros.
# occupied_columns = nothing has a = b = 1 on every column, i.e. in the N sector all configurations with
# <= 1 particle per column at equal weight; a list has a = 0, b = 1 on those columns and a = 1, b = 0 on
# the others (1:2:Lx is Tao-Thouless-like). noise*randn is added to every entry so no amplitude is exactly zero.
# kwargs go to PEPS (contract_cutoff, sample_cutoff, ...)
function column_exclusion_peps(::Type{T}, hilbert; bond_dim::Int=5, occupied_columns=nothing, noise=1e-2, kwargs...) where {T<:Number}
	Lx, Ly = size(hilbert)
	bond_dim >= 2 || error("column_exclusion_peps needs bond_dim >= 2 for the column counter")
	peps = PEPS(T, hilbert; bond_dim, kwargs...)
	up, dn = 1, 2
	for i in 1:Lx, j in 1:Ly
		old_inds = inds(peps[i,j])
		s = hilbert[i,j]
		yl = j > 1 ? commonind(peps[i,j], peps[i,j-1]) : nothing
		yr = j < Ly ? commonind(peps[i,j], peps[i,j+1]) : nothing
		xls = [commonind(peps[i,j], peps[k,j]) for k in (i-1,i+1) if 1 <= k <= Lx]
		a, b = isnothing(occupied_columns) ? (1.0, 1.0) : (i in occupied_columns ? (0.0, 1.0) : (1.0, 0.0))

		t = ITensor(noise .* randn(T, dims(old_inds)...), old_inds...)
		# counter transitions (left -> right, occupation); the last site of a column weighs the end states with (a, b)
		for (cl, cr, occ) in ((1,1,dn), (1,2,up), (2,2,dn))
			(isnothing(yl) && cl != 1) && continue
			amp = isnothing(yr) ? (cr == 1 ? a : b) : 1.0
			idx = [s => occ; isnothing(yl) ? [] : [yl => cl]; isnothing(yr) ? [] : [yr => cr]; [xl => 1 for xl in xls]]
			t[idx...] += amp
		end
		peps[i,j] = t
	end
	return peps
end

include("../synth-dims/long-range-ttn.jl")
include("../exact-diag/two-dimensions.jl")
using Test,PyPlot

function rebuild_ed_ham(ttn_ham,lattice_params::Dict)
	full_basis = lattice_params["full_basis"]
	rebuilt_ham = spzeros(ComplexF64,size(full_basis)[2],size(full_basis)[2])
	Lx = lattice_params["Lx"]
	make_smaller_lattice = [Lx,Lx]

	interacting_pairs = []
	hopping_pairs = []
	for t1 in TTNKit.terms(ttn_ham)
		coeff = TTNKit.coefficient(t1)
		if length(t1) == 1
			continue
		end
					
		both_terms = TTNKit.terms(t1)

		which_sites = TTNKit.site.(both_terms)
		if any(which_sites[1] .> Lx) || any(which_sites[2] .> Lx)
			continue
		end
		ops_here = TTNKit.which_op.(both_terms)
		if ops_here == ["Adag","A"]
			end_site,start_site = linear_index.(which_sites,make_smaller_lattice[1],make_smaller_lattice[2])
			#println("Moving from $start_site to $end_site with coeff $coeff")
			append!(hopping_pairs,[(start_site,end_site,coeff)])
		elseif ops_here == ["Adag * A","Adag * A"]
			s1,s2 = linear_index.(which_sites,make_smaller_lattice[1],make_smaller_lattice[2])
			#println("Interacting at sites $s1 and $s2 with coeff $coeff")
			append!(interacting_pairs,[(s1,s2,coeff)])
		end
	end

	for j in 1:size(full_basis)[2]
		#println(round(100*j/size(full_basis)[2],digits=2))
		basis = full_basis[:,j]
					
		# check for hopping
		for (start_site,end_site,coeff) in hopping_pairs
			if start_site in basis && !(end_site in basis)
				where_start_site = findfirst(x -> x == start_site,basis)
				new_basis = copy(basis)
				new_basis[where_start_site] = end_site
				sort!(new_basis,rev=true)
				new_basis_index = find_basis_index(new_basis)
				rebuilt_ham[new_basis_index,j] += coeff
			end
		end

		# add interaction terms
		for (s1,s2,coeff) in interacting_pairs
			if s1 in basis && s2 in basis
				rebuilt_ham[j,j] += coeff
			end
		end
	end

	return rebuilt_ham
end

open_cores = 5
if typeof(open_cores) != String
    BLAS.set_num_threads(open_cores)	
    display(BLAS.get_config())
end

cols = ["b","r","g","m","c"]
es_count = 3
for phys in [true,false]
	for virt in [true,false]
		if phys && virt
			what_type = "Torus"
		elseif phys && !virt
			what_type = "Cylinder Phys"
		elseif !phys && virt
			what_type = "Cylinder Virt"
		else
			what_type = "OBC"
		end
		fig = figure()
		xlabel("Filling")
		ylabel("Energy")
		title("ED vs TTN Energy Levels $what_type 4x4 N=2")
		all_ed_nrgs = zeros(es_count+1,37)
		all_ttn_nrgs = zeros(es_count+1,10)
		ttn_count = 0
		for (idx,nu) in enumerate(range(0.0,2.0,length=37))

			params_dict = Dict([("if_save_data",true),("lattice_size",(4,4)),("mdim",100),("filling",nu),("if_periodic_phys",phys),("if_periodic_virt",virt)])

			nrgtol = 5E-5
			onsite_strength = get(params_dict, "onsite_strength", 0.0)
			anis = get(params_dict, "hopping_anisotropy", 1.0)
			mdim = get(params_dict, "mdim", 10)
			longrange_dist = get(params_dict, "lr", 0)
			alpha = get(params_dict, "alpha", nothing)
			phys_edge_length,virt_edge_length = get(params_dict, "lattice_size", [4,4])
			if 2 <= phys_edge_length < 4
				layer_count = 2
			elseif 4 <= phys_edge_length < 8
				layer_count = 4
			else
				layer_count = 6
			end
			num_particles = 2#get(params_dict, "particles", Int(phys_edge_length/2))

			if_per_phys = get(params_dict, "if_periodic_phys", true)
			if_per_virt = get(params_dict, "if_periodic_virt", false)
			if isnothing(alpha)
				filling = get(params_dict, "filling", 1.0)
				if filling == 0.0
					alpha = 0.0
					mag_off = true
				else
					phys_shift,synth_shift = !if_per_phys,!if_per_virt
					alpha = num_particles/(filling*(phys_edge_length - phys_shift)*(virt_edge_length - synth_shift))
					mag_off = false
				end
			else
				mag_off = alpha == 0.0
			end
			#check_fluxes(alpha,phys_edge_length,virt_edge_length,if_per_phys,if_per_virt)
				
			expan = TTNKit.DefaultExpander(1.0)#TTNKit.NoExpander()
			ts = 1.0
			tot_sites = 2^layer_count
			syms = get(params_dict, "syms", true)

			nswps = 10

			save_data = get(params_dict, "if_save_data", true)
			if_cluster = any([occursin("local",pwd()),occursin("Local",pwd()),occursin("geraghty",pwd())])
			if_continuous_saving = get(params_dict,"if_continuous_saving",if_cluster || layer_count >= 7)
			save_data ? nothing : if_continuous_saving = false

				
			if_cliff = false
			trunc = get(params_dict,"trunc",1e-3)
			sc_type = "flat"


			loc = get_folder_location("cluster-data/synth-dims/excited-states/testing")
				
			filename_dict = Dict([("layers",layer_count),("lr",longrange_dist),("particles",num_particles),("alpha",round(alpha,digits=4)),("if_periodic_phys",if_per_phys),("if_periodic_synth",if_per_virt)])
			if [phys_edge_length,virt_edge_length] != [sqrt(2^layer_count),sqrt(2^layer_count)]
				filename_dict["make_smaller_lattice"] = phys_edge_length
			end

			datafile_name = make_parameters_filename(filename_dict)
					#end
			model_paras = (hopping_anisotropy=anis,
									syms=syms,
									cutoff=0.0,
									twist_angle=0.0,
									if_continuous_saving=if_continuous_saving,
									nrgtol=nrgtol,
									if_densmat=true,
									restricted_size=(phys_edge_length,virt_edge_length),
									centralflux_strength=0.0,
									if_pinning_pot=false,
									if_periodic_phys=if_per_phys,
									if_periodic_virt=if_per_virt,
									if_nn_int=false,
									nn_int_strength=0.0,
									chem_strength=0.0,
									no_magF=mag_off,
									scaling=sc_type,
									scaling_dist=longrange_dist,
									onsite_strength=onsite_strength,
									which_dir="virt",
									cliff=false,
									trunc=1e-8,
									if_change=false,
									change=0.0,
									if_gpu=false,
									noise=[0.0],
									if_save_data=save_data,
									if_save_fig=false,
									if_sweep=true,
									sweep_type="dmrg",
									expander=expan,
									max_occ=1,
									mdim=mdim,
									weight=100.0,
									num_sweeps=nswps,
									phi=alpha,
									output_level=0,
									name="ttn-"*datafile_name,
									location=loc)
					
			metadata_dict = merge(named_tuple_to_dict(model_paras),filename_dict)

			file_list = find_data_file(filename_dict,"ttn",loc; output_level=0)
			if length(file_list) != 0
				ttn_count += 1
				data,metadata = read_data_jld2(file_list[1],loc; output_level=0)
				ham = metadata["ham"]
				ttn_nrgs = [metadata["observer"].nrg[end]]
				for i in 1:es_count
					append!(ttn_nrgs,metadata["observer_$i"].nrg[end])
				end
				all_ttn_nrgs[:,ttn_count] = ttn_nrgs
			else
				net = build_HH_net(layer_count; syms=syms, max_occ=1)
				ham = long_range_HH_ham(net,ts,alpha; model_paras...)
			end
			
			#= finding ground state TTN
			if true
				println("Starting Script using $num_particles particles on $tot_sites sites with $(!mag_off) Mag Field, Bond Dim = $mdim, and Long Range Dist = $longrange_dist")
				starting = time()
				net = build_HH_net(layer_count; syms=syms, max_occ=1)
				ham = long_range_HH_ham(net,ts,alpha; model_paras...)
				og_ttn, hamilt, dm_sp, rezobs, runtime, dens = find_ground_state(layer_count,num_particles; ttn_net=net,ham_op=ham,model_paras...,metadata=merge(metadata_dict,Dict([("ham",ham),("net",net),("t_strength",1.0)])))
				total_time = time() - starting
				println("Running time = $total_time")
			end

			# finding excited states TTN
			if true
				all_states, hamilt, all_obs, all_densmats, all_runtimes = find_excited_states(layer_count,es_count,num_particles,dm_sp.ttn; ham_op=ham,model_paras...,metadata=merge(metadata_dict,Dict([("ham",ham),("t_strength",ts)])))
			end=#

			# finding all states ED
			if true
				lattice_params = Dict([("Lx",phys_edge_length),("Ly",virt_edge_length),("N",num_particles),("full_basis",n_particle_basis(num_particles,phys_edge_length,phys_edge_length; output_level=1)),("if_periodic_x",if_per_phys),("if_periodic_y",if_per_virt)])
				ed_ham = rebuild_ed_ham(ham,lattice_params)
				x0_ed = rand(Float64,size(lattice_params["full_basis"])[2])
				rez_ed = eigsolve(ed_ham,x0_ed,es_count+1,:SR,Lanczos())
				sorted_indices_ed = sortperm(rez_ed[1])
				states_ed = rez_ed[2][sorted_indices_ed][1:es_count+1]
				nrgs_ed = rez_ed[1][sorted_indices_ed][1:es_count+1]
			end

			all_ed_nrgs[:,idx] = nrgs_ed

			#println("ED Energies:")
			#display(nrgs_ed)

			#ttn_nrgs = [[rezobs.nrg[end]]; [ob_here.nrg[end] for ob_here in all_obs]]
			#println("TTN Energies:")
			#display(ttn_nrgs)

			#println("NRG Diffs:")
			#nrg_diffs = abs.(nrgs_ed - ttn_nrgs)
			#display(nrg_diffs)
			#for i in 1:es_count+1
			#	scatter(nu,nrg_diffs[i],c=cols[i],label="$(i-1)")
			#end
		end

		for i in 1:es_count+1
			if i == 1
				plot(collect(range(0.0,2.0,length=37)),all_ed_nrgs[i,:],c=cols[i],label="ED")
				scatter(collect(range(0.0,2.0,length=10)),all_ttn_nrgs[i,:],c=cols[i],label="TTN")
			else
				plot(collect(range(0.0,2.0,length=37)),all_ed_nrgs[i,:],c=cols[i])
				scatter(collect(range(0.0,2.0,length=10)),all_ttn_nrgs[i,:],c=cols[i])
			end
		end
		legend()
	end
end










































"fin"
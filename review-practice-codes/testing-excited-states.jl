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

#=open_cores = 5
if typeof(open_cores) != String
    BLAS.set_num_threads(open_cores)	
    display(BLAS.get_config())
end=#


if false
	cols = ["b","r","g","m","c"]
	es_count = 10
	#for phys in [true,false]
		#for virt in [true,false]
			phys = true
			virt = true
			if phys && virt
				what_type = "Torus"
			elseif phys && !virt
				what_type = "Cylinder Phys"
			elseif !phys && virt
				what_type = "Cylinder Virt"
			else
				what_type = "OBC"
			end
			#fig = figure()
			#xlabel("Filling")
			#ylabel("Energy")
			#title("ED vs TTN Energy Levels $what_type 4x4 N=2")
			all_ed_nrgs = zeros(es_count+1,37)
			all_ttn_nrgs = zeros(es_count+1,10)
			ttn_count = 0
			#for (idx,nu) in enumerate(range(0.0,2.0,length=37))
			nu = 0.5

				params_dict_ttn = Dict([("es_count",es_count),("if_find_data",false),("if_save_data",false),("layers",4),("mdim",50),("filling",nu),("if_periodic_phys",phys),("if_periodic_synth",virt)])
				params_dict_ed = Dict([("nev",es_count+1),("if_find_data",false),("if_save_data",false),("Lx",4),("N",2),("filling",nu),("if_periodic_x",phys),("if_periodic_y",virt)])

				if true
					all_states, hamilt, all_obs, all_densmats, all_runtimes = run_synth_dims_generic(params_dict_ttn)
				end

				# finding all states ED
				if true
					lattice_params,hamilt_params,running_args = get_normal_model_params_ed(params_dict_ed)
					lattice_params["full_basis"] = n_particle_basis(lattice_params)
					states_ed,nrgs_ed,rhos_ed,ed_ham = find_eigenstates(running_args.nev,lattice_params,hamilt_params; running_args...)

					#=model_paras_ttn = get_normal_model_params(params_dict_ttn)
					phys_edge_length,virt_edge_length = model_paras[:restricted_size]
					if_periodic_phys,if_periodic_synth = model_paras[:if_periodic_phys],model_paras[:if_periodic_synth]
					num_particles = model_paras[:particles]
					net = build_HH_net(model_paras)
					model_paras[:net] = net
					ham = long_range_HH_ham(named_tuple_to_dict(model_paras))
					ttn_ham_rebuilt = rebuild_ed_ham(ham,lattice_params)

					x0_ed_ttnham = rand(Float64,size(lattice_params["full_basis"])[2])
					rez_ed_ttnham = eigsolve(ttn_ham_rebuilt,x0_ed,es_count+1,:SR,Lanczos())
					sorted_indices_ed_ttnham = sortperm(rez_ed[1])
					states_ed_ttnham = rez_ed[2][sorted_indices_ed][1:es_count+1]
					nrgs_ed_ttnham = rez_ed[1][sorted_indices_ed][1:es_count+1]=#


				end

				#all_ed_nrgs[:,idx] = nrgs_ed

				#println("ED from TTN-Ham Energies:")
				#display(nrgs_ed_ttnham)

				println("ED Energies:")
				display(nrgs_ed)

				ttn_nrgs = [obs.nrg[end] for obs in all_obs]
				println("TTN Energies:")
				display(ttn_nrgs)

				#println("NRG Diffs:")
				#nrg_diffs = abs.(nrgs_ed - ttn_nrgs)
				#display(nrg_diffs)
				#for i in 1:es_count+1
				#	scatter(nu,nrg_diffs[i],c=cols[i],label="$(i-1)")
				#end
			#end

			#=for i in 1:es_count+1
				if i == 1
					plot(collect(range(0.0,2.0,length=37)),all_ed_nrgs[i,:],c=cols[i],label="ED")
					scatter(collect(range(0.0,2.0,length=10)),all_ttn_nrgs[i,:],c=cols[i],label="TTN")
				else
					plot(collect(range(0.0,2.0,length=37)),all_ed_nrgs[i,:],c=cols[i])
					scatter(collect(range(0.0,2.0,length=10)),all_ttn_nrgs[i,:],c=cols[i])
				end
			end
			legend()=#
		#end
	#end

end


if true
	params_dict = Dict([("es_count",0),("particles",2),("if_save_data",false),("layers",4),("mdim",10),("filling",0.5),("if_periodic_phys",true),("if_periodic_synth",true)])

	#og_ttn, hamilt, gs_sp, gs_obs, gs_runtime, gs_dens = run_synth_dims_generic(params_dict)
	#wavefunc = gs_sp.ttn

	# need a two different wavefunctions that have a non-zero overlap
	model_paras = get_normal_model_params(params_dict)
	net1 = build_HH_net(model_paras)
	net2 = build_HH_net(model_paras)

	psi1 = TTNKit.RandomTreeTensorNetwork(net1)

	states = fill("0", 2^(model_paras[:layers]))
	old_ttn = TTNKit.ProductTreeTensorNetwork(net2,states)
	psi2 = initialize_ttn(old_ttn,model_paras[:mdim],model_paras[:particles]; model_paras...)

	if true
		@testset "Up flow for overlap" begin
			upflow = TTNKit.bottom_overlap_environments(psi1,psi2)

			tensorlist = vcat(psi1[(params_dict["layers"],1)],(upflow[params_dict["layers"]][1])...,dag(psi2[(params_dict["layers"],1)]))
			overlap_fromflow = TTNKit.ITensors.scalar(TTNKit.contract(tensorlist))
			overlap_direct = TTNKit.inner(psi1,psi2)
			@test isapprox(overlap_fromflow,overlap_direct,atol=1e-10)
		end
	end

	if true
		@testset "Overlap environments" begin

			overlap_direct = TTNKit.inner(psi1,psi2)

			full_envs = TTNKit.build_overlap_environments(psi1,psi2)
			for ll in 1:params_dict["layers"]-1
				for nn in 1:2^(params_dict["layers"]-ll)
					top,bot_left,bot_right = full_envs[ll][nn]
					overlap_fromenvs = TTNKit.ITensors.scalar(TTNKit.contract(top,psi1[(ll,nn)],bot_left,bot_right,dag(psi2[(ll,nn)])))
					@test isapprox(overlap_fromenvs,overlap_direct,atol=1e-10)
				end
			end

		end
	end

end








































"fin"
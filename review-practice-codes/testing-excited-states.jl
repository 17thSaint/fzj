include("../synth-dims/long-range-ttn.jl")
include("../exact-diag/two-dimensions.jl")
include("extra-testing-functions.jl")
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


if true
	cols = ["b","r","g","m","c","k"]
	es_count = 2
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
			
			fillings_ttn = range(0.5,3.0,length=10)
			fillings_ed = range(0.5,3.0,length=37)

			for (idx,nu) in enumerate(fillings_ed)

				params_dict_ed = Dict([("nev",es_count+1),("if_check_fluxes",false),("if_find_data",false),("if_save_data",false),("Lx",4),("N",2),("filling",nu),("if_periodic_x",phys),("if_periodic_y",virt)])

				lattice_params,hamilt_params,running_args = get_normal_model_params_ed(params_dict_ed)
				lattice_params["full_basis"] = n_particle_basis(lattice_params)
				states_ed,nrgs_ed,rhos_ed,ed_ham = find_eigenstates(running_args.nev,lattice_params,hamilt_params; running_args...)

				all_ed_nrgs[:,idx] = nrgs_ed
			end

			for i in 1:es_count+1
				if i == 1
					plot(collect(fillings_ed),all_ed_nrgs[i,:],c=cols[i],label="ED")
				else
					plot(collect(fillings_ed),all_ed_nrgs[i,:],c=cols[i])
				end
			end

			for (idx,nu) in enumerate(fillings_ttn)

				params_dict_ttn = Dict([("if_old_excited",false),("if_check_fluxes",false),("es_count",es_count),("if_find_data",false),("if_save_data",false),("layers",4),("mdim",50),("filling",nu),("if_periodic_phys",phys),("if_periodic_synth",virt)])

				all_states, hamilt, all_obs, all_densmats, all_runtimes = run_synth_dims_generic(params_dict_ttn)

				all_ttn_nrgs[:,idx] = [obs.nrg[end] for obs in all_obs]

				for i in 1:es_count+1
					scatter(nu,all_ttn_nrgs[i,idx],c=cols[i])
				end
				legend()
			end

				#println("ED from TTN-Ham Energies:")
				#display(nrgs_ed_ttnham)

				#println("ED Energies:")
				#display(nrgs_ed)

				#ttn_nrgs = [obs.nrg[end] for obs in all_obs]
				#all_ttn_nrgs[:,idx] = ttn_nrgs
				#println("TTN Energies:")
				#display(ttn_nrgs)

				#println("NRG Diffs:")
				#nrg_diffs = abs.(nrgs_ed - ttn_nrgs)
				#display(nrg_diffs)
				#for i in 1:es_count+1
				#	scatter(nu,nrg_diffs[i],c=cols[i],label="$(i-1)")
				#end

		#end
	#end

end


if false

	if_local_all = false

	params_dict = Dict([("particles",2),("layers",4),("mdim",10)])
	
	max_occ = 1
	net1 = simple_boson_network(params_dict["layers"],true,max_occ)
	lat1 = TTNKit.physical_lattice(net1)
	psi1 = make_randomconfig_ttn(net1,params_dict["particles"],max_occ)

	#net2 = simple_boson_network(params_dict["layers"],true,max_occ)
	psi2 = make_parton_ttn(net1,params_dict["particles"],params_dict["mdim"],max_occ)

	psi3 = make_parton_ttn(net1,params_dict["particles"],params_dict["mdim"],max_occ)

	if false || if_local_all
		@testset "Up flow for overlap" begin
			upflow = TTNKit.bottom_overlap_environments(psi1,psi2)

			tensorlist = vcat(psi1[(params_dict["layers"],1)],(upflow[params_dict["layers"]][1])...,prime(dag(psi2[(params_dict["layers"],1)])))
			overlap_fromflow = TTNKit.ITensors.scalar(TTNKit.contract(tensorlist))
			overlap_direct = TTNKit.inner(psi1,psi2)
			@test isapprox(overlap_fromflow,overlap_direct,atol=1e-10)
		end
	end

	if false || if_local_all
		@testset "Overlap environments with existing inner product" begin

			overlap_direct = TTNKit.inner(psi1,psi2)

			full_envs = TTNKit.build_overlap_environments(psi1,psi2)
			for ll in 1:params_dict["layers"]-1
				for nn in 1:2^(params_dict["layers"]-ll)
					overlap_fromenvs = TTNKit.ITensors.scalar(TTNKit.contract(psi1[(ll,nn)],full_envs[ll][nn]))
					@test isapprox(overlap_fromenvs,overlap_direct,atol=1e-10)
				end
			end

		end
	end

	if false || if_local_all
		@testset "Overlap environments equal slow inner function" begin
			for ll in 1:params_dict["layers"]-1
				for nn in 1:2^(params_dict["layers"]-ll)
					which_site = (ll,nn)
					oldinner_overlap = dag(noprime(TTNKit.inner(psi1,psi2,which_site)))

					full_envs = TTNKit.build_overlap_environments(psi1,psi2)
					flow_env = full_envs[which_site[1]][which_site[2]]

					@test isapprox(oldinner_overlap[1],flow_env[1],atol=1e-10)
				end
			end
		end
	end

	if false || if_local_all
		@testset "Overlap at top node" begin
			overlap_direct = TTNKit.inner(psi1,psi2)

			full_envs = TTNKit.build_overlap_environments(psi1,psi2)

			tensor_list = [psi1[(params_dict["layers"],1)],full_envs[params_dict["layers"]][1],prime(dag(psi2[(params_dict["layers"],1)]))]
			overlap_fromenvs = TTNKit.ITensors.scalar(TTNKit.contract(tensor_list))

			@test isapprox(overlap_direct,overlap_fromenvs,atol=1e-10)
		end
	end

	if false
		ham = long_range_HH_ham(net1,1.0,0.0; scaling="flat")
		ham_tpo = TTNKit.TPO(ham,lat1)

		rez = TTNKit.dmrg(psi1, [psi2], ham_tpo; maxdims=params_dict["mdim"],number_of_sweeps=1,if_old_excited=false)
	end

end








































"fin"
#using Pkg
#Pkg.activate(".")
using Test

function find_center()
	all_folders = split(pwd(),"/")
	if "fzj" in all_folders
		return "fzj"
	elseif "local" in all_folders
		return all_folders[findfirst(x -> all_folders[x] == "local",1:length(all_folders))+1]
	elseif "Local" in all_folders
		return all_folders[findfirst(x -> all_folders[x] == "Local",1:length(all_folders))+1]
	else
		println("Not sure where the center is: $(pwd())")
	end
end

function include_other_files(all_files,output_level=0)
	center = find_center()
	get_to_fzj = split(pwd(),center)[1]
	if typeof(all_files) == String
		all_files = [all_files]
	end
	for file in all_files
		occursin("main-git",pwd()) ? include(get_to_fzj * center * "/main-git/" * file) : include(get_to_fzj * center * "/" * file)
		output_level > 0 ? println("Included $file") : nothing
	end
end

include_other_files(["other-funcs/data-storage-funcs.jl","synth-dims/long-range-ttn.jl","exact-diag/two-dimensions.jl","synth-dims/old-long-range-ttn.jl"])

if_all = false
if_plot = false

if false || if_all
    @testset "Torus 4x4 N=4" begin
        ttn_dataloc = get_folder_location("cluster-data/synth-dims/torus")
        ed_dataloc = get_folder_location("cluster-data/exact-diag/torus")

        Lx,Ly = 4,4
        N = 4

        ttn_files = find_data_file(Dict([("layers",Int(log(2,Lx*Ly))),("particles",N),("if_periodic_phys",true),("if_periodic_virt",true),("hopping_anisotropy",1.0),("onsite_strength",0.0)]),"ttn",ttn_dataloc; output_level=0)
        ed_files = find_data_file(Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",0.0)]),"ed",ed_dataloc; output_level=0)

        ttn_data = Dict()
        for f in ttn_files
            data,metadata = read_data_jld2(f,ttn_dataloc; output_level=0)
			ttn_data[metadata["alpha"] == 0.0 ? "sf" : "fqh"] = (data["ttn"],data["densmat"],metadata["observer"].nrg)
		end

		ed_data = Dict()
		for f in ed_files
			data,metadata = read_data_jld2(f,ed_dataloc; output_level=0)
			ed_data[metadata["alpha"] == 0.0 ? "sf" : "fqh"] = (data["state"],data["densmat"],data["nrg"])
		end

		# test that the final energy of the TTN is with 0.1% of the ED
		@test abs((ttn_data["sf"][3][end] - ed_data["sf"][3][1]) / ed_data["sf"][3][1]) < 1e-3
		@test abs((ttn_data["fqh"][3][end] - ed_data["fqh"][3][1]) / ed_data["fqh"][3][1]) < 1e-3

		# test that the density matrix is the same for the Superfluid state
		densmat_diff_sf = ttn_data["sf"][2] .- ed_data["sf"][2][1]
		if if_plot
			fig = figure()
			imshow(real.(densmat_diff_sf))
			colorbar()
			title("TTN - ED Superfluid")
		end
		@test maximum(abs.(densmat_diff_sf)) < 1e-8

		ed_fqh_densmat = ed_data["fqh"][2][1]
		reordered_densmat = zeros(ComplexF64,size(ed_fqh_densmat))
		for i in 1:Lx*Ly
			for j in 1:Lx*Ly
				s1_coord = coordinate(i,Lx,Ly)
				s2_coord = coordinate(j,Lx,Ly)
				s1_old_lin = s1_coord[2] + (s1_coord[1]-1)*Lx
				s2_old_lin = s2_coord[2] + (s2_coord[1]-1)*Lx
				reordered_densmat[i,j] = ed_fqh_densmat[s2_old_lin,s1_old_lin]
			end
		end

		# test that the density matrix is the same for the FQH state but include reordering from TTN density matrix function
		densmat_diff_fqh = ttn_data["fqh"][2] .- reordered_densmat#ed_data["fqh"][2][1]
		if if_plot
			fig = figure()
			imshow(imag.(densmat_diff_fqh))
			#plot(real.(densmat_diff_fqh))
			colorbar()
			title("TTN - ED FQH")
		end
		@test maximum(abs.(densmat_diff_fqh)) < 1e-7#

		# test greens function for the Superfluid state
		ttn_phys_corr_sf = physical_correlation(ttn_data["sf"][2],Lx,Ly; plot_title="TTN Superfluid", if_plot=false)		
		ed_phys_corr_sf = physical_correlation(ed_data["sf"][2][1],Lx,Ly; plot_title="ED Superfluid", if_plot=false)
		if if_plot
			plot_physical_correlation(abs.(ttn_phys_corr_sf .- ed_phys_corr_sf); plot_title="TTN - ED Superfluid")
		end
		@test isapprox(ttn_phys_corr_sf,ed_phys_corr_sf,atol=1e-5)

		ttn_syn_corr_sf = synthetic_correlation(ttn_data["sf"][2],Lx,Ly; plot_title="TTN Superfluid", if_plot=false)
		ed_syn_corr_sf = synthetic_correlation(ed_data["sf"][2][1],Lx,Ly; plot_title="ED Superfluid", if_plot=false)
		if if_plot
			plot_synthetic_correlation(abs.(ttn_syn_corr_sf .- ed_syn_corr_sf); plot_title="TTN - ED Superfluid")
		end
		@test isapprox(ttn_syn_corr_sf,ed_syn_corr_sf,atol=1e-5)

		# test greens function for the FQH state
		ttn_phys_corr_fqh = physical_correlation(ttn_data["fqh"][2],Lx,Ly; plot_title="TTN FQH", if_plot=false)
		ed_phys_corr_fqh = physical_correlation(ed_data["fqh"][2][1],Lx,Ly; plot_title="ED FQH", if_plot=false)
		if if_plot
			plot_physical_correlation(abs.(ttn_phys_corr_fqh .- ed_phys_corr_fqh); plot_title="TTN - ED FQH")
		end
		@test isapprox(ttn_phys_corr_fqh,ed_phys_corr_fqh,atol=1e-5)

		ttn_syn_corr_fqh = synthetic_correlation(ttn_data["fqh"][2],Lx,Ly; plot_title="TTN FQH", if_plot=false)
		ed_syn_corr_fqh = synthetic_correlation(ed_data["fqh"][2][1],Lx,Ly; plot_title="ED FQH", if_plot=false)
		if if_plot
			plot_synthetic_correlation(abs.(ttn_syn_corr_fqh .- ed_syn_corr_fqh); plot_title="TTN - ED FQH")
		end
		@test isapprox(ttn_syn_corr_fqh,ed_syn_corr_fqh,atol=1e-5)#

    end
end

if false || if_all
	@testset "OBC 4x4 N=4" begin
		ttn_dataloc = get_folder_location("cluster-data/synth-dims/obc")
		ed_dataloc = get_folder_location("cluster-data/exact-diag/obc")

		Lx,Ly = 4,4
		N = 4

		ttn_files = find_data_file(Dict([("layers",Int(log(2,Lx*Ly))),("particles",N),("if_periodic_phys",false),("if_periodic_virt",false),("hopping_anisotropy",1.0),("onsite_strength",0.0)]),"ttn",ttn_dataloc; output_level=0)
		ed_files = find_data_file(Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",false),("if_periodic_y",false),("hopping_anisotropy",1.0),("interaction_strength",0.0)]),"ed",ed_dataloc; output_level=0)

		ttn_data = Dict()
		for f in ttn_files
			data,metadata = read_data_jld2(f,ttn_dataloc; output_level=0)
			ttn_data[metadata["alpha"] == 0.0 ? "sf" : "fqh"] = (data["ttn"],data["densmat"],metadata["observer"].nrg)
		end

		ed_data = Dict()
		for f in ed_files
			data,metadata = read_data_jld2(f,ed_dataloc; output_level=0)
			ed_data[metadata["alpha"] == 0.0 ? "sf" : "fqh"] = (data["state"],data["densmat"],data["nrg"])
		end

		# test that the final energy of the TTN is with 0.1% of the ED
		@test abs((ttn_data["sf"][3][end] - ed_data["sf"][3][1]) / ed_data["sf"][3][1]) < 1e-3
		@test abs((ttn_data["fqh"][3][end] - ed_data["fqh"][3][1]) / ed_data["fqh"][3][1]) < 1e-3

		# test that the density matrix is the same for the Superfluid state
		densmat_diff_sf = ttn_data["sf"][2] .- ed_data["sf"][2][1]
		if if_plot
			fig = figure()
			imshow(real.(densmat_diff_sf))
			colorbar()
			title("TTN - ED Superfluid")
		end
		@test maximum(abs.(densmat_diff_sf)) < 1e-8#

		ed_fqh_densmat = ed_data["fqh"][2][1]
		reordered_densmat = zeros(ComplexF64,size(ed_fqh_densmat))
		for i in 1:Lx*Ly
			for j in 1:Lx*Ly
				s1_coord = coordinate(i,Lx,Ly)
				s2_coord = coordinate(j,Lx,Ly)
				s1_old_lin = s1_coord[2] + (s1_coord[1]-1)*Lx
				s2_old_lin = s2_coord[2] + (s2_coord[1]-1)*Lx
				reordered_densmat[i,j] = ed_fqh_densmat[s2_old_lin,s1_old_lin]
			end
		end

		# test that the density matrix is the same for the FQH state but include reordering from TTN density matrix function
		densmat_diff_fqh = ttn_data["fqh"][2] .- reordered_densmat#ed_data["fqh"][2][1]
		if if_plot
			fig = figure()
			imshow(imag.(densmat_diff_fqh))
			#plot(real.(densmat_diff_fqh))
			colorbar()
			title("TTN - ED FQH")
		end
		@test maximum(abs.(densmat_diff_fqh)) < 1e-8#

		# test greens function for the Superfluid state
		ttn_phys_corr_sf = physical_correlation(ttn_data["sf"][2],Lx,Ly; plot_title="TTN Superfluid", if_plot=false)
		ed_phys_corr_sf = physical_correlation(ed_data["sf"][2][1],Lx,Ly; plot_title="ED Superfluid", if_plot=false)
		if if_plot
			plot_physical_correlation(abs.(ttn_phys_corr_sf .- ed_phys_corr_sf); plot_title="TTN - ED Superfluid")
		end
		@test isapprox(ttn_phys_corr_sf,ed_phys_corr_sf,atol=1e-5)

		ttn_syn_corr_sf = synthetic_correlation(ttn_data["sf"][2],Lx,Ly; plot_title="TTN Superfluid", if_plot=false)
		ed_syn_corr_sf = synthetic_correlation(ed_data["sf"][2][1],Lx,Ly; plot_title="ED Superfluid", if_plot=false)
		if if_plot
			plot_synthetic_correlation(abs.(ttn_syn_corr_sf .- ed_syn_corr_sf); plot_title="TTN - ED Superfluid")
		end
		@test isapprox(ttn_syn_corr_sf,ed_syn_corr_sf,atol=1e-5)

		# test greens function for the FQH state
		ttn_phys_corr_fqh = physical_correlation(ttn_data["fqh"][2],Lx,Ly; plot_title="TTN FQH", if_plot=false)
		ed_phys_corr_fqh = physical_correlation(ed_data["fqh"][2][1],Lx,Ly; plot_title="ED FQH", if_plot=false)
		if if_plot
			plot_physical_correlation(abs.(ttn_phys_corr_fqh .- ed_phys_corr_fqh); plot_title="TTN - ED FQH")
		end
		@test isapprox(ttn_phys_corr_fqh,ed_phys_corr_fqh,atol=1e-5)

		ttn_syn_corr_fqh = synthetic_correlation(ttn_data["fqh"][2],Lx,Ly; plot_title="TTN FQH", if_plot=false)
		ed_syn_corr_fqh = synthetic_correlation(ed_data["fqh"][2][1],Lx,Ly; plot_title="ED FQH", if_plot=false)
		if if_plot
			plot_synthetic_correlation(abs.(ttn_syn_corr_fqh .- ed_syn_corr_fqh); plot_title="TTN - ED FQH")
		end
		@test isapprox(ttn_syn_corr_fqh,ed_syn_corr_fqh,atol=1e-5)#

	end
end

if false || if_all
	@testset "hopping anisotropy comparison TTN Torus" begin
		ttn_dataloc = get_folder_location("cluster-data/synth-dims/torus")
		ed_dataloc = get_folder_location("cluster-data/exact-diag/torus")

		Lx,Ly = 4,4
		N = 4

		ttn_files = find_data_file(Dict([("layers",Int(log(2,Lx*Ly))),("particles",N),("if_periodic_phys",true),("if_periodic_virt",true),("onsite_strength",0.0)]),"ttn",ttn_dataloc; output_level=0)
		ed_files = find_data_file(Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",0.0)]),"ed",ed_dataloc; output_level=0)

		ttn_data = Dict()
		for f in ttn_files
			data,metadata = read_data_jld2(f,ttn_dataloc; output_level=0)
			ttn_data[string(metadata["hopping_anisotropy"])] = (data["densmat"],metadata["observer"].nrg)
		end

		ed_data = Dict()
		for f in ed_files
			data,metadata = read_data_jld2(f,ed_dataloc; output_level=0)
			ed_data[string(metadata["hopping_anisotropy"])] = (data["densmat"][1],data["nrg"][1])
		end

		# test that the energies match-up
		ed_nrgs = []
		ttn_nrgs = []
		for (key,val) in ed_data
			push!(ed_nrgs,val[2])
			push!(ttn_nrgs,ttn_data[key][2][end])
		end
		if if_plot
			fig = figure()
			xs = parse.(Float64,keys(ed_data))
			scatter(xs,ttn_nrgs,label="TTN")
			scatter(xs,ed_nrgs,label="ED")
			legend()
			xlabel("Hopping Anisotropy")
			ylabel("Energy")
			title("Energy comparison")
		end
		@test isapprox(ed_nrgs,ttn_nrgs,atol=1e-6)

		# test that the density matrices match-up
		for (key,val) in ed_data
			densmat_diff = abs.(ttn_data[key][1]) .- abs.(val[1])
			if if_plot
				fig = figure()
				imshow(abs.(densmat_diff))
				colorbar()
				title("TTN - ED: $key")#
			end
			@test maximum(abs.(densmat_diff)) < 1e-7
		end
	end
end

if false || if_all
	@testset "hopping anisotropy comparison TTN OBC" begin
		ttn_dataloc = get_folder_location("cluster-data/synth-dims/obc")
		ed_dataloc = get_folder_location("cluster-data/exact-diag/obc")

		Lx,Ly = 4,4
		N = 4

		ttn_files = find_data_file(Dict([("layers",Int(log(2,Lx*Ly))),("particles",N),("if_periodic_phys",false),("if_periodic_virt",false),("onsite_strength",0.0)]),"ttn",ttn_dataloc; output_level=0)
		ed_files = find_data_file(Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",false),("if_periodic_y",false),("interaction_strength",0.0)]),"ed",ed_dataloc; output_level=0)

		ttn_data = Dict()
		for f in ttn_files
			data,metadata = read_data_jld2(f,ttn_dataloc; output_level=0)
			ttn_data[string(metadata["hopping_anisotropy"])] = (data["densmat"],metadata["observer"].nrg)
		end

		ed_data = Dict()
		for f in ed_files
			data,metadata = read_data_jld2(f,ed_dataloc; output_level=0)
			ed_data[string(metadata["hopping_anisotropy"])] = (data["densmat"][1],data["nrg"][1])
		end

		# test that the energies match-up
		ed_nrgs = []
		ttn_nrgs = []
		for (key,val) in ed_data
			push!(ed_nrgs,val[2])
			push!(ttn_nrgs,ttn_data[key][2][end])
		end
		if if_plot
			fig = figure()
			xs = parse.(Float64,keys(ed_data))
			scatter(xs,ttn_nrgs,label="TTN")
			scatter(xs,ed_nrgs,label="ED")
			legend()
			xlabel("Hopping Anisotropy")
			ylabel("Energy")
			title("Energy comparison")
		end
		@test isapprox(ed_nrgs,ttn_nrgs,atol=1e-6)

		# test that the density matrices match-up
		for (key,val) in ed_data
			densmat_diff = abs.(ttn_data[key][1]) .- abs.(val[1])
			if if_plot
				fig = figure()
				imshow(abs.(densmat_diff))
				colorbar()
				title("TTN - ED: $key")#
			end
			@test maximum(abs.(densmat_diff)) < 1e-7
		end
	end
end

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


if false || if_all
	#@testset "equivalence of hamiltonians btw TTN and ED" begin
		#for Lx in [4,6]

			nev = 5
				Lx = 6
				lr_dist = Lx-1
				int_stren = 10.0
				wd = "virt"
				if Lx <= 4
					layer_count = 4
				else
					layer_count = 6
				end
				make_smaller_lattice = [Lx,Lx]
				N = Int(Lx/2)
				if_periodic_phys = true
				if_periodic_virt = false

				full_basis = n_particle_basis(N,Lx,Lx; output_level=1)
				lattice_params = Dict([("Lx",Lx),("Ly",Lx),("N",N),("if_periodic_x",if_periodic_phys),("if_periodic_y",if_periodic_virt),("twist_angle",0.0),("full_basis",full_basis)])

				tx = 1.0
				ty = 1.0
				hopping_anisotropy = 1.0
				us = zeros(Float64,Lx)
				for i in 1:lr_dist+1
					us[i] = int_stren
				end
				filling = 0.5
				alpha = N / (Lx*Lx*filling)
				hamilt_params = Dict("alpha"=>alpha,
									"tx"=>tx,
									"ty"=>ty,
									"hopping_anisotropy"=>hopping_anisotropy,
									"U"=>us,
									"which_dir"=>wd,
									"interaction_cutoff"=>1e-5)
				
				ed_ham = buildHam(lattice_params,hamilt_params)

				model_paras = (hopping_anisotropy=hopping_anisotropy,
								restricted_size=make_smaller_lattice,
								if_periodic_phys=if_periodic_phys,
								if_periodic_virt=if_periodic_virt,
								scaling="flat",
								scaling_dist=lr_dist,
								which_dir=wd,
								onsite_strength=int_stren)
				net = build_HH_net(layer_count; syms=true, max_occ=2)

				#ttn_ham_old = long_range_HH_ham_old(net,1.0,alpha; model_paras...)
				#rebuilt_ham_old = rebuild_ed_ham(ttn_ham_old,lattice_params)

				ttn_ham = long_range_HH_ham(net,1.0,alpha; model_paras...,hopping_old=true)
				rebuilt_ham = rebuild_ed_ham(ttn_ham,lattice_params)

				#=fig = figure()
				imshow(real.(Matrix(rebuilt_ham - rebuilt_ham_old)))
				colorbar()
				title("Real Diff")

				fig = figure()
				imshow(imag.(Matrix(rebuilt_ham - rebuilt_ham_old)))
				colorbar()
				title("Imag Diff")=#

				#println("Do the TTN methods make the same Hamiltonian? ",rebuilt_ham == rebuilt_ham_old)
				println("Does the New TTN method match ED? ",rebuilt_ham == ed_ham)
				
				x0 = rand(Float64,size(lattice_params["full_basis"])[2])
        		#everything = eigen(Matrix(rebuilt_ham))
        		#rez = (everything.values,everything.vectors)
				rez = eigsolve(rebuilt_ham,x0,nev,:SR,Lanczos())
				sorted_indices = sortperm(rez[1])
    			states = rez[2][sorted_indices][1:nev]
    			nrgs = rez[1][sorted_indices][1:nev]
				display(nrgs)
				get_occupancy(states[1],lattice_params; plot_title="Correct")

				#=x0_old = rand(Float64,size(lattice_params["full_basis"])[2])
        		rez_old = eigsolve(rebuilt_ham_old,x0_old,nev,:SR,Lanczos())
				sorted_indices_old = sortperm(rez_old[1])
    			states_old = rez_old[2][sorted_indices_old][1:nev]
    			nrgs_old = rez_old[1][sorted_indices_old][1:nev]
				display(nrgs_old)
				get_occupancy(states_old[1],lattice_params; plot_title="Old")=#

				x0_ed = rand(Float64,size(lattice_params["full_basis"])[2])
        		#everything_ed = eigen(Matrix(ed_ham))
        		#rez_ed = (everything_ed.values,everything_ed.vectors)
				rez_ed = eigsolve(ed_ham,x0_ed,nev,:SR,Lanczos())
				sorted_indices_ed = sortperm(rez_ed[1])
    			states_ed = rez_ed[2][sorted_indices_ed][1:nev]
    			nrgs_ed = rez_ed[1][sorted_indices_ed][1:nev]
				display(nrgs_ed)
				get_occupancy(states_ed[1],lattice_params; plot_title="ED, Per Phys=$if_periodic_phys, Per Virt=$if_periodic_virt, LR=$int_stren, Dir=$wd")


				
				#@test ed_ham == rebuilt_ham
		#end
	
	#end
end

#=layers = 4
Lx = Int(sqrt(2^layers))
lr_dist = 0
hopping_anisotropy = 1.0
dataloc = get_folder_location("cluster-data/synth-dims")
params_dict = Dict([("layers",layers),("particles",Int(Lx/2)),("if_periodic_phys",true),("alpha",0.25)])
lattice_params = Dict([("Lx",Lx),("Ly",Lx),("N",params_dict["particles"]),("if_periodic_x",true),("if_periodic_y",false),("full_basis",n_particle_basis(params_dict["particles"],Lx,Lx; output_level=0),("twist_angle",0.0))])
files = find_data_file(params_dict,"ttn",dataloc; output_level=0)

all_hams = Dict()
for f in files
	data,metadata = read_data_jld2(f,dataloc; output_level=0)
	all_hams[string(metadata["onsite_strength"])] = metadata["ham"]
end=#


































"fin"
#using Pkg
#Pkg.activate(".")
using Test

include("../other-funcs/include-other-files.jl")

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
	Ly = lattice_params["Ly"]
	make_smaller_lattice = [Lx,Ly]
	if_synth_rectangle = lattice_params["if_synth_rectangle"]

	interacting_pairs = []
	hopping_pairs = []
	for t1 in TTNKit.terms(ttn_ham)
		coeff = TTNKit.coefficient(t1)
		if length(t1) == 1
			continue
		end
					
		both_terms = TTNKit.terms(t1)

		which_sites = TTNKit.site.(both_terms)
		if if_synth_rectangle
			for i in 1:2
				which_sites[i] = (which_sites[i][2],which_sites[i][1])
			end
		end


		#=if any(which_sites[1] .> make_smaller_lattice[1]) || any(which_sites[2] .> make_smaller_lattice[2])
			display(which_sites)
			continue
		end=#
		ops_here = TTNKit.which_op.(both_terms)
		if ops_here == ["Adag","A"]
			end_site,start_site = linear_index.(which_sites,make_smaller_lattice[1],make_smaller_lattice[2])
			#println("Moving from $start_site to $end_site with coeff $coeff")
			append!(hopping_pairs,[(start_site,end_site,coeff)])
		elseif ops_here == ["Adag * A","Adag * A"]
			s1,s2 = linear_index.(which_sites,make_smaller_lattice[1],make_smaller_lattice[2])
			#println("Interacting at sites $s1 and $s2 with coeff $coeff from $which_sites")
			append!(interacting_pairs,[(s1,s2,coeff)])
		end
	end

	for j in 1:size(full_basis)[2]
		println(round(100*j/size(full_basis)[2],digits=3))
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

function rebuild_ttn_ham(ed_ham,lattice_params::Dict)
	lx,ly = lattice_params["Lx"],lattice_params["Ly"]
	rebuilt_ham = TTNKit.OpSum()
	hopsites = []
	for i in 1:size(ed_ham)[1]
		for j in 1:size(ed_ham)[2]
			if ed_ham[i,j] != 0.0
				if i == j
					basis = lattice_params["full_basis"][:,j]
					println("not including interaction terms")
				else
					basis_start = lattice_params["full_basis"][:,i]
					basis_end = lattice_params["full_basis"][:,j]
					
					allsites = vcat(basis_start,basis_end)
					which_indices = findall(x -> length(findall(y -> y == x,allsites)) == 1,allsites)
					
					starting_site,ending_site = coordinate(allsites[which_indices[1]],lx,ly),coordinate(allsites[which_indices[2]],lx,ly)
					if !in((allsites[which_indices[1]],allsites[which_indices[2]],ed_ham[i,j]),hopsites)
						append!(hopsites,[(allsites[which_indices[1]],allsites[which_indices[2]],ed_ham[i,j])])
						rebuilt_ham += (ed_ham[i,j],"Adag",starting_site,"A",ending_site)
					end
				end
			end
		end
	end
	return rebuilt_ham
end

function rebuild_1deff_to_ed_ham(mps_ham,lattice_params::Dict)
	full_basis = lattice_params["full_basis"]
	rebuilt_ham = spzeros(ComplexF64,size(full_basis)[2],size(full_basis)[2])
	Lx = lattice_params["Lx"]
	Ly = lattice_params["Ly"]
	make_smaller_lattice = [Lx,Ly]
	hopping_pairs = []
	for t1 in TTNKit.terms(mps_ham)
		coeff = TTNKit.coefficient(t1)
					
		all_terms = TTNKit.terms(t1)

		if length(all_terms) == 1
			this_term = all_terms[1]
			synth_sites = [parse(Int,ss[end]) for ss in split.(split(TTNKit.which_op(this_term)," * "),"")]
			phys_site = TTNKit.site(this_term)
			end_site,start_site = linear_index.([(phys_site,synth_sites[1]),(phys_site,synth_sites[2])],make_smaller_lattice[1],make_smaller_lattice[2])
			#println("Found Synth Hopping from $synth_sites on Phys $phys_site with coeff $coeff")
			append!(hopping_pairs,[(start_site,end_site,coeff)])
		else
			phys_sites = TTNKit.site.(all_terms)
			first_operator = TTNKit.which_op(all_terms[1])
			synth_site = parse(Int,split(first_operator,"")[end])
			end_site,start_site = first_operator == "Cr$synth_site" ? linear_index.([(phys_sites[1],synth_site),(phys_sites[2],synth_site)],make_smaller_lattice[1],make_smaller_lattice[2]) : linear_index.([(phys_sites[2],synth_site),(phys_sites[1],synth_site)],make_smaller_lattice[1],make_smaller_lattice[2])
			#println("Found Phys Hopping from $phys_sites on Synth $synth_site with coeff $coeff")
			append!(hopping_pairs,[(start_site,end_site,coeff)])
		end
			

		#=which_sites = TTNKit.site.(both_terms)
		ops_here = TTNKit.which_op.(both_terms)
		if ops_here == ["Adag","A"]
			end_site,start_site = linear_index.(which_sites,make_smaller_lattice[1],make_smaller_lattice[2])
			#println("Moving from $start_site to $end_site with coeff $coeff")
			append!(hopping_pairs,[(start_site,end_site,coeff)])
		elseif ops_here == ["Adag * A","Adag * A"]
			s1,s2 = linear_index.(which_sites,make_smaller_lattice[1],make_smaller_lattice[2])
			#println("Interacting at sites $s1 and $s2 with coeff $coeff from $which_sites")
			append!(interacting_pairs,[(s1,s2,coeff)])
		end=#
	end

	for j in 1:size(full_basis)[2]
		println(round(100*j/size(full_basis)[2],digits=3))
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

	end

	return rebuilt_ham
end

function get_lattice_params_from_ttn_modelparas(model_params::Dict)
	layer_count = model_params[:layers]
	if layer_count % 2 != 0
		Lx = Int(sqrt(2^(layer_count+1)))
		Ly = Int(sqrt(2^(layer_count-1)))
	else
		Lx = Int(sqrt(2^layer_count))
		Ly = Int(sqrt(2^layer_count))
	end
	if model_params[:restricted_size] != [Lx,Ly]
		Lx,Ly = model_params[:restricted_size]
	end

	N = model_params[:particles]
	if_periodic_phys = model_params[:if_periodic_phys]
	if_periodic_virt = model_params[:if_periodic_synth]
	full_basis = n_particle_basis(N,Lx,Lx; output_level=1)
	lattice_params = Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",if_periodic_phys),("if_periodic_y",if_periodic_virt),("twist_angle",0.0),("full_basis",full_basis)])
	return lattice_params
end

function get_lattice_params_from_1deff_modelparas(model_params::Dict)
	lattice_params::Dict{String,Any} = Dict([("Lx",model_params[:L]),("Ly",model_params[:nflavors]),("N",model_params[:nbosons]),("if_periodic_x",model_params[:if_periodic_phys]),("if_periodic_y",model_params[:if_periodic_synth]),("twist_angle",[0.0,0.0])])
	full_basis = n_particle_basis(model_params[:nbosons],model_params[:L],model_params[:nflavors]; output_level=0)
	lattice_params["full_basis"] = full_basis
	return lattice_params
end

if false
	cols = ["b","r","g","m","c"]
    if 10 > length(cols)
        cols = repeat(cols,ceil(Int,10/length(cols)))
    end
	
	anis = 0.7
	intstren = 10.0
	params_dict = Dict([("hopping_anisotropy",anis),("make_smaller_lattice",[6,6]),("if_check_fluxes",true),("particles",3),("layers",6),("filling",0.5),("onsite_strength",intstren),("lr","all"),("if_periodic_phys",false),("if_periodic_synth",true)])
	nev = 10
	model_paras = get_normal_model_params(params_dict)

	this_alpha = 2*model_paras[:particles]/(prod(model_paras[:restricted_size]))#model_paras[:alpha]
	net = build_HH_net_old(model_paras[:layers]; syms=true, max_occ=1)
	ttn_ham = long_range_HH_ham_old(net,1.0,this_alpha; model_paras...)

	lattice_params = get_lattice_params_from_ttn_modelparas(model_paras)
	rebuilt_ham = rebuild_ed_ham(ttn_ham,lattice_params)

	x0 = rand(Float64,size(lattice_params["full_basis"])[2])
	rez = eigsolve(rebuilt_ham,x0,nev,:SR,Lanczos())
	states = rez[2]
	nrgs = rez[1]

	#occs = get_occupancy(states[1],lattice_params; if_plot=true)
	#display(occs)

end

if false || if_all
	@testset "equivalence of hamiltonians btw TTN and ED" begin
		for Lx in [4,6]

			nev = 5
				#Lx = 6
				lr_dist = Lx-1
				int_stren = 0.0
				wd = "virt"
				if Lx <= 4
					layer_count = 4
				else
					layer_count = 6
				end
				make_smaller_lattice = [Lx,Lx]
				N = Int(Lx/2)
				if_periodic_phys = false
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
				x_shift,yshift = !if_periodic_phys, !if_periodic_virt
				alpha = N / ((Lx-x_shift)*(Lx-yshift)*filling)
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

				ttn_ham = long_range_HH_ham(net,1.0,alpha; model_paras...,hopping_old=false)
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
				#println("Does the old TTN method match ED? ",rebuilt_ham_old == ed_ham)
				
				#=x0 = rand(Float64,size(lattice_params["full_basis"])[2])
        		#everything = eigen(Matrix(rebuilt_ham))
        		#rez = (everything.values,everything.vectors)
				rez = eigsolve(rebuilt_ham,x0,nev,:SR,Lanczos())
				sorted_indices = sortperm(rez[1])
    			states = rez[2][sorted_indices][1:nev]
    			nrgs = rez[1][sorted_indices][1:nev]
				display(nrgs)
				get_occupancy(states[1],lattice_params; plot_title="Correct")=#

				#=x0_old = rand(Float64,size(lattice_params["full_basis"])[2])
        		rez_old = eigsolve(rebuilt_ham_old,x0_old,nev,:SR,Lanczos())
				sorted_indices_old = sortperm(rez_old[1])
    			states_old = rez_old[2][sorted_indices_old][1:nev]
    			nrgs_old = rez_old[1][sorted_indices_old][1:nev]
				display(nrgs_old)
				get_occupancy(states_old[1],lattice_params; plot_title="Old")

				x0_ed = rand(Float64,size(lattice_params["full_basis"])[2])
        		#everything_ed = eigen(Matrix(ed_ham))
        		#rez_ed = (everything_ed.values,everything_ed.vectors)
				rez_ed = eigsolve(ed_ham,x0_ed,nev,:SR,Lanczos())
				sorted_indices_ed = sortperm(rez_ed[1])
    			states_ed = rez_ed[2][sorted_indices_ed][1:nev]
    			nrgs_ed = rez_ed[1][sorted_indices_ed][1:nev]
				display(nrgs_ed)
				get_occupancy(states_ed[1],lattice_params; plot_title="ED, Per Phys=$if_periodic_phys, Per Virt=$if_periodic_virt, LR=$int_stren, Dir=$wd")=#


				
				@test ed_ham == rebuilt_ham
		end
	
	end
end

if false
	lxs = [4,6,8,8]
	lys = [4,4,4,8]
	ns = [2,3,4,8]
	cols = ["b","r","g","m","c"]
	if 10 > length(cols)
		cols = repeat(cols,ceil(Int,10/length(cols)))
	end

	thermlim_nrgs = []
	for (idx,N) in enumerate(ns)
		if N < 8
			pdict = Dict([("Lx",lxs[idx]),("Ly",lys[idx]),("N",N),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",0.0),("lr",0),("filling",0.5),("nev",10),("if_save_data",false)])
			lattice_params,hamilt_params,running_args = get_normal_model_params_ed(pdict)    
			basis_dataloc = running_args.basis_dataloc
			filename_dict = make_filename_dict(lattice_params,hamilt_params)
			if_exists,found_data = running_args.if_find_data ? check_data_exists(filename_dict,"ed"; location=running_args.dataloc,output_level=false) : (false,nothing)

			if if_exists
				nrgs = found_data[1]["nrg"]
			else
				full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
				lattice_params["full_basis"] = full_basis
				states,nrgs,rhos,hh = find_eigenstates(running_args.nev,lattice_params,hamilt_params; running_args...)
			end
		else
			numlayers = Int(log(2,lxs[idx]*lys[idx]))
			pdict = Dict([("hopping_anisotropy",1.0),("particles",N),("layers",numlayers),("onsite_strength",0.0)])
			f = find_data_file(pdict,"ttn",get_folder_location("cluster-data/synth-dims/excited-states"); output_level=0)
			data,metadata = read_data_jld2(f[1],get_folder_location("cluster-data/synth-dims/excited-states"); output_level=0)
			nrgs = [metadata["observer"].nrg[end],metadata["observer_1"].nrg[end],metadata["observer_2"].nrg[end],metadata["observer_3"].nrg[end]]
		end

		for i in 1:length(nrgs)
			scatter(lxs[idx]*lys[idx],nrgs[i] - nrgs[1],color=cols[i])
		end
		xlabel("System Size")
		ylabel("Energy - E1")

		if idx == 1
			append!(thermlim_nrgs,nrgs[8])
		elseif idx == 2
			append!(thermlim_nrgs,nrgs[3])
		else
			append!(thermlim_nrgs,nrgs[2])
		end		
	end
	
end

if false
	partcount = 2
	layers = 3
	lx,ly = Int(sqrt(2^(layers+1))),Int(sqrt(2^(layers-1)))
	if_periodic = false
	alpha = partcount / ((lx-!if_periodic)*(ly-!if_periodic))
	stren = 10.0


	us_synth = stren == 0.0 ? zeros(Float64,lx) : stren .* ones(Float64,lx)
	us_synth[1] = stren == 0.0 ? 1.0 : stren
	hamilt_params_synth = Dict("alpha"=>alpha,
                        "flux_direction"=>"x",
                        "tx"=>1.0,
                        "ty"=>1.0,
                        "hopping_anisotropy"=>1.0,
                        "U"=>us_synth,
                        "which_dir"=>"virt",
                        "interaction_cutoff"=>1e-8)

	pdict_synthrect = Dict([("hopping_anisotropy",1.0),("es_count",0),("if_synth_rectangle",true),("particles",partcount),("layers",layers),("mdim",50),("if_save_data",false),("alpha",alpha),("onsite_strength",stren),("lr",stren == 0.0 ? 0 : "all"),("if_periodic_phys",if_periodic),("if_periodic_synth",if_periodic)])
	model_paras_synth = get_normal_model_params(pdict_synthrect)
	net_synth = build_HH_net(model_paras_synth)
	ham_synth = long_range_HH_ham(net_synth,model_paras_synth[:ts],model_paras_synth[:alpha]; model_paras_synth...)
	lattice_params_synth = Dict([("if_synth_rectangle",true),("Lx",ly),("Ly",lx),("N",partcount),("if_periodic_x",pdict_synthrect["if_periodic_synth"]),("if_periodic_y",pdict_synthrect["if_periodic_phys"]),("twist_angle",0.0),("full_basis",n_particle_basis(partcount,ly,lx; output_level=0))])
	rebuilt_ham_synth = rebuild_ed_ham(ham_synth,lattice_params_synth)
	ed_ham_synth = buildHam(lattice_params_synth,hamilt_params_synth)

	println("For Synth Rectangle: ",rebuilt_ham_synth == ed_ham_synth)

	us_phys = stren == 0.0 ? zeros(Float64,ly) : stren .* ones(Float64,ly)
	us_phys[1] = stren == 0.0 ? 1.0 : stren
	hamilt_params_phys = Dict("alpha"=>alpha,
                        "flux_direction"=>"x",
                        "tx"=>1.0,
                        "ty"=>1.0,
                        "hopping_anisotropy"=>1.0,
                        "U"=>us_phys,
                        "which_dir"=>"virt",
                        "interaction_cutoff"=>1e-8)

	pdict_phys = Dict([("hopping_anisotropy",1.0),("es_count",0),("if_synth_rectangle",false),("particles",partcount),("layers",layers),("mdim",50),("if_save_data",false),("alpha",alpha),("onsite_strength",stren),("lr",stren == 0.0 ? 0 : "all"),("if_periodic_phys",if_periodic),("if_periodic_synth",if_periodic)])
	model_paras_phys = get_normal_model_params(pdict_phys)
	net_phys = build_HH_net(model_paras_phys)
	ham_phys = long_range_HH_ham(net_phys,model_paras_phys[:ts],model_paras_phys[:alpha]; model_paras_phys...)
	lattice_params_phys = Dict([("if_synth_rectangle",false),("Lx",lx),("Ly",ly),("N",partcount),("if_periodic_x",pdict_phys["if_periodic_phys"]),("if_periodic_y",pdict_phys["if_periodic_synth"]),("twist_angle",0.0),("full_basis",n_particle_basis(partcount,lx,ly; output_level=0))])
	rebuilt_ham_phys = rebuild_ed_ham(ham_phys,lattice_params_phys)
	ed_ham_phys = buildHam(lattice_params_phys,hamilt_params_phys)

	println("For Phys Rectangle: ",rebuilt_ham_phys == ed_ham_phys)


end

if false
	layers = 2
	N = 2
	lx,ly = Int(sqrt(2^layers)),Int(sqrt(2^layers))
	if_per = true
	tw1 = 0.32
	tw2 = 0.78

	pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",N),("tw1",tw1),("tw2",tw2),("if_periodic_x",if_per),("if_periodic_y",if_per),("hopping_anisotropy",1.0)])
	lattice_params,hamilt_params,running_args = get_normal_model_params_ed(pdict_ed)
	lattice_params["full_basis"] = n_particle_basis(N,lx,ly; output_level=0)
	lattice_params["if_synth_rectangle"] = false
	ham_ed = buildHam(lattice_params,hamilt_params; running_args...)
	#edham_ttn = rebuild_ttn_ham(ham_ed,lattice_params)

	pdict_ttn = Dict([("hopping_anisotropy",1.0),("tw1",tw1),("tw2",tw2),("es_count",0),("particles",N),("layers",layers),("mdim",50),("if_save_data",false),("filling",0.5),("if_periodic_phys",if_per),("if_periodic_synth",if_per)])
	model_paras_ttn = get_normal_model_params(pdict_ttn)
	net = build_HH_net(model_paras_ttn)
	sumham_ttn = long_range_HH_ham(net,model_paras_ttn[:ts],model_paras_ttn[:alpha]; model_paras_ttn...)
	ham_ttn = rebuild_ed_ham(sumham_ttn,lattice_params)

	#ttnsum_matches = all(ed_term in TTNKit.terms(sumham_ttn) for ed_term in TTNKit.terms(edham_ttn))
	#println("TTN Hamilts match: $ttnsum_matches")


	diffmat = round.(ham_ed .- ham_ttn,digits=6)
	display(diffmat)
end

# comparing ED with 1D effective model
if false
	include("../synth-dims/oneD-effective-LR.jl")

	lx = 4
	ly = 4
	N = 2
	if_per = true

	mps_params = Dict([("Lphys",lx),("Lsynth",ly),("if_check_fluxes",false),("particles",N),("if_periodic_phys",if_per),("if_periodic_synth",if_per),("filling",0.5)])
	model_paras = get_1deff_model_params(mps_params)
	mps_ham = hamiltonian(model_paras)

	latparas = get_lattice_params_from_1deff_modelparas(model_paras)
	hamilt_params = Dict([("tx",model_paras[:hopping_anisotropy]),("disorder_strength",0.0),("ty",1.0),("hopping_anisotropy",model_paras[:hopping_anisotropy]),("alpha",model_paras[:alpha]),("U",zeros(ly)),("interaction_cutoff",1e-5),("which_dir","virt"),("flux_direction","x")])

	rebuilt_ham = rebuild_1deff_to_ed_ham(mps_ham,latparas)
	ed_ham = buildHam(latparas,hamilt_params; output_level=0)

	diffham = round.(ed_ham - rebuilt_ham,digits=6)
	display(diffham)

end

































"fin"
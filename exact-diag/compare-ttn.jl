using Pkg
Pkg.activate(".")
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

include_other_files(["other-funcs/data-storage-funcs.jl","synth-dims/long-range-ttn.jl","exact-diag/two-dimensions.jl"])

if_all = true
if_plot = false

if false || if_all
    @testset "Torus 4x4 N=4" begin
        ttn_dataloc = get_folder_location("cluster-data/synth-dims/torus")
        ed_dataloc = get_folder_location("cluster-data/exact-diag")

        Lx,Ly = 4,4
        N = 4

        ttn_files = find_data_file(Dict([("layers",Int(log(2,Lx*Ly))),("particles",N),("if_periodic_phys",true),("if_periodic_virt",true)]),"ttn",ttn_dataloc; output_level=0)
        ed_files = find_data_file(Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",true),("if_periodic_y",true)]),"ed",ed_dataloc; output_level=0)

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

if true || if_all
	@testset "OBC 4x4 N=4" begin
		ttn_dataloc = get_folder_location("cluster-data/synth-dims/obc")
		ed_dataloc = get_folder_location("cluster-data/exact-diag")

		Lx,Ly = 4,4
		N = 4

		ttn_files = find_data_file(Dict([("layers",Int(log(2,Lx*Ly))),("particles",N),("if_periodic_phys",false),("if_periodic_virt",false)]),"ttn",ttn_dataloc; output_level=0)
		ed_files = find_data_file(Dict([("Lx",Lx),("Ly",Ly),("N",N),("if_periodic_x",false),("if_periodic_y",false)]),"ed",ed_dataloc; output_level=0)

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









































"fin"
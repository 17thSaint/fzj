#####################################################
#=

This file contains the simple observable functions for ED

Depends on:
    other-funcs/basic-2d-stuff.jl
    two-dimensions.jl
    observables.jl
    hatsugai-mbcn.jl
    plottings.jl

=#
######################################################

#using Pkg
#Pkg.activate(".")
#Pkg.activate(@__DIR__)
using JLD2

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


include_other_files(["other-funcs/basic-2d-stuff.jl","other-funcs/basic-2d-observables.jl","exact-diag/two-dimensions.jl","exact-diag/observables.jl","exact-diag/hatsugai-mbcn.jl"])
#include_other_files(["exact-diag/time-evolution.jl"])
#include_other_files(["other-funcs/basic-2d-plottings.jl","exact-diag/plottings.jl"])

function make_filename_dict(lattice_params::Dict,hamilt_params::Dict)
    if hamilt_params["U"][2] == 0.0
        intstren = 0.0
    else
        intstren = hamilt_params["U"][1]
    end
    if typeof(hamilt_params["alpha"]) == Vector{Float64}
        alpha_vals = filter(x -> x != 0.0,hamilt_params["alpha"])
        if length(alpha_vals) == 0
            alpha_val = 0.0
        else
            alpha_val = alpha_vals[1]
        end
    else
        alpha_val = hamilt_params["alpha"]
    end

    fdict::Dict{String,Any} = Dict([("Lx",lattice_params["Lx"]),("Ly",lattice_params["Ly"]),("N",lattice_params["N"]),("alpha",alpha_val),("hopping_anisotropy",hamilt_params["tx"]/hamilt_params["ty"]),("interaction_strength",intstren),("if_periodic_x",lattice_params["if_periodic_x"]),("if_periodic_y",lattice_params["if_periodic_y"])])
    if lattice_params["twist_angle"] != [0.0,0.0]
        fdict["twist_angle1"] = lattice_params["twist_angle"][1]
        fdict["twist_angle2"] = lattice_params["twist_angle"][2]
    end
    if hamilt_params["disorder_strength"] != 0.0
        fdict["disorder_strength"] = hamilt_params["disorder_strength"]
    end
    if haskey(hamilt_params,"if_pinning") && hamilt_params["if_pinning"]
        fdict["if_pinning"] = hamilt_params["if_pinning"]
        fdict["pinning_strength"] = hamilt_params["pinning_strength"]
    end
    if hamilt_params["scaling_type"] != "flat"
		fdict["scaling"] = hamilt_params["scaling_type"]
		if hamilt_params["scaling_type"] == "gaussian"
			fdict["sigma"] = hamilt_params["sigma"]
		elseif hamilt_params["scaling_type"] == "exp"
			fdict["corr_length"] = hamilt_params["corr_length"]
		elseif hamilt_params["scaling_type"] == "rydberg"
			fdict["blockade_radius"] = hamilt_params["blockade_radius"]
		else
			error("ULR Scaling Type Not Recognized: $(hamilt_params["scaling_type"])")
		end
	end
    if haskey(hamilt_params,"periodic_potential_strength") && hamilt_params["periodic_potential_strength"] != 0.0
        fdict["periodic_potential_strength"] = hamilt_params["periodic_potential_strength"]
    end

    return fdict
end

function get_lattice_params_from_metadata(metadata::Dict)
    lat_paras = Dict([("Lx",metadata["Lx"]),("Ly",metadata["Ly"]),("N",metadata["N"]),("if_periodic_x",metadata["if_periodic_x"]),("if_periodic_y",metadata["if_periodic_y"]),("full_basis",metadata["full_basis"]),("twist_angle",metadata["twist_angle"])])
    if isnothing(lat_paras["full_basis"])
        lat_paras["full_basis"] = n_particle_basis(lat_paras; output_level=0)
    end
    return lat_paras
end

function get_normal_params_from_lattham(lattice_params::Dict,hamilt_params::Dict,other_params::Dict)
    # get interaction strength
    if hamilt_params["U"][2] == 0.0
        intstren = 0.0
    else
        intstren = hamilt_params["U"][1]
    end

    nev = get(other_params,"nev",10)
    if_save_data = other_params["if_save_data"]
    if_find_data = get(other_params,"if_find_data",false)

    new_params = Dict("Lx"=>lattice_params["Lx"],
                "Ly"=>lattice_params["Ly"],
                "N"=>lattice_params["N"],
                "if_periodic_x"=>lattice_params["if_periodic_x"],
                "if_periodic_y"=>lattice_params["if_periodic_y"],
                "twist_angle"=>lattice_params["twist_angle"],
                "interaction_strength"=>intstren,
                "filling"=>0.5,
                "tx"=>hamilt_params["tx"],
                "ty"=>hamilt_params["ty"],
                "hopping_anisotropy"=>hamilt_params["hopping_anisotropy"],
                "if_reading"=>false,
                "lr"=>"all",
                "nev"=>nev,
                "if_save_data"=>if_save_data,
                "if_find_data"=>if_find_data,
                "periodic_potential_strength"=>hamilt_params["periodic_potential_strength"])

    for (k,v) in other_params
        new_params[k] = v
    end

    return new_params
end

function get_quick_running_args(nev::Int)
    # set running operation parameters    
    running_args::NamedTuple = (nev=nev,
                    if_exact=false,
                    if_function=false,
                    if_reading=false,
                    if_densmat=false,
                    if_find_data=false,
                    if_save_data=false,
                    dataloc="",
                    output_level=0)

    return running_args
end

function get_normal_model_params_ed(params_dict::Dict)
    opl::Int64 = get(params_dict, "output_level", 1)

    # set lattice parameters
    Lx::Int64 = get(params_dict, "Lx", 4)
    Ly::Int64 = get(params_dict, "Ly", Lx)
    N::Int64 = get(params_dict, "N", 2)
    if_periodic_x::Bool = get(params_dict, "if_periodic_x", false)
    if_periodic_y::Bool = get(params_dict, "if_periodic_y", false)
    twist_angle::Vector{Float64} = [get(params_dict, "tw1", 0.0),get(params_dict, "tw2", 0.0)]
    #=if typeof(twist_angle) == String
        println("Twist angle string is $twist_angle")
        tw_str = split(twist_angle,"c")
        tws = tryparse.(Float64,tw_str)
        println("Twist angles found are $tws")
        for (idx,twa) in enumerate(tws)
            if isnothing(twa)
                parts = parse.(Float64,split(tw_str[idx],"p"))
                tws[idx] = sum(parts .* [10.0^(-(i-1)) for i in 1:length(parts)])
            end
        end
        println("Twist angles are $tws")
        twist_angle = tws
    end=#
    expected_dimHilb::Int64 = binomial(Lx*Ly,N)


    opl > 0 ? println("Using ",N," particles with density ",round(N/(Lx*Ly),digits=3)) : nothing

    # build lattice parameters dictionary
    lattice_params::Dict{String,Any} = Dict("Lx"=>Lx,
                        "Ly"=>Ly,
                        "N"=>N,
                        "if_periodic_x"=>if_periodic_x,
                        "if_periodic_y"=>if_periodic_y,
                        "twist_angle"=>twist_angle)

    #println("Finished Building Lattice Parameters")
    # build long range interaction parameters
    stren::Float64 = get(params_dict,"interaction_strength", 0.0)
    lr_dist::Union{String,Int64} = get(params_dict,"lr", "all")
    lr_dist == "all" ? lr_dist = Ly-1 : nothing
    scaling_type::String = get(params_dict,"scaling_type","flat")
    other_params_dict::Dict{String,Any} = Dict([("scaling",scaling_type)])
    if scaling_type != "flat"
        corr_length = get(params_dict,"corr_length",Ly)
        sigma::Float64 = get(params_dict, "sigma", 1.0)
        blockade_radius::Float64 = get(params_dict, "blockade_radius", 1.0)
        magnetic_spacing::Float64 = get(params_dict, "magnetic_spacing", 1.0)
        other_params_dict["corr_length"] = corr_length
        other_params_dict["sigma"] = sigma
        other_params_dict["blockade_radius"] = blockade_radius
        other_params_dict["magnetic_spacing"] = magnetic_spacing
    end
    us::Vector{Float64} = long_range_scaling(lr_dist,Ly,stren; dict_to_symbols(other_params_dict)...)
    interaction_length = scaling_type == "flat" ? lr_dist : corr_length

    # get hopping anisotropy values
    tx::Union{Float64,Nothing} = get(params_dict,"tx",nothing)
    ty::Union{Float64,Nothing} = get(params_dict,"ty",nothing)
    hopping_anisotropy::Float64 = get(params_dict,"hopping_anisotropy",1.0)
    if isnothing(tx) && isnothing(ty)
        if hopping_anisotropy < 1.0
            ty = 1.0 / hopping_anisotropy
            tx = 1.0
        else
            tx = 1.0 * hopping_anisotropy
            ty = 1.0
        end
    end
    #=println("Using Alberto's Hopping Anisotropy")
    tx = hopping_anisotropy
    ty = 1/hopping_anisotropy=#
    
    # build magnetic field parameters
    alpha = get(params_dict,"alpha",nothing)
    if_bc_shift::Bool = get(params_dict,"if_bc_shift",true)
    x_shift,y_shift = if_bc_shift ? (!if_periodic_x, !if_periodic_y) : (0.0,0.0)
    if isnothing(alpha)
        filling::Float64 = get(params_dict,"filling",0.5)
        alpha = N / (filling * (Lx - x_shift) * (Ly - y_shift))
        filling == 0.0 ? alpha = 0.0 : nothing
    end

    # build hamiltonian parameters dictionary and check fluxes for periodicity
    int_cutoff::Float64 = get(params_dict,"interaction_cutoff",1e-5)
    which_dir::String = get(params_dict,"which_dir","virt")

    if typeof(alpha) == Vector{Float64}
        flux_dir::Union{String,Vector{String}} = ["x","y"]
    else
        flux_dir = get(params_dict,"flux_direction","y")
        if if_periodic_y && !if_periodic_x
            flux_dir = "y"
        elseif !if_periodic_y && if_periodic_x
            flux_dir = "x"
        end
    end
    if_check_fluxes::Bool = get(params_dict,"if_check_fluxes",true)
    if_check_fluxes ? flux_dir = check_fluxes(alpha,Lx,Ly,if_periodic_x,if_periodic_y,flux_dir; output_level=opl) : nothing
    if typeof(alpha) != Vector{Float64}
        alpha = [alpha * (flux_dir == "x"), alpha * (flux_dir == "y")]
    end

    pinning_strength::Float64 = get(params_dict,"pinning_strength",1e-3)

    disorder_strength::Float64 = get(params_dict,"disorder_strength",0.0)
    if_pinning::Bool = get(params_dict,"if_pinning",false)

    periodic_potential_strength::Float64 = get(params_dict,"periodic_potential_strength",0.0)

    hamilt_params::Dict{String,Any} = Dict("alpha"=>alpha,
                        "flux_direction"=>flux_dir,
                        "if_pinning"=>if_pinning,
                        "tx"=>tx,
                        "ty"=>ty,
                        "hopping_anisotropy"=>hopping_anisotropy,
                        "disorder_strength"=>disorder_strength,
                        "periodic_potential_strength"=>periodic_potential_strength,
                        "pinning_strength"=>pinning_strength,
                        "U"=>us,
                        "scaling_type"=>scaling_type,
                        "corr_length"=>interaction_length,
                        "which_dir"=>which_dir,
                        "interaction_cutoff"=>int_cutoff)

    
    # set running operation parameters
    filename_dict::Dict{String,Any} = make_filename_dict(lattice_params,hamilt_params)
    filename::String = join(["ed",make_parameters_filename(filename_dict)],"-")
    nev::Int64 = get(params_dict,"nev",1)
    if_save_data::Bool = get(params_dict, "if_save_data", true)
    if if_periodic_x && if_periodic_y
		dataloc::String = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
	elseif if_periodic_x || if_periodic_y
		dataloc = get_folder_location("cluster-data/exact-diag")
	elseif !if_periodic_x && !if_periodic_y
		dataloc = get_folder_location("cluster-data/exact-diag/obc")
    end
    if scaling_type != "flat"
        dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/ulr-length")
	end
    if periodic_potential_strength != 0.0
        dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/periodic-potential")
    end
    if if_pinning
        dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/pinned-scaling")
    end
    dataloc = get(params_dict, "dataloc", dataloc)


    if occursin("geraghty1",dataloc)
        basis_dataloc::String = "/p/project/netenesyquma/geraghty1/data/data-ed/basis-files"
    else
        basis_dataloc = get_folder_location("cluster-data/exact-diag")#dataloc
    end
    if_exact::Bool = get(params_dict, "if_exact", false)
    if_densmat::Bool = get(params_dict, "if_densmat", false)
    if_find_data::Bool = get(params_dict, "if_find_data", true)
    if_function::Bool = get(params_dict, "if_function", false)
    if_reading::Bool = get(params_dict, "if_reading", false)
    running_args::NamedTuple = (nev=nev,
                    filename=filename,
                    if_exact=if_exact,
                    if_function=if_function,
                    if_reading=if_reading,
                    if_densmat=if_densmat,
                    if_find_data=if_find_data,
                    if_save_data=if_save_data,
                    dataloc=dataloc,
                    basis_dataloc=basis_dataloc,
                    output_level=opl)
    
                  
    opl > 0 && println("Finished Building Model")
    return lattice_params,hamilt_params,running_args

end

function run_normal_ed(params_dict::Dict; kwargs...)
    output_level::Int64 = get(kwargs,:output_level,1)

    lattice_params,hamilt_params,running_args = get_normal_model_params_ed(params_dict)    
    basis_dataloc = running_args.basis_dataloc

    # build filename dictionary
    output_level > 0 && display(running_args.filename)
    filename_dict = make_filename_dict(lattice_params,hamilt_params)
    if_exists::Bool,found_data::Union{Vector{Dict},Nothing} = running_args.if_find_data ? check_data_exists(filename_dict,"ed"; location=running_args.dataloc,output_level=output_level-1,if_exact=true,file_type="jld2") : (false,nothing)

    # some old data has bad naming with int_stren = 1.0 even though rest of Us is zeros
    if params_dict["interaction_strength"] == 1.0 && if_exists
        if found_data[2]["U"] != hamilt_params["U"]
            println("Found Data has wrong Interaction Potential")
            if_exists = false
        end
    end

    # check if data exists and rerun if need more eigenstates
    if if_exists
        start_time = time()
        println("Found existing data: ",found_data[2]["filename"])
        if running_args.nev > found_data[2]["nev"]
            running_args.output_level > 0 ? println("Asking for more eigenstates than in file, rerunning") : nothing
            start_time = time()
            full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
            running_args.output_level > 0 ? println("Made basis in ",time()-start_time) : nothing
            lattice_params["full_basis"] = full_basis 
            states,nrgs,rhos = rerun_eigenstates(running_args.nev,lattice_params,hamilt_params,found_data[2],found_data[1]; running_args...)
        elseif running_args.nev < found_data[2]["nev"]
            running_args.output_level > 0 ? println("Asking for fewer eigenstates than in file, using existing data") : nothing
            states = found_data[1]["state"][1:running_args.nev]
            nrgs = found_data[1]["nrg"][1:running_args.nev]
            rhos = found_data[2]["if_densmat"] ? found_data[1]["densmat"][1:running_args.nev] : nothing
            full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
            running_args.output_level > 0 ? println("Made basis in ",time()-start_time) : nothing
            lattice_params["full_basis"] = full_basis 
        else
            states = found_data[1]["state"]
            nrgs = found_data[1]["nrg"]
            rhos = found_data[1]["densmat"]
            full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
            running_args.output_level > 0 ? println("Made basis in ",time()-start_time) : nothing
            lattice_params["full_basis"] = full_basis 
        end
    else

        # make basis only if data doesn't exist
        start_time = time()
        full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
        running_args.output_level > 0 ? println("Made basis in ",time()-start_time) : nothing
        lattice_params["full_basis"] = full_basis 

        # run exact diagonalization for find eigenstates
        if running_args.nev == 1
            states,nrgs,rhos = find_ground_state(lattice_params,hamilt_params; running_args...)
        else
            if running_args.if_function
                states,nrgs,rhos = find_eigenstates(running_args.nev,lattice_params,hamilt_params; running_args...)
            else
                states,nrgs,rhos,hh = find_eigenstates(running_args.nev,lattice_params,hamilt_params; running_args...)
                hamilt_params["H"] = hh
            end
        end
    end

    filepath = running_args.dataloc * "/" * running_args.filename
    filepath = make_sure_file_type(filepath,"jld2")

    if running_args.if_function
        return states,nrgs,rhos,filepath,if_exists,lattice_params,hamilt_params
    else
        return states,nrgs,rhos,filepath,if_exists,lattice_params,hamilt_params
    end

end

#= run data collection with for loops
if false
    
    
    #which_one = args_dict["which_one"]
    #starting_val = (which_one-1)*10 + 1
    #ending_val = which_one*10
    lx,ly,n = 8,4,4
    intstren = 300.0
    #xi = 1.0
    #dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    #anis = 1e-4
    #intstrens = range(0.0,10.0,length=11)
    #intlens = range(0.0,ly,length=11)
    #cols = ["b","g","r","c","y","orange","purple","pink","brown","gray"]
    #for (idx,intstren) in enumerate(intstrens)
    #for (idx,xi) in enumerate(intlens)
    #for (idx,anis) in enumerate(anises)

    
    #=BLAS.set_num_threads(5)
    args_dict = make_args_dict(ARGS)
    intstren = args_dict["interaction_strength"]
    lx = args_dict["Lx"]
    ly = args_dict["Ly"]
    n = args_dict["N"]=#
    #xi = args_dict["corr_length"]


    #dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/pinned-scaling")
    #tws = range(0.0,1.0,length=11)
    #tws2 = range(0.0,1.0,length=3)
    #for (idx,tw1) in enumerate(tws)
    #for (idx2,tw2) in enumerate(tws)
    #for tw1 in tws
    #for ii in 1:1
        #if tw1 == 0.0 && tw2 == 0.0
        #    continue
        #end
        #println("Working on Twist Angle: $(round(tw1,digits=3)) and $(round(tw2,digits=3))")
        #params_dict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
        
        #other_params_dict = make_args_dict(ARGS)
        #intstren = 3.2#other_params_dict["onsite_strength"]
        #xi = 2.0#other_params_dict["corr_length"]
        params_dict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("if_pinning",true),("pinning_strength",1e-2),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",true)])

        #println("Starting from here")

        if true
            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)
        end


        #=fourpt = four_point(states[1],lattice_params; if_plot=false)
        fourpt_2 = four_point(states[2],lattice_params; if_plot=false)
        datadict = Dict([("fourpt_momentum",fourpt),("fourpt_momentum_1",fourpt_2)])
        modify_data(datadict,filepath,"metadata")=#

        #=d,m = read_data(filepath; output_level=0)
        if !haskey(m,"fourpt_momentum")
            fourpt = four_point(states[1],lattice_params; if_plot=false)
            fourpt_2 = four_point(states[2],lattice_params; if_plot=false)
            datadict = Dict([("fourpt_momentum",fourpt),("fourpt_momentum_1",fourpt_2)])
            modify_data(datadict,filepath,"metadata")
        end=#

        #=if idx == 1
            for i in 1:params_dict["nev"]
                scatter(intstren,nrgs[i] - nrgs[1],c=cols[i])
            end
        else
            overlaps = zeros(Float64,params_dict["nev"],params_dict["nev"])
            for i in 1:params_dict["nev"]
                for j in 1:params_dict["nev"]
                    overlaps[i,j] = abs2(dot(states[i],prev_states[j]))
                end
            end
            display(overlaps)

            for i in 1:params_dict["nev"]
                all_overlaps = overlaps[i,:]
                if all(all_overlaps .< 0.01) 
                    scatter(intstren,nrgs[i] - nrgs[1],c="k")
                elseif sort(all_overlaps)[end] / sort(all_overlaps)[end-1] > 10 && maximum(all_overlaps) > 0.1
                    tracked_index = findfirst(x -> all_overlaps[x] == maximum(all_overlaps),1:length(all_overlaps))
                    scatter(intstren,nrgs[i] - nrgs[1],c=cols[tracked_index])
                    cols[i] = cols[tracked_index]
                elseif sum(all_overlaps .> 0.1) > 1
                    scatter(intstren,nrgs[i] - nrgs[1],c="m")
                else
                    # idk whatever
                    scatter(intstren,nrgs[i] - nrgs[1],c="k")
                end
            end
        end

        global prev_states = states=#


        
        #=if params_dict["if_pinning"]
            nrg1_pin = nrgs[1]
            nrg2_pin = nrgs[2]
        else
            nrg1_none = nrgs[1]
            nrg2_none = nrgs[2]
        end=#
        
        #scatter3D(tw1,tw2,nrgs[1] - nrgs[1],c="b")
        #scatter3D(tw1,tw2,nrgs[2] - nrgs[1],c="g")
        #scatter3D(tw1,tw2,nrgs[3] - nrgs[1],c="r")

        #make_density_correlations(states[1],lattice_params; if_save=true,filepath=filepath)
        #=if idx == 1
            get_occupancy(states[1],lattice_params; fix_colorbar=false,plot_title=" ULR=$intstren")
            fig = figure()
        end=#

        #plot_spectrum(intlens,nrgs,idx,params_dict["nev"],"Interaction Length",true; plot_title="")
        #plot_spectrum(intstrens,nrgs,idx,params_dict["nev"],"Interaction Strength",true; plot_title="")
        #plot_spectrum(tws,nrgs,idx,params_dict["nev"],"Theta_x / 2pi",false; plot_title=" V=$intstren")
        #plot_spectrum(anises,nrgs,idx,params_dict["nev"],"Physical Hopping",true; plot_title=" with PP $(lx)x$(ly) N=$n ULR=$intstren")

        #=if true
            if idx == 1
                occs = get_occupancy(states[1],lattice_params; if_plot=true,plot_title=" $(lx)x$(ly) N=$n ULR=$intstren")
                fig = figure()
            else
                occs = get_occupancy(states[1],lattice_params; if_plot=false)
            end
            #display(occs)
            #display(nrgs)
            #println("GS gap is $(nrgs[2] - nrgs[1])")
        end

        
        ftd = ft_density([pi,0],occs)
        #fts[idx,idx2] = real(ftd)
        scatter(intstren,ftd,c="b")
        
        xlabel("Interaction Strength")
        ylabel("FT Density at k=(pi,0)")
        xscale("log")

        if idx == length(intstrens)
            get_occupancy(states[1],lattice_params; if_plot=true,plot_title=" ULR=$intstren")
        end=#

        #=if idx == length(intstrens)
            get_occupancy(states[1],lattice_params; fix_colorbar=false,plot_title=" ULR=$intstren")
        end=#

        
        #gamma1,gamma2,omega = get_hatsugaifull(states[1],states[2],ref_multiplets; if_save=false,filepath=filepath,ref_multis_filenames=[rm1_name,rm2_name])

        #omegas[idx,idx2] = omega
        #gammas1[idx,idx2] = gamma1
        #gammas2[idx,idx2] = gamma2

    #end
    #end

    #=fig = figure()
    imshow(fts; origin="lower",extent=[minimum(tws),maximum(tws),minimum(tws),maximum(tws)])
    colorbar()
    xlabel(L"\theta_x / 2\pi")
    ylabel(L"\theta_y / 2\pi")
    title("FT Density at k=(pi,0) for $(lx)x$(ly) N=$n ULR=$intstren")=#

    #end

    #xlabel(L"\theta_x / 2\pi")
    #ylabel(L"\theta_y / 2\pi")
    #zlabel("Energy Gap")
    #title("Energy Gap Spectrum for $(lx)x$(ly) with $n particles")

    #plot_gamma(collect(tws),collect(tws),gammas1,1)
    #plot_gamma(collect(tws),collect(tws),gammas2,2)
    #plot_omega(collect(tws),collect(tws),omegas)

    #=fig = figure()
    imshow(abs.(omegas))
    colorbar()
    title("Omega Magnitude")

    fig = figure()
    imshow(angle.(omegas) .+ pi; cmap="hsv")
    colorbar()
    title("Omega Phase")

    fig = figure()
    imshow(abs.(gammas1))
    colorbar()
    title("Gamma1 Magnitude")

    fig = figure()
    imshow(abs.(gammas2))
    colorbar()
    title("Gamma2 Magnitude")=#

end=#

#= testing time evolution
if false
    #lx,ly,n = 4,4,2
    anis = 1e-4
    intstren = 0.0
    ppstren = 0.0
    end_tx = 1.0

    args_dict = make_args_dict(ARGS)
    dt_global = args_dict["dt"]
    lx = args_dict["Lx"]
    ly = args_dict["Ly"]
    n = args_dict["particles"]

    params_dict_i = Dict([("output_level",1),("periodic_potential_strength",ppstren),("tx",anis),("ty",1.0),("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
    states_i,nrgs_i,rhos_i,filepath_i,if_found_i,lattice_params_i,hamilt_params_i = run_normal_ed(params_dict_i; output_level=1)

    params_dict_f = Dict([("output_level",1),("periodic_potential_strength",ppstren),("tx",end_tx),("ty",1.0),("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
    states_f,nrgs_f,rhos_f,filepath_f,if_found_f,lattice_params_f,hamilt_params_f = run_normal_ed(params_dict_f; output_level=1)
    
    starting_gs = states_i[1]

    speccount = 1
    dataloc = get_folder_location("cluster-data/exact-diag/time-evo/dt-benchmark")
    time_running_args = (nev=speccount,output_level=1,if_instant_gs=false,if_save_data=true,dataloc=dataloc,)

    tmax_global = 10.0
    ramptime = 2.0
    #dt_global = 0.05

    tevo_params = Dict([ ("tx",(linear_ramp,params_dict_i["tx"],end_tx,ramptime)),("dt",dt_global),("tmax",tmax_global) ])
    tevo_gs,tevo_dict,intspec,saving_args = run_timeevo(starting_gs,tevo_params,lattice_params_i,hamilt_params_i; time_running_args...)
    
    final_fidelity = abs2(dot(tevo_gs[:,end-1],states_f[1]))
    modify_data(Dict([("final_fidelity",final_fidelity)]),joinpath(saving_args[:dataloc],saving_args[:filename]),"metadata"; output_level=2)

end=#



























"fin"
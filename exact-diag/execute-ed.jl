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


include_other_files(["other-funcs/basic-2d-stuff.jl","other-funcs/basic-2d-observables.jl","exact-diag/two-dimensions.jl","exact-diag/observables.jl","exact-diag/hatsugai-mbcn.jl","exact-diag/plottings.jl"])

function make_filename_dict(lattice_params::Dict,hamilt_params::Dict)
    if hamilt_params["U"][2] == 0.0
        intstren = 0.0
    else
        intstren = hamilt_params["U"][1]
    end
    fdict = Dict([("Lx",lattice_params["Lx"]),("Ly",lattice_params["Ly"]),("N",lattice_params["N"]),("alpha",hamilt_params["alpha"]),("hopping_anisotropy",hamilt_params["tx"]/hamilt_params["ty"]),("interaction_strength",intstren),("if_periodic_x",lattice_params["if_periodic_x"]),("if_periodic_y",lattice_params["if_periodic_y"])])
    if lattice_params["twist_angle"] != [0.0,0.0]
        fdict["twist_angle1"] = lattice_params["twist_angle"][1]
        fdict["twist_angle2"] = lattice_params["twist_angle"][2]
    end
    if hamilt_params["disorder_strength"] != 0.0
        fdict["disorder_strength"] = hamilt_params["disorder_strength"]
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

function get_normal_model_params_ed(params_dict::Dict)

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


    # set running operation parameters
    nev::Int64 = get(params_dict,"nev",1)
    if_save_data::Bool = get(params_dict, "if_save_data", true)
    if if_periodic_x && if_periodic_y
		dataloc::String = get_folder_location("cluster-data/exact-diag/torus")
	elseif if_periodic_x || if_periodic_y
		dataloc = get_folder_location("cluster-data/exact-diag")
	elseif !if_periodic_x && !if_periodic_y
		dataloc = get_folder_location("cluster-data/exact-diag/obc")
	end
    dataloc = get(params_dict, "dataloc", dataloc)
    if occursin("geraghty1",dataloc)
        basis_dataloc::String = "/p/project/netenesyquma/geraghty1/data/data-ed/basis-files"
    else
        basis_dataloc = get_folder_location("cluster-data/exact-diag")#dataloc
    end
    opl::Int64 = get(params_dict, "output_level", 1)
    if_exact::Bool = get(params_dict, "if_exact", false)
    if_densmat::Bool = get(params_dict, "if_densmat", false)
    if_find_data::Bool = get(params_dict, "if_find_data", true)
    if_function::Bool = get(params_dict, "if_function", false)
    running_args::NamedTuple = (nev=nev,
                    if_exact=if_exact,
                    if_function=if_function,
                    if_densmat=if_densmat,
                    if_find_data=if_find_data,
                    if_save_data=if_save_data,
                    dataloc=dataloc,
                    basis_dataloc=basis_dataloc,
                    output_level=opl)


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
        corr_length::Int64 = get(params_dict,"corr_length",Ly)
        sigma::Float64 = get(params_dict, "sigma", 1.0)
        blockade_radius::Float64 = get(params_dict, "blockade_radius", 1.0)
        other_params_dict["corr_length"] = corr_length
        other_params_dict["sigma"] = sigma
        other_params_dict["blockade_radius"] = blockade_radius
    end
    us::Vector{Float64} = long_range_scaling(lr_dist,Ly,stren; dict_to_symbols(other_params_dict)...)

    # get hopping anisotropy values
    hopping_anisotropy::Float64 = get(params_dict,"hopping_anisotropy",1.0)
    if hopping_anisotropy < 1.0
		ty::Float64 = 1.0 / hopping_anisotropy
		tx::Float64 = 1.0
	else
		tx = 1.0 * hopping_anisotropy
		ty = 1.0
	end
    #=println("Using Alberto's Hopping Anisotropy")
    tx = hopping_anisotropy
    ty = 1/hopping_anisotropy=#
    
    # build magnetic field parameters
    alpha::Union{Float64,Nothing} = get(params_dict,"alpha",nothing)
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
    flux_dir::String = get(params_dict,"flux_direction","x")
    if if_periodic_y && !if_periodic_x
        flux_dir = "y"
    elseif !if_periodic_y && if_periodic_x
        flux_dir = "x"
    end
    if_check_fluxes::Bool = get(params_dict,"if_check_fluxes",true)
    if_check_fluxes ? flux_dir = check_fluxes(alpha,Lx,Ly,if_periodic_x,if_periodic_y,flux_dir; output_level=opl) : nothing

    disorder_strength::Float64 = get(params_dict,"disorder_strength",0.0)
    hamilt_params::Dict{String,Any} = Dict("alpha"=>alpha,
                        "flux_direction"=>flux_dir,
                        "tx"=>tx,
                        "ty"=>ty,
                        "hopping_anisotropy"=>hopping_anisotropy,
                        "disorder_strength"=>disorder_strength,
                        "U"=>us,
                        "which_dir"=>which_dir,
                        "interaction_cutoff"=>int_cutoff)

    opl > 0 && println("Finished Building Model")
    return lattice_params,hamilt_params,running_args

end

function run_normal_ed(params_dict::Dict; kwargs...)
    output_level::Int64 = get(kwargs,:output_level,1)

    lattice_params,hamilt_params,running_args = get_normal_model_params_ed(params_dict)    
    basis_dataloc = running_args.basis_dataloc

    # build filename dictionary
    filename_dict::Dict{String,Any} = make_filename_dict(lattice_params,hamilt_params)
    filename::String = join(["ed",make_parameters_filename(filename_dict)],"-")
    output_level > 0 && display(filename)
    if_exists::Bool,found_data::Union{Vector{Dict},Nothing} = running_args.if_find_data ? check_data_exists(filename_dict,"ed"; location=running_args.dataloc,output_level=output_level-1) : (false,nothing)

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
            end
        end
    end

    filepath = running_args.dataloc * "/" * filename

    if running_args.if_function
        return states,nrgs,rhos,hh,filepath,if_exists,lattice_params,hamilt_params
    else
        return states,nrgs,rhos,filepath,if_exists,lattice_params,hamilt_params
    end

end

# run data collection with for loops
if false
    
    lx,ly,n = 8,3,3
    #for (idx,n) in enumerate([2,3,4,5])
    intstrens = range(100.0,1000.0,length=10)
    #other_intstrens = range(2.0,10.0,length=37)
    #intstrens = sort([intstrens; other_intstrens])
    #all_nrgs = zeros(Float64,length(thetas))
    #anises = range(1.0,5.0,length=20)
    #nus = range(0.4,0.6,length=100)
    #alphas = range(0.16,2*3/(6*5),length=30)
    #for (idx,alpha) in enumerate(alphas)
    #for (idx,ly) in enumerate(lys)
    #for (idx,nu) in enumerate(nus)
    #for (idx,anis) in enumerate(anises)
    for (idx,intstren) in enumerate(intstrens)
    #for (idx2,sigma) in enumerate(sigmas)
    #for lrd in [0,1]
    #intstren = 0.0

    #= set number of open cores
    open_cores = 5#get(params_dict, "open_cores", 5)
    if typeof(open_cores) != String
        BLAS.set_num_threads(open_cores)
        display(BLAS.get_config())
    end=#

    #intstren = 0.0
    tw2 = 0.0
    tw1 = 0.0
    #tws = range(0.5,1.5,length=21)
    #tws2 = range(0.7,0.8,length=11)
    #omegas::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws))
    #gammas1::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws))
    #gammas2::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws))
    #ref_multiplets,rm1_name,rm2_name = get_reference_multiplets(lx,ly,n; interaction_strength=intstren)
    #cps = zeros(Float64,length(tws))    tws[args_dict["which_twist_angle"]]
    #for (idx,tw1) in enumerate(tws)
    #for (idx2,tw2) in enumerate(tws)
    #for tw1 in tws
    #for ii in 1:1
        #if tw1 == 0.0 && tw2 == 0.0
        #    continue
        #end
        #println("Working on Twist Angle: $(round(tw1,digits=3)) and $(round(tw2,digits=3))")
        params_dict = Dict([("output_level",1),("if_check_fluxes",false),("Lx",lx),("Ly",ly),("N",n),("tw1",tw1),("tw2",tw2),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
        #params_dict = make_args_dict(ARGS)

        states,nrgs,rhos,filepath,if_found = run_normal_ed(params_dict; output_level=1)

        plot_spectrum(intstrens,nrgs,idx,params_dict["nev"],"Interaction Strength",true; plot_title="Theta_y = pi")
        #plot_spectrum(tws,nrgs,idx,params_dict["nev"],"Theta_x / 2pi",false; plot_title=" V=$intstren")

        #=if !if_found
            gamma1,gamma2,omega = get_hatsugaifull(states[1],states[2],ref_multiplets; if_save=true,filepath=filepath,ref_multis_filenames=[rm1_name,rm2_name])
        elseif if_found
            d,m = read_data_jld2(filepath)
            omega = m["omega"]
            gamma1 = m["gamma1"]
            gamma2 = m["gamma2"]
        else
            println("unclear finding")
        end
            
        omegas[idx,idx2] = omega
        gammas1[idx,idx2] = gamma1
        gammas2[idx,idx2] = gamma2=#

    #end
    end

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

end





























"fin"
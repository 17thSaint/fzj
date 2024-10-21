#####################################################
#=

This file contains data comparisons for both TTNs and ED

Depends on:
    synth-dims/long-range-ttn.jl
    review-practice-codes/observables.jl
    review-practice-codes/plottings.jl

    exact-diag/execute-ed.jl
    exact-diag/observables.jl
    exact-diag/plottings.jl

=#
######################################################

include("../other-funcs/include-other-files.jl")
include_other_files(["synth-dims/long-range-ttn.jl","review-practice-codes/observables.jl","review-practice-codes/plottings.jl","exact-diag/execute-ed.jl","exact-diag/observables.jl","exact-diag/plottings.jl"])

function plot_nrg_vs_intstren_fromdata_ttn(layers::Int64,which_strens::Union{String,Vector{Float64}}="all"; kwargs...)
    hanis = get(kwargs, :hopping_anisotropy, 1.0)
    levels_count = get(kwargs, :levels_count, 3)
    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    lx,ly = get_lattice_dims_from_layers(layers)
    n = get(kwargs, :particles, lx)
    pdict = Dict([("layers",layers),("particles",n),("hopping_anisotropy",hanis),("if_periodic_phys",true),("if_periodic_synth",true)])

    if_check = which_strens == "all" ? false : true

    all_files = find_data_file(pdict,"ttn",dataloc)

    nrgs::Dict{String,Vector{Float64}} = Dict([("1",Float64[])])
    for i in 2:levels_count
        nrgs[string(i)] = []
    end
    intstrens = []
    for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        if_done,all_checks = check_nrg_convergence(m)
        if if_check && !(m["onsite_strength"] in which_strens)
            continue
        end
        if !any(collect(values(all_checks)))
            continue
        end
        append!(intstrens,[m["onsite_strength"]])
        for i in 1:levels_count
            keyname = i == 1 ? "observer" : "observer_$(i-1)"
            if keyname in keys(m) && all_checks[string(i-1)]
                append!(nrgs[string(i)],[m[keyname].nrg[end]])
            else
                append!(nrgs[string(i)],[0.0])
            end
        end
    end

    cols = ["b","r","k"]
    fig = figure()
    for i in 1:levels_count
        scatter(intstrens,nrgs[string(i)] .- nrgs["1"],c=cols[i],label="E$(i-1)")
    end
    xlabel("Interaction Strength")
    ylabel("Energy Difference")
    legend()
    ylim([-0.1,0.7])

    title("Energy Spectrum for "*L"\rho_{1D}=1.0"*" TTNs $(lx)x$ly lattice")

end

function plot_nrg_vs_intstren_fromdata_ed(lx::Int64,ly::Int64,which_strens::Union{String,Vector{Float64}}="all"; kwargs...)
    hanis = get(kwargs, :hopping_anisotropy, 1.0)
    levels_count = get(kwargs, :levels_count, 5)
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    n = get(kwargs, :particles, lx)
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",hanis),("if_periodic_x",true),("if_periodic_y",true)])

    if_check = which_strens == "all" ? false : true

    all_files = find_data_file(pdict,"ed",dataloc)

    nrgs::Dict{String,Vector{Float64}} = Dict([("1",Float64[])])
    for i in 2:levels_count
        nrgs[string(i)] = []
    end
    intstrens = []
    for f in all_files
        filename_dict = get_params_dict_from_filename(f)
        if haskey(filename_dict,"twist_angle1") || haskey(filename_dict,"twist_angle2")
            continue
        end
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        if if_check && !(m["U"][end] in which_strens)
            continue
        end
        append!(intstrens,[m["U"][end]])
        for i in 1:levels_count
            if length(d["nrg"]) < i
                append!(nrgs[string(i)],[0.0])
            else
                append!(nrgs[string(i)],[d["nrg"][i]])
            end
        end
    end

    cols = ["b","r","k","k","k"]
    fig = figure()
    for i in 1:levels_count
        scatter(intstrens,nrgs[string(i)] .- nrgs["1"],c=cols[i],label="E$(i-1)")
    end
    xlabel("Interaction Strength")
    ylabel("Energy Difference")
    legend()
    ylim([-0.1,0.7])

    title("Energy Spectrum for "*L"\rho_{1D}=1.0"*" ED $(lx)x$ly lattice")

end


# look at finite size scaling of commensurate filling interaction strength spectrum
if true
    plot_nrg_vs_intstren_fromdata_ttn(6; particles=4)
    plot_nrg_vs_intstren_fromdata_ed(8,8; particles=4)
end































"fin"
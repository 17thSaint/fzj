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

using JLD2,LaTeXStrings
include("../other-funcs/include-other-files.jl")
include_other_files(["synth-dims/long-range-ttn.jl","review-practice-codes/observables.jl","other-funcs/basic-2d-observables.jl","exact-diag/execute-ed.jl","exact-diag/observables.jl"])
#include_other_files(["synth-dims/oneD-effective-LR.jl","synth-dims/plottings-oneD.jl"])
include_other_files(["other-funcs/basic-2d-plottings.jl","review-practice-codes/plottings.jl","exact-diag/plottings.jl"])

function plot_nrg_vs_intstren_fromdata_ttn(layers::Int64,which_strens::Union{String,Vector{Float64}}="all"; kwargs...)
    hanis = get(kwargs, :hopping_anisotropy, 1.0)
    if_gap = get(kwargs, :if_gap, true)
    levels_count = get(kwargs, :levels_count, 3)
    if_plot = get(kwargs, :if_plot, true)
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
        if_done,all_checks = check_nrg_convergence(m,false)
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

    if if_plot
        shift_nrg = if_gap ? nrgs["1"] : zeros(length(nrgs["1"]))
        cols = ["b","r","k"]
        fig = figure()
        for i in 1:levels_count
            scatter(intstrens,nrgs[string(i)] .- shift_nrg,c=cols[i],label="E$(i-1)")
        end
        xlabel("Interaction Strength")
        ylabel("Energy Difference")
        legend()
        if_gap && ylim([-0.1,0.7])

        title("Energy Spectrum for "*L"\rho_{1D}=1.0"*" TTNs $(lx)x$ly lattice")
    end

    return intstrens,nrgs
end

function plot_nrg_vs_intstren_fromdata_ed(lx::Int64,ly::Int64,which_strens::Union{String,Vector{Float64}}="all"; kwargs...)
    hanis = get(kwargs, :hopping_anisotropy, 1.0)
    if_gap = get(kwargs, :if_gap, true)
    levels_count = get(kwargs, :levels_count, 5)
    if_plot = get(kwargs, :if_plot, true)
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    n = get(kwargs, :particles, lx)
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",hanis),("if_periodic_x",true),("if_periodic_y",true)])

    if_check = which_strens == "all" ? false : true

    all_files = find_data_file(pdict,"ed",dataloc)

    display(all_files)

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

    if if_plot
        shift_nrg = if_gap ? nrgs["1"] : zeros(length(nrgs["1"]))
        cols = ["b","r","k","k","k"]
        fig = figure()
        for i in 1:levels_count
            scatter(intstrens,nrgs[string(i)] .- shift_nrg,c=cols[i],label="E$(i-1)")
        end
        xlabel("Interaction Strength")
        ylabel("Energy Difference")
        legend()
        if_gap && ylim([-0.1,0.7])
        xlim([-0.1,2.1])

        title("Energy Spectrum for "*L"\rho_{1D}=1.0"*" ED $(lx)x$ly lattice")
    end
    
    return intstrens,nrgs
end

function get_lattice_params_from_ttn(ttn_metadata::Dict)

    N = ttn_metadata["particles"]
    Lx,Ly = get_lattice_dims_from_layers(ttn_metadata["layers"])
    if_periodic_x = ttn_metadata["if_periodic_phys"]
    if_periodic_y = ttn_metadata["if_periodic_synth"]
    twist_angle = ttn_metadata["twist_angle"]
    full_basis = n_particle_basis(N,Lx,Ly)

    lattice_params::Dict{String,Any} = Dict("Lx"=>Lx,
                        "Ly"=>Ly,
                        "N"=>N,
                        "full_basis"=>full_basis,
                        "if_periodic_x"=>if_periodic_x,
                        "if_periodic_y"=>if_periodic_y,
                        "twist_angle"=>twist_angle)

    return lattice_params
end

function focking_vector_element(wavefunc::TTN.TreeTensorNetwork,local_config::Vector{Int})
    # make states with local configuration
    all_states = ["0" for i in 1:TTN.number_of_sites(wavefunc.net)]
    all_states[local_config] .= "1"

    # build TTN with local configuration
    local_ttn = TTN.ProductTreeTensorNetwork(wavefunc.net,all_states)

    # make sure ortho_center is at top
    TTN.move_ortho!(local_ttn,[TTN.number_of_layers(local_ttn),1])
    TTN.move_ortho!(wavefunc,[TTN.number_of_layers(wavefunc),1])

    # find the overlap
    local_overlap = TTN.inner(local_ttn,wavefunc)

    return local_overlap
end

function focking_vector(wavefunc::TTN.TreeTensorNetwork,full_basis::Matrix{Int64}; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)

    allcoeffs = zeros(ComplexF64,size(full_basis,2))
    for i in 1:size(full_basis,2)
        opl > 0 && println("Working on $i / $(size(full_basis,2))")

        # get the local configuration
        local_config = full_basis[:,i]

        #= make states with local configuration
        all_states = ["0" for i in 1:lx*ly]
        all_states[local_config] .= "1"

        # build TTN with local configuration
        local_ttn = TTN.ProductTreeTensorNetwork(wavefunc.net,all_states)

        # find the overlap
        local_overlap = TTN.inner(local_ttn,wavefunc)=#
        local_overlap = focking_vector_element(wavefunc,local_config)
        allcoeffs[i] = local_overlap
    end

    return allcoeffs
end

function focking_matrix_element(wavefunc::TTN.TreeTensorNetwork,tpo::TTN.MPOWrapper,local_config_left::Vector{Int},local_config_right::Vector{Int}; kwargs...)
    # make states with local configuration
    all_states_left = ["0" for i in 1:TTN.number_of_sites(wavefunc.net)]
    all_states_left[local_config_left] .= "1"

    all_states_right = ["0" for i in 1:TTN.number_of_sites(wavefunc.net)]
    all_states_right[local_config_right] .= "1"

    # build TTN with local configuration
    local_ttn_left = TTN.ProductTreeTensorNetwork(wavefunc.net,all_states_left)
    local_ttn_right = TTN.ProductTreeTensorNetwork(wavefunc.net,all_states_right)

    # make sure ortho_center is at top
    TTN.move_ortho!(local_ttn_left,[TTN.number_of_layers(local_ttn_left),1])
    TTN.move_ortho!(local_ttn_right,[TTN.number_of_layers(local_ttn_right),1])
    TTN.move_ortho!(wavefunc,[TTN.number_of_layers(wavefunc),1])

    # find the overlap
    local_expval = calculate_mpo_expectation(local_ttn_right,local_ttn_left,tpo; kwargs...)

    return local_expval
end

function focking_matrix(wavefunc::TTN.TreeTensorNetwork,tpo::TTN.MPOWrapper,full_basis::Matrix{Int64}; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)

    fock_operator = spzeros(ComplexF64,size(full_basis,2),size(full_basis,2))
    for f in 1:size(full_basis,2)
        opl > 0 && println("Working on $f / $(size(full_basis,2))")
        local_config_left = full_basis[:,f]
        for fp in 1:size(full_basis,2)
            local_config_right = full_basis[:,fp]

            local_element = focking_matrix_element(wavefunc,tpo,local_config_left,local_config_right)
            fock_operator[f,fp] = local_element
        end
    end

    return fock_operator
end

function focking_matrix(which_function::Function,wavefunc::TTN.TreeTensorNetwork,full_basis::Matrix{Int64}; kwargs...)

    tree_operator = which_function(wavefunc, kwargs[:momentum1], kwargs[:momentum2]; kwargs...)

    fockmatrix_operator = focking_matrix(wavefunc,tree_operator,full_basis; kwargs...)

    return fockmatrix_operator
end


#= test 4pt MPO on non-regular TTN
if false
    lx,ly,n = 5,3,3
    intstren = 0.0
    pdict = Dict([("hopping_anisotropy",1.0),("Ly",ly),("if_check_fluxes",false),("Lx",lx),("es_count",0),("expander_fraction",1e-5),("particles",n),("mdim",100),("if_save_data",false),("filling",0.5),("if_find_data",false),("onsite_strength",intstren),("lr","all"),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_states, hamilt, all_obs, all_densmats, all_runtimes = run_synth_dims_generic(pdict)

    m1 = [0.0,1/ly]
    m2 = [0.0,0/ly]
    #fourpt_rez = four_point(all_states,m1,m2,[lx,ly])
    fourpt_rez = four_point(all_states,[lx,ly]; if_plot=true, plot_title=" TTN")

    pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("nev",10),("hopping_anisotropy",1.0),("if_check_fluxes",false),("if_save_data",false),("if_find_data",false)])
    states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict_ed; output_level=1)

    #fourpt_ed = ft_fourpt(states[1],m1,m2,lattice_params)
    fourpt_ed = four_point(states[1],lattice_params; if_plot=true,plot_title=" ED")

    println("TTN fourpt: ",fourpt_rez)
    println("ED fourpt: ",fourpt_ed)
    println("If equal: ",isapprox(fourpt_rez,fourpt_ed,atol=1e-5))
end=#

#= plot rho1D NRG spectrum transition scaling in size
if false
    cols = ["b","r","k"]
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    all_files = find_data_file(Dict([("hopping_anisotropy",1.0),("if_periodic_x",true),("if_periodic_y",true)]),"ed",dataloc; file_type="jld2")
    
    g2s = []
    g3s = []
    ulrs = []
    lxs = []
    firstgaps = Dict()
    for f in all_files
        fparas = get_params_dict_from_filename(f)
        if fparas["Ly"] == 2*fparas["Lx"] && fparas["N"] == fparas["Lx"]
            d,m = read_data_jld2(joinpath(dataloc,f); output_level=0)
            lx = m["Lx"]
            #scatter(m["U"][end],0.0,c=cols[lx-2])
            #scatter(m["U"][end],d["nrg"][2] - d["nrg"][1],c=cols[lx-2])
            #scatter(m["U"][end],d["nrg"][3] - d["nrg"][1],c=cols[lx-2])

            append!(g2s,[d["nrg"][2] - d["nrg"][1]])
            append!(g3s,[d["nrg"][3] - d["nrg"][1]])
            append!(ulrs,[m["U"][end]])
            append!(lxs,[lx])

            if m["U"][end] == 0.0
                firstgaps[string(lx)] = d["nrg"][3] - d["nrg"][1]
            end
        end
    end

    for (idx,lx) in enumerate(lxs)
        col = cols[lx - 2]
        scatter(ulrs[idx],0.0,c=col)
        scatter(ulrs[idx],g2s[idx]/firstgaps[string(lx)],c=col)
        scatter(ulrs[idx],g3s[idx]/firstgaps[string(lx)],c=col)
    end
    scatter(0.0,0.0,c=cols[1],label="N=3")
    scatter(0.0,0.0,c=cols[2],label="N=4")
    scatter(0.0,0.0,c=cols[3],label="N=5")

    xlabel("Interaction Strength")
    ylabel("E - E0")
    title("Energy Spectrum Finite Size Scaling rho_1D = 1.0")
    legend()

end=#

#= see if excited TTNs restrict_size matches ED
if false
    lx,ly,n = 3,3,2
    layers = Int(ceil(log(2,lx*ly)))

    pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",0.0),("nev",10),("hopping_anisotropy",1.0),("if_check_fluxes",false),("if_save_data",false),("if_find_data",false)])
    states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict_ed; output_level=1)
    display(nrgs)

    pdict_ttn = Dict([("hopping_anisotropy",1.0),("make_smaller_lattice",[lx,ly]),("es_count",2),("if_check_fluxes",false),("expander_fraction",1.0),("particles",n),("layers",layers),("mdim",200),("if_save_data",false),("filling",0.5),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_states, hamilt, all_obs, all_densmats, all_runtimes = run_synth_dims_generic(pdict_ttn)

end=#

#= energy scaling Lx with fixed Ly
if false
    intstren = 300.0
    ns = [3,4,5,6,7,8]
    splittings = zeros(Float64,length(ns))
    gaps = zeros(Float64,length(ns))
    for n in ns
        lx,ly = 2*n,4
        
        if n < 6    
            pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("nev",10),("hopping_anisotropy",1.0),("if_save_data",false),("if_find_data",true)])
            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict_ed; output_level=1)
            splittings[n-2] = nrgs[2] - nrgs[1]
            gaps[n-2] = nrgs[3] - nrgs[1]
        else
            dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
            pdict_ttn = Dict([("Lx",lx),("Ly",ly),("particles",n),("onsite_strength",intstren),("hopping_anisotropy",1.0),("if_periodic_phys",true),("if_periodic_synth",true)])
            all_files = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
            length(all_files) == 0 && continue
            d,m = read_data(joinpath(dataloc_ttn,all_files[1]); output_level=0)
            splittings[n-2] = abs(m["energies_1"][end] - m["energies"][end])
            if haskey(m,"energies_2")
                gaps[n-2] = abs(m["energies_2"][end] - m["energies"][end])
            end
        end
    end

    fig = figure()
    plot(2 .* ns,splittings,"p",c="b")
    xlabel("Lx")
    ylabel("E2 - E1")
    title("Topological Degeneracy Splitting Scaling ULR=$intstren")
    yscale("log")

    fig = figure()
    plot(2 .* ns,gaps,"p",c="r")
    xlabel("Lx")
    ylabel("E3 - E1")
    title("Topological Gap Scaling ULR=$intstren")
    ylim([-0.05,1.1*maximum(gaps)])
end=#

# max fourpt diag transition scaling Lx with fixed Ly
if true
    ns = [4,6,7,8]
    cols = ["b","r","k","g"]
    ly = 4
    maxdiags_dict = Dict{String,Vector{Float64}}()
    visibs_dict = Dict{String,Vector{Float64}}()
    intstrens_dict = Dict{String,Vector{Float64}}()
    for (idx,n) in enumerate(ns)
        lx = 2*n
        intstrens = []
        maxdiags = []
        visibs = []
        if n < 6
            dataloc_ed = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
            pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
            all_files_ed = find_data_file(pdict_ed,"ed",dataloc_ed; file_type="jld2")
            for f in all_files_ed
                d,m = read_data(joinpath(dataloc_ed,f); output_level=0)
                if haskey(m,"fourpt_momentum")
                    append!(intstrens,[m["U"][end]])
                    append!(maxdiags,[maximum(diag(m["fourpt_momentum"]))])
                    append!(visibs,visibility_fourpt(m["fourpt_momentum"]))
                elseif haskey(m,"fourpt_momentum_diag")
                    append!(intstrens,[m["U"][end]])
                    append!(maxdiags,[maximum(m["fourpt_momentum_diag"])])
                end
            end
        else
            dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
            pdict_ttn = Dict([("Lx",lx),("Ly",ly),("particles",n),("hopping_anisotropy",1.0),("if_periodic_phys",true),("if_periodic_synth",true)])
            all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
            for f in all_files_ttn
                d,m = read_data(joinpath(dataloc_ttn,f); output_level=0)
                if haskey(m,"fourpt_momentum")
                    append!(intstrens,[m["onsite_strength"]])
                    if haskey(m,"fourpt_momentum_1")
                        append!(maxdiags,[0.5*maximum(diag(m["fourpt_momentum"] + m["fourpt_momentum_1"]))])
                        append!(visibs,visibility_fourpt(m["fourpt_momentum"]))
                    else
                        append!(maxdiags,[maximum(diag(m["fourpt_momentum"]))])
                        append!(visibs,visibility_fourpt(m["fourpt_momentum"]))
                    end
                end
            end                    
        end
        println("N=$n, Intstrens: ",intstrens)
        length(intstrens) == 0 && continue
        
        #=laughlin = maxdiags[findfirst(x -> intstrens[x] == 0.0,1:length(intstrens))]
        scatter(intstrens,maxdiags ./ laughlin,c=cols[idx],label="N=$n")
        ylabel("Max["*L"a_{k}^{\dagger}a_{k}^{\dagger}a_{k}a_{k}"*"]")=#

        laughlin = visibs[findfirst(x -> intstrens[x] == 0.0,1:length(intstrens))]
        scatter(intstrens,visibs ./ laughlin,c=cols[idx],label="N=$n")
        ylabel("Visibility, min/max (normalized)")
        
        xlabel("Interaction Strength")
        title("k-DW Transition Experimental Finite Size Scaling")
        legend()
        maxdiags_dict[string(n)] = maxdiags
        intstrens_dict[string(n)] = intstrens
        visibs_dict[string(n)] = visibs
    end
end#

#= look at finite size scaling of commensurate filling interaction strength spectrum
if true
    eight_xs, eight_ys = plot_nrg_vs_intstren_fromdata_ttn(6; particles=8, if_gap=true,if_plot=false)
    four_xs, four_ys = plot_nrg_vs_intstren_fromdata_ed(4,8; particles=4, if_gap=true,if_plot=false)
    three_xs, three_ys = plot_nrg_vs_intstren_fromdata_ed(3,8; particles=3, if_gap=true,if_plot=false)
    five_xs, five_ys = plot_nrg_vs_intstren_fromdata_ed(5,8; particles=5, if_gap=true,if_plot=false)

    eight_23gap = eight_ys["3"] .- eight_ys["2"]
    four_23gap = four_ys["3"] .- four_ys["2"]
    three_23gap = three_ys["3"] .- three_ys["2"]
    five_23gap = five_ys["3"] .- five_ys["2"]

    cols = get_cmap("plasma")
    divisor = 7

    fig = figure()
    plot(three_xs,three_23gap,"-p",c=cols(4/divisor),label="ED 3x8")
    plot(four_xs,four_23gap,"-p",c=cols(3/divisor),label="ED 4x8")
    plot(five_xs,five_23gap,"-p",c=cols(2/divisor),label="ED 5x8")
    plot(eight_xs,eight_23gap,"-p",c=cols(1/divisor),label="TTN 8x8")
    xlim([-0.1,2.1])
    ylim([-0.1,0.7])
    legend()
    xlabel("Interaction Strength")
    ylabel("E3 - E2")
end=#

#= check the NRGs of synth-rectangle TTNs vs ED
if false
    ttn_intstrens,ttn_nrgs = plot_nrg_vs_intstren_fromdata_ttn(5; particles=4, if_plot=false)
    ed_intstrens,ed_nrgs = plot_nrg_vs_intstren_fromdata_ed(4,8; particles=4, if_plot=false)
    
    for level in 1:length(keys(ttn_nrgs))
        fig = figure()
        plot(ed_intstrens,ed_nrgs[string(level)],c="r",label="ED")
        scatter(ttn_intstrens,ttn_nrgs[string(level)],c="b",label="TTN")
        xlabel("Interaction Strength")
        ylabel("Energy Difference")
        legend()
        title("Energy Spectrum for "*L"\rho_{1D}=1.0"*" E$(level-1)")
    end
end

# build phase diagram for ULR vs rho1D at finite and infinite ULR from flatness
function build_phase_diagram_ulr_rho1d_flatness()
    configs = [(8,3,3),(8,4,4),(4,6,3),(3,8,3),(8,5,5),(8,7,7)]

    ulrs::Vector{Float64} = Float64[]
    flatnesses::Vector{Float64} = Float64[]
    oneDrhos::Vector{Float64} = Float64[]

    for (idx,config) in enumerate(configs)
        lx,ly,n = config
        append!(flatnesses,[twist_flatness_1deff(lx,ly,n; if_plot_spectrum=true,plot_title=" $(lx)x$ly N=$n")])
        append!(ulrs,[400.0])
        append!(oneDrhos,[n/lx])
    end

    #=for (lx,ly,n) in configs
        if (lx,ly,n) == (8,7,7)
            continue
        end
        local_strens,local_flats = twist_flatness_ed(lx,ly,n; if_plot=true, max_intstren=500.0)
        append!(ulrs,local_strens)
        append!(flatnesses,local_flats)
        append!(oneDrhos,ones(Float64,length(local_strens)) .* (n / lx))
    end=#

    bin_count = 100
    data_dict = bin_values(flatnesses,bin_count)
    bv = [data_dict[val] for val in flatnesses]
    min_nrgs2, max_nrgs2 = 0.0,1.0
    normalized_bv = [(val - minimum(bv)) / (maximum(bv) - minimum(bv)) * (max_nrgs2 - min_nrgs2) + min_nrgs2 for val in bv]

    #=fig = figure()
    scatter(oneDrhos, ulrs, c=normalized_bv, cmap="viridis")
    colorbar()
    xlabel(L"\rho_{1D}")
    ylabel("ULR")
    title("Flatness Phase Diagram")
    ylim([-0.1,4.1])=#

    plot_phasediag_ulrrho1d_flatness(oneDrhos,ulrs,bv,500.0)

end

# build_phase_diagram_ulr_rho1d for ft dd
if false
    configs = [(8,3,3),(4,6,3),(3,8,3),(8,4,4),(8,5,5)]

    ulrs::Vector{Float64} = Float64[]
    ftdds::Vector{Float64} = Float64[]
    oneDrhos::Vector{Float64} = Float64[]

    for (lx,ly,n) in configs
        local_strens,local_ftdd = findall_ft_dd(lx,ly,n; if_plot=false)
        append!(ulrs,local_strens)
        append!(ftdds,1 ./ local_ftdd)
        append!(oneDrhos,ones(Float64,length(local_strens)) .* (n / lx))
    end

    #=for (idx,config) in enumerate(configs)
        lx,ly,n = config
        append!(ftdds,[findall_ft_dd_1deff(lx,ly,n)])
        append!(ulrs,[900.0])
        append!(oneDrhos,[n/lx])
    end=#

    bin_count = 100
    data_dict = bin_values(ftdds,bin_count)
    bv = [data_dict[val] for val in ftdds]
    min_nrgs2, max_nrgs2 = minimum(ftdds),maximum(ftdds)
    normalized_bv = [(val - minimum(bv)) / (maximum(bv) - minimum(bv)) * (max_nrgs2 - min_nrgs2) + min_nrgs2 for val in bv]

    #=fig = figure()
    scatter(oneDrhos, ulrs, c=normalized_bv, cmap="plasma")
    colorbar()
    xlabel(L"\rho_{1D}")
    ylabel("ULR")
    title("FT-DD Ratio Phase Diagram")
    yscale("log")=#

    plot_phasediag_ulrrho1d_flatness(oneDrhos,ulrs,normalized_bv,1000.0)
    title("Inverse FT DD Max Phase Diagram")
end

# check perturbative 1Deff against ED # it seems this only works for rho1D = 1/2 and 1.0
if false
    lx,ly,n = 8,3,3

    intstrens = vcat(range(1000.0,10000.0,length=21),range(0.0,100.0,length=21))
    cols = get_colors(3)
    for (idx,intstren) in enumerate(intstrens)
        if intstren > 100.0
            pdict_mps = Dict([("Lphys",lx),("Lsynth",ly),("particles",n),("interaction_strength",intstren),("if_remapping",false),("es_count",2),("nrgtol",1e-6),("mdim",200),("if_periodic_phys",true),("if_periodic_synth",true),("filling",0.5),("if_find_data",false),("if_save_data",false)])
            psis,rhos,nrgs,mparas,if_found = run_normal_1deffmps(pdict_mps)
            result = ft_densitydensity_correlation(pi/2,psis[1])
        else
            pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
            psis,nrgs,rhos,filepath_ed,if_exists_ed,latparas,hamparas = run_normal_ed(pdict_ed)
            result = ft_densitydensity_correlation(pi/2,psis[1],latparas)
        end

        scatter(intstren,result,c="b")
        xscale("log")
        xlabel("Interaction Strength")
        ylabel("FT-DD pi/2")

    end
end

# checking ED vs full infinite limit 1Deff # need larger intstren for (8,5,5)
if false
    cols = get_colors(3)
    configs = [(8,3,3),(8,4,4),(4,6,3),(3,8,3),(8,5,5)]
    for (lx,ly,n) in configs
        #lx,ly,n = 8,3,3
        pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
        all_files_ed = find_data_file(pdict_ed,"ed",get_folder_location("cluster-data/exact-diag/torus"))
        filter!(x -> !occursin("twist",x),all_files_ed)
        fig = figure()
        intstrens_ed = Float64[]
        for f in all_files_ed
            d,m = read_data_jld2(get_folder_location("cluster-data/exact-diag/torus") * "/" * f; output_level=0)
            nrgs_ed = d["nrg"]
            append!(intstrens_ed,[m["U"][end]])
            for i in 1:3
                scatter(m["U"][end],nrgs_ed[i],c=cols[i])
            end
        end

        pdict_mps = Dict([("Lphys",lx),("Lsynth",ly),("nbosons",n),("if_periodic_phys",true),("if_periodic_synth",true)])
        all_files = find_data_file(pdict_mps,"mps",get_folder_location("cluster-data/synth-dims/excited-states"))
        d,m = read_data_jld2(get_folder_location("cluster-data/synth-dims/excited-states") * "/" * all_files[1]; output_level=0)
        nrgs_mps = zeros(Float64,3)
        for i in 1:3
            obs_string = i == 1 ? "observer" : "observer_$(i-1)"
            nrgs_mps[i] = m[obs_string].energies[end]
        end
        sort!(nrgs_mps)
        for i in 1:3
            plot([0.0,maximum(intstrens_ed)],[nrgs_mps[i],nrgs_mps[i]],c=cols[i],label="E$i")
        end
        legend()

        xlabel("Interaction Strength")
        ylabel("Energy")
        title("ED (scatter) vs 1Deff (line) for $(lx)x$ly N=$n")
    end

    
    
end

# energy spectrum scaling with Lx/Ly = 2.0
if false
    configs = [(6,3,3),(8,4,4),(10,5,5)]
    for config in configs
        lx,ly,n = config
        if lx == 10
            filename = "../cluster-data/exact-diag/torus/ed-Ly-5-interaction_strength-0.0-Lx-10-if_pinning-true-N-5-alpha-0.2-if_periodic_x-true-if_periodic_y-true-hopping_anisotropy-1.0.jld2"
            d,m = read_data_jld2(filename; output_level=0)
            nrgs = d["nrg"]
        else
            pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",0.0),("filling",0.5),("lr","all"),("if_find_data",true),("if_save_data",false),("nev",10)])
            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)
        end

        scatter(lx,nrgs[1] - nrgs[1],c="b")
        scatter(lx,nrgs[2] - nrgs[1],c="g")
        scatter(lx,nrgs[3] - nrgs[1],c="r")

        xlabel("Lx")
        ylabel("Energy Gap")
        title("Energy Spectrum Finite Size Scaling")
    end

    ttn_configs = [(16,8,8)]
    for config in ttn_configs
        lx,ly,n = config
        layers = get_layers_from_latticesize(lx,ly)
        pdict = Dict([("layers",layers),("particles",n),("hopping_anisotropy",1.0),("if_periodic_phys",true),("if_periodic_synth",true),("onsite_strength",0.0)])
        dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
        all_files = find_data_file(pdict,"ttn",dataloc)
        d,m = read_data_jld2(dataloc * "/" * all_files[1]; output_level=0)
        nrgs::Dict{String,Float64} = Dict([("1",m["observer"].nrg[end])])

        for i in 2:3
            keyname = "observer_$(i-1)"
            if keyname in keys(m)
                nrgs[string(i)] = m[keyname].nrg[end]
            else
                nrgs[string(i)] = 0.0
            end
        end

        scatter(lx,nrgs["1"] - nrgs["1"],c="b")
        scatter(lx,nrgs["2"] - nrgs["1"],c="g")
        scatter(lx,nrgs["3"] - nrgs["1"],c="r")

        xlabel("Lx")
        ylabel("Energy Gap")
        title("Energy Spectrum Finite Size Scaling")
    end
end

# check the FT-D scaling with Lx/Ly = 2.0
if false
    limit_vals = Float64[]
    limit_sizes = Int64[]

    # ED section
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    ks = range(0,2*pi,length=50)
    configs = [(10,5,5),(8,4,4)]#(6,4,3),(8,4,4),(10,4,5)]
    for (lx,ly,n) in configs
        if_pinning = true
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
        all_files = find_data_file(pdict,"ed",dataloc; output_level=0)

        filter!(x -> !occursin("twist_angle1",x),all_files)
        filter!(x -> !occursin("mk",x),all_files)
        filter!(x -> occursin("if_pinning",x),all_files)

        display(all_files)

        cdw_sfs = Float64[]
        intstrens = Float64[]
        if_plot = false
        for f in all_files

            filename_dict = get_params_dict_from_filename(f)
            if lx == 10 && filename_dict["interaction_strength"] == 0.8
                continue
            end
            d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

            if lx == 8 && round(m["U"][end],digits=0) != m["U"][end]
                continue
            end

            append!(intstrens,m["U"][end])

            latparas = get_lattice_params_from_metadata(m)
            occs = get_occupancy(d["state"][1],latparas; if_plot=if_plot,plot_title="$(lx)x$(ly) n=$n ULR=$(m["U"][end])")
            ftd = ft_density([pi,0],occs)
            append!(cdw_sfs,[abs(ftd)])

            if m["U"][end] == 1000.0
                push!(limit_vals,abs(ftd))
                push!(limit_sizes,lx)
            end

        end

        #fig = figure()
        scatter(intstrens,cdw_sfs,label="$(lx)x$(ly)")
        xlabel("Interaction Strength")
        ylabel("CDW SF at k=(pi,0)")
        title("CDW SF at "*L"\rho_{1D}=1/2")
        xscale("log")
        legend()
    end

    # TTN section
    dataloc = get_folder_location("cluster-data/synth-dims/torus")
    configs = [(7,8),(7,6)]
    for (layers,particles) in configs
        pdict = Dict([("layers",layers),("particles",particles),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
        all_files = find_data_file(pdict,"ttn",dataloc)
        display(all_files)

        intstrens = Float64[]
        ftvals = ComplexF64[]
        for f in all_files
            d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

            #= checking convergence of energies
            all_convs = [false,false,false]
            for i in 1:3
                observerkey = i == 1 ? "observer" : "observer_$(i-1)"
                if observerkey in keys(m)
                    all_convs[i] = abs(m[observerkey].nrg[end] - m[observerkey].nrg[end-1]) < m["nrgtol"]
                end
            end
            println("For $(m["onsite_strength"]) all convs are ",all_convs)=#

            lx,ly = 2*particles,particles
            occs = get_occupancy(d["densmat"]; plot_title="Intstren = $(m["onsite_strength"])",if_plot=false)[1:lx,1:ly]
            ftval = ft_density([pi,0.0],occs)
            push!(intstrens,m["onsite_strength"])
            push!(ftvals,ftval)

            if m["onsite_strength"] == 1000.0
                push!(limit_vals,abs(ftval))
                push!(limit_sizes,particles*2)
            end
        end
        scatter(intstrens,abs.(ftvals),label="$(2*particles)x$particles")
        legend()
    end


    fig = figure()
    scatter(limit_sizes,limit_vals)
    xlabel("Lx")
    ylabel("CDW SF at k=(pi,0)")
    title("CDW SF at "*L"\rho_{1D}=1/2"*", ULR=1000.0")
    yscale("log")

end

# finite size scaling of FT-DD at ULR = 0.0
if false
    intstren = 100.0

    # ED version with 8x4 and 10x5
    configs = [(8,4,4),(10,5,5)]
    for config in configs
        lx,ly,n = config
        dataloc = get_folder_location("cluster-data/exact-diag/torus")
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
        all_files = find_data_file(pdict,"ed",dataloc; output_level=0)
        filter!(x -> !occursin("twist_angle1",x),all_files)
        filter!(x -> !occursin("mk",x),all_files)
        filter!(x -> !occursin("pinning",x),all_files)
        f = all_files[1]
        d,m = read_data_jld2(dataloc * "/" * f; output_level=1)
        dds = m["densitydensity"]
        ft_dds_stripe = abs(ft_densitydensity([pi,0],dds))
        scatter(lx,ft_dds_stripe,c="r")
    end

    # TTN version with 16x8 and 12x6
    configs = [8,6]
    for n in configs
        dataloc = get_folder_location("cluster-data/synth-dims/torus")
        pdict = Dict([("layers",7),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
        all_files = find_data_file(pdict,"ttn",dataloc)
        f = all_files[1]
        d,m = read_data_jld2(dataloc * "/" * f; output_level=1)
        dds = m["densitydensity"]
        ft_dds_stripe = abs(ft_densitydensity([pi,0],dds))
        scatter(n*2,ft_dds_stripe,c="b")
    end

    xlabel("Lx")
    ylabel("FT-DD at k=(" * L"\pi" * ",0)")
    title("Finite Size Scaling of FT-DD at ULR=$intstren")
    yscale("log")
end

# testing 2pt MPO FT constructor against known density matrix FT
if false
    lx,ly,n = 4,4,2
    params_dict = Dict([("hopping_anisotropy",1.0),("if_check_fluxes",false),("es_count",0),("expander_fraction",0.5),("particles",n),("layers",Int(log(2,lx*ly))),("mdim",100),("if_save_data",false),("alpha",0.0),("onsite_strength",0.0),("lr",0),("if_periodic_phys",true),("if_periodic_synth",true)])
    #psi, hamilthere, obs, rho, rt = run_synth_dims_generic(params_dict)

    ks = range(0,1.0,length=16)
    vals = zeros(Float64,length(ks),length(ks))
    for (idx,kx) in enumerate(ks)
        for (idx2,ky) in enumerate(ks)
            println("Working on kx = $kx, ky = $ky")
            rho = two_point_mpo(psi; momentum = [kx,ky])
            vals[idx2,idx] = real(calculate_mpo_expectation(psi,rho))
        end
    end
    fig = figure()
    imshow(vals,origin="lower",extent=[0,1,0,1])
    colorbar()
    xlabel("kx")
    ylabel("ky")

    rho_real = density_matrix(psi)
    dm_vals = zeros(Float64,length(ks),length(ks))
    for (idx,kx) in enumerate(ks)
        for (idx2,ky) in enumerate(ks)
            println("Working on kx = $kx, ky = $ky")
            dm_vals[idx2,idx] = real(ft_density_matrix(rho_real,[kx,ky],lx,ly))
        end
    end
    fig = figure()
    imshow(dm_vals,origin="lower",extent=[0,1,0,1])
    colorbar()
    xlabel("kx")
    ylabel("ky")

end

# test zigzag curve
if false
    lx,ly = 8,4
    zc = zigzag_curve(lx,ly)
    display(zc)
end=#

#= read 4pt momentum MPO
if false
    lx,ly,n = 16,8,8
    layers = Int(log(2,lx*ly))
    intstren = 0.0

    #=pdict_ttn = Dict([("particles",n),("cutoff",0.0),("if_check_fluxes",false),("max_occ",1),("expander_fraction",100),("flux_direction","synth"),("layers",layers),("mdim",200),("if_save_data",false),("if_find_data",false),("filling",0.5),("onsite_strength",intstren),("lr","all"),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_results = run_synth_dims_generic(pdict_ttn)
    psi = all_results[1]=#

    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)
    f_ttn = all_files_ttn[1]
    d_ttn,m_ttn = read_data(dataloc_ttn * "/" * f_ttn; output_level=0)
    
    fourpt_mpo_1 = m_ttn["fourpt_momentum"]
    #fourpt_mpo_2 = m_ttn["fourpt_momentum_1"]
    #fourpt_mpo_mix = 0.5 .* (fourpt_mpo_1 + fourpt_mpo_2)

    #plot_four_point(fourpt_mpo_mix; plot_title="TTN GS1+GS2 $(lx)x$ly N=$n ULR=$intstren")
    plot_four_point(fourpt_mpo_1; plot_title="TTN GS1 $(lx)x$ly N=$n ULR=$intstren")
    #plot_four_point(fourpt_mpo_2; plot_title="TTN GS2 $(lx)x$ly N=$n ULR=$intstren")

    #=dataloc_ed = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren)])
    all_files_ed = find_data_file(pdict_ed,"ed",dataloc_ed; file_type="jld2")
    display(all_files_ed)
    f_ed = all_files_ed[1]
    d_ed,m_ed = read_data(joinpath(dataloc_ed,f_ed))

    fourpt_ed_1 = m_ed["fourpt_momentum"]

    plot_four_point(fourpt_ed_1; plot_title="ED GS1 $(lx)x$ly N=$n ULR=$intstren")=#

    #plot_four_point(abs.(fourpt_mpo_1 - fourpt_ed_1); plot_title="Diff TTN/ED $(lx)x$ly N=$n ULR=$intstren")

end=#

#= do compare both pts momentum MPO ED
if false
    lx,ly,n = 4,4,2
    layers = Int(log(2,lx*ly))
    intstren = 300.0

    
    #=dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    filter!(x -> !occursin("if_synth_rectangle",x),all_files)
    display(all_files)
    f = all_files[1]
    d,m = read_data(dataloc * "/" * f; output_level=0)
    psi = [d["ttn"],d["ttn_1"]]=#



    #

    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_check_fluxes",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",5),("if_find_data",true),("if_save_data",false)])
    #states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)

    println("Energy gap is $(nrgs[2] - nrgs[1])")
    
    pdict = Dict([("particles",n),("es_count",1),("expander_fraction",100),("if_check_fluxes",false),("layers",layers),("mdim",200),("if_save_data",false),("if_find_data",false),("filling",0.5),("onsite_strength",intstren),("lr","all"),("if_periodic_phys",true),("if_periodic_synth",true)])
    #all_results = run_synth_dims_generic(pdict)
    psi = all_results[1]



    #psi_fock = focking_vector(psis[1],lattice_params["full_basis"]; output_level=0)
    #psi_fock1 = focking_vector(psis[2],lattice_params["full_basis"]; output_level=0)

    #
    

    #=twopt_mpo = two_point(psi; output_level=0,if_plot=true,plot_title="MPO")
    twopt_fock = two_point(psi_fock,lattice_params; output_level=0,if_plot=true,plot_title="Fock")

    println("MPO version")
    display(twopt_mpo)
    println("Fock version")
    display(twopt_fock)=#

    m1 = [0.0,1/ly]
    m2 = [0.0,3/ly]

    #=overlapmat = zeros(Float64,3,2)
    overlapmat[1,1] = abs2(adjoint(psi_fock) * states[1])
    overlapmat[2,1] = abs2(adjoint(psi_fock) * states[2])
    overlapmat[3,1] = abs2(adjoint(psi_fock) * states[3])
    overlapmat[1,2] = abs2(adjoint(psi_fock1) * states[1])
    overlapmat[2,2] = abs2(adjoint(psi_fock1) * states[2])
    overlapmat[3,2] = abs2(adjoint(psi_fock1) * states[3])
    println("Overlap between TTN Fock and ED")
    display(overlapmat)=#


    #=fourpt_mpo = focking_matrix(four_point_mpowrapped,psi[1],lattice_params["full_basis"]; output_level=1,momentum1 = m1, momentum2 = m2)
    fourpt_ed = ft_fourpt_matrix(states[1],m1,m2,lattice_params; output_level=0)
    println("Does MPO match with ED: ",round.(fourpt_mpo,digits=10) == round.(fourpt_ed,digits=10))
    println("Does Transpose Match: ",round.(transpose(fourpt_mpo),digits=10) == round.(fourpt_ed,digits=10))=#

    # ttn1 with ortho_center at top has arrows pointing all away from oc
    # this could be the dagger already actually
    # this is where the transpose is
    # need to get single mpo to match primes and flows with MPO expectation part

    #thismpo = single_point_mpo(psi,"Adag"; coeff_kwargs=(Lx=lx,Ly=ly,), mapping=collect(1:lx*ly))

    fourpt_mpo_value = four_point(psi,m1,m2; output_level=0)
    fourpt_ed_value = ft_fourpt(states[1:2],m1,m2,lattice_params; output_level=0)
    println("The measured values are MPO=$(fourpt_mpo_value) and ED=$(fourpt_ed_value)")

    #=fourpt_mpo_value = four_point(psis; output_level=0, if_plot=true, plot_title="MPO $(lx)x$(ly) N=$n ULR=$intstren")
    fourpt_ed_value = four_point(states[1:2],lattice_params; output_level=0, if_plot=true, plot_title="ED $(lx)x$(ly) N=$n ULR=$intstren")
    println("MPO values are ")
    display(fourpt_mpo_value)
    println("ED values are ")
    display(fourpt_ed_value)=#

    #=expval = zeros(ComplexF64,2,2)
    expval[1,1] = adjoint(psi_fock) * fourpt_ed * psi_fock
    expval[2,1] = adjoint(psi_fock) * fourpt_ed * psi_fock1
    expval[1,2] = adjoint(psi_fock1) * fourpt_ed * psi_fock
    expval[2,2] = adjoint(psi_fock1) * fourpt_ed * psi_fock1
    println("The calculated value with ED is $(abs.(eigvals(expval)))")

    expval_mpo = zeros(ComplexF64,2,2)
    expval_mpo[1,1] = adjoint(psi_fock) * fourpt_mpo * psi_fock
    expval_mpo[2,1] = adjoint(psi_fock) * fourpt_mpo * psi_fock1
    expval_mpo[1,2] = adjoint(psi_fock1) * fourpt_mpo * psi_fock
    expval_mpo[2,2] = adjoint(psi_fock1) * fourpt_mpo * psi_fock1
    println("The calculated value with Fock MPO is $(abs.(eigvals(expval_mpo)))")
    
    expval_mpot = zeros(ComplexF64,2,2)
    expval_mpot[1,1] = adjoint(psi_fock) * transpose(fourpt_mpo) * psi_fock
    expval_mpot[2,1] = adjoint(psi_fock) * transpose(fourpt_mpo) * psi_fock1
    expval_mpot[1,2] = adjoint(psi_fock1) * transpose(fourpt_mpo) * psi_fock
    expval_mpot[2,2] = adjoint(psi_fock1) * transpose(fourpt_mpo) * psi_fock1
    println("The calculated value with Transpose Fock MPO is $(abs.(eigvals(expval_mpot)))")=#


end=#

#= fix transpose MPO
if false
    lx,ly,n = 4,4,2
    layers = Int(log(2,lx*ly))
    intstren = 0.0

    pdict = Dict([("particles",n),("expander_fraction",100),("if_check_fluxes",false),("layers",layers),("mdim",200),("if_save_data",false),("if_find_data",false),("filling",0.5),("onsite_strength",intstren),("lr","all"),("if_periodic_phys",true),("if_periodic_synth",true)])
    #all_results = run_synth_dims_generic(pdict)
    psi = all_results[1]

    mapss = zigzag_curve(lx,ly)
    coeff_kwargs = (Lx=lx,Ly=ly,)

    m1 = [0.0,1/ly]
    m2 = [0.0,3/ly]

    fourpt_val = four_point(psi)
    #println("4pt MPO value is $(fourpt_val)")
    display(fourpt_val)

end=#

#= do multi-state 4pt momentum MPO
if false
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))
    intstren = 0.0

    dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    filter!(x -> !occursin("if_synth_rectangle",x),all_files)
    display(all_files)
    f = all_files[1]
    d,m = read_data(dataloc * "/" * f; output_level=0)
    psi = d["ttn"]

    

end=#

#= test 4pt momentum with ED
if false
    lx,ly,n = 10,5,5
    intstren = 0.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("dataloc",dataloc),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",false)])
    states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)

    fourpt_vals = four_point(states[1],lattice_params; plot_title="ED $(lx)x$(ly) N=$n ULR=$intstren",if_plot=true)
    datadict = Dict([("fourpt_momentum",fourpt_vals)])
    modify_data(datadict,filepath,"metadata")
end=#

#= read 4pt data ED
if false
    println("Starting ED")
    lx,ly,n = 8,4,4
    intstren = 0.0
    
    #=dataloc_ed = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren)])
    all_files = find_data_file(pdict,"ed",dataloc_ed; file_type="jld2")
    display(all_files)
    f = all_files[1]
    d_ed,m_ed = read_data(joinpath(dataloc_ed,f))
    lattice_params = get_lattice_params_from_metadata(m_ed)
    states = d_ed["state"]=#

    m1 = [0.0,7/ly]
    m2 = [0.0,7/ly]

    #pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
    #states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict_ed; output_level=1)

    fourpt_val_ed = ft_fourpt_alberto(states[1],m1,m2,lattice_params)
    println("Final Value is $fourpt_val_ed")
    #fourpt_vals_ed = four_point(states[1],lattice_params; if_plot=false)
end=#

#= test 2pt momentum with ED
if false
    println("Starting ED")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",false)])
    states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)

    ks = [n/ly for n in 0:lx-1]
    ed_twopt_vals_1 = zeros(Float64,length(ks))
    ed_twopt_vals_2 = zeros(Float64,length(ks))
    for (idx,ky) in enumerate(ks)
        ed_twopt_vals_1[idx] = abs(ft_twopt(states[1],[0.0,ky],[0.0,ky],lattice_params; output_level=1))
        ed_twopt_vals_2[idx] = abs(ft_twopt(states[2],[0.0,ky],[0.0,ky],lattice_params; output_level=1))
    end
    scatter(ks .* ly,ed_twopt_vals_1,c="r",label="GS1")
    scatter(ks .* ly,ed_twopt_vals_2,c="b",label="GS2")
    legend()

    datadict = Dict([("mom_occs",ed_twopt_vals_1),("mom_occs_1",ed_twopt_vals_2)])
    modify_data(datadict,filepath,"metadata")

    #ed_twopt_vals = ft_twopt_alberto(states[1:2],mp,mp,lattice_params)

    #rho_ed1 = density_matrix(states[1],lattice_params)
    #rho_ed2 = density_matrix(states[2],lattice_params)
end=#

#= compare 4pt operator on 8x4 direct with Alberto
if false
    lx,ly,n = 8,4,4
    intstren = 0.0

    #=params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_check_fluxes",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",true)])
    lattice_params,hamilt_params,running_args = get_normal_model_params_ed(params_dict)    
    basis_dataloc = running_args.basis_dataloc
    full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
    lattice_params["full_basis"] = full_basis=#

    fourpt_table = zeros(ComplexF64,lx*ly,lx*ly,lx*ly,lx*ly)
    bf = open("fourPointTable.dat", "r")

    for (idx,line) in enumerate(eachline(bf))

        linevect = split(line,"\t")

        coord1 = [parse(Int,linevect[1]) + 1,parse(Int,linevect[2]) + 1]
        coord2 = [parse(Int,linevect[3]) + 1,parse(Int,linevect[4]) + 1]
        coord3 = [parse(Int,linevect[5]) + 1,parse(Int,linevect[6]) + 1]
        coord4 = [parse(Int,linevect[7]) + 1,parse(Int,linevect[8]) + 1]
        lin1 = linear_index(coord1,lx,ly)
        lin2 = linear_index(coord2,lx,ly)
        lin3 = linear_index(coord3,lx,ly)
        lin4 = linear_index(coord4,lx,ly)

        fourpt_val = parse(Float64,linevect[9]) + im*parse(Float64,linevect[10])

        fourpt_table[lin1,lin2,lin3,lin4] = fourpt_val
    end
    close(bf)

    #pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",false)])
    #states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)


    #=coord1 = [7,4]
    coord2 = [1,1]
    coord3 = [6,4]
    coord4 = [1,3]
    lin1 = linear_index(coord1,lx,ly)
    lin2 = linear_index(coord2,lx,ly)
    lin3 = linear_index(coord3,lx,ly)
    lin4 = linear_index(coord4,lx,ly)=#

    mismatch_sites = []
    for lin1 in 1:lx*ly
        for lin2 in 1:lx*ly
            println("Working on lin1 = $lin1, lin2 = $lin2")
            for lin3 in 1:lx*ly
                for lin4 in 1:lx*ly
                    bigop = four_point_operator(lin1,lin2,lin3,lin4,lattice_params)
                    expval = adjoint(states[1]) * bigop * states[1]
                    if !isapprox(expval,fourpt_table[lin1,lin2,lin3,lin4],atol=1e-3)
                        println("Mismatch at ",coordinate(lin1,lx,ly) .- [1,1]," ",coordinate(lin2,lx,ly) .- [1,1]," ",coordinate(lin3,lx,ly) .- [1,1]," ",coordinate(lin4,lx,ly) .- [1,1])
                        println("Expected ",fourpt_table[lin1,lin2,lin3,lin4])
                        println("Got ",expval)
                        push!(mismatch_sites,[[lin1,lin2,lin3,lin4]])
                    end
                end
            end
        end
    end
    println("All Matching! Woohoo")
    
    
end=#

#= test 4pt ED with Alberto's real space data
if false

    lx,ly,n = 8,4,4
    intstren = 0.0

    #=fourpt_table = zeros(ComplexF64,lx*ly,lx*ly,lx*ly,lx*ly)
    bf = open("fourPointTable.dat", "r")

    for (idx,line) in enumerate(eachline(bf))

        linevect = split(line,"\t")

        coord1 = [parse(Int,linevect[1]) + 1,parse(Int,linevect[2]) + 1]
        coord2 = [parse(Int,linevect[3]) + 1,parse(Int,linevect[4]) + 1]
        coord3 = [parse(Int,linevect[5]) + 1,parse(Int,linevect[6]) + 1]
        coord4 = [parse(Int,linevect[7]) + 1,parse(Int,linevect[8]) + 1]
        lin1 = linear_index(coord1,lx,ly)
        lin2 = linear_index(coord2,lx,ly)
        lin3 = linear_index(coord3,lx,ly)
        lin4 = linear_index(coord4,lx,ly)

        fourpt_val = parse(Float64,linevect[9]) + im*parse(Float64,linevect[10])

        fourpt_table[lin1,lin2,lin3,lin4] = fourpt_val
    end
    close(bf)

    dataloc_ed = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren)])
    all_files = find_data_file(pdict,"ed",dataloc_ed; file_type="jld2")
    display(all_files)
    f = all_files[1]
    d_ed,m_ed = read_data(joinpath(dataloc_ed,f))
    lattice_params = get_lattice_params_from_metadata(m_ed)
    states = d_ed["state"]=#

    m1 = [0.0,1/ly]
    m2 = [0.0,0/ly]

    fourpt_vals_alberto = ft_fourpt_alberto(states[1],m1,m2,lattice_params; if_readrealspace=true)
    #fourpt_vals_patrick = ft_fourpt_alberto(states[1],m1,m2,lattice_params; if_readrealspace=false)

    #println("Patrick's Value is $fourpt_vals_patrick")
    println("Alberto's Value is $fourpt_vals_alberto")
end=#

#= compare momoccs and 4pt ED and TTN 8x4 
if false
    lx,ly,n = 8,4,4
    intstren = 300.0
    layers = Int(log(2,lx*ly))

    dataloc_ed = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren)])
    all_files_ed = find_data_file(pdict_ed,"ed",dataloc_ed; file_type="jld2")
    display(all_files_ed)
    f_ed = all_files_ed[1]
    d_ed,m_ed = read_data(joinpath(dataloc_ed,f_ed))
    
    #=momoccs_ed = m_ed["mom_occs"]
    plot(collect(0:lx-1),momoccs_ed,"-p",c="b",label="ED")
    momoccs_ed_1 = m_ed["mom_occs_1"]
    plot(collect(0:lx-1),momoccs_ed_1,"-p",c="r",label="ED GS2")
    xlabel("Momentum")
    ylabel("Occupancy")
    title("Momentum Occupancies for $(lx)x$(ly) N=$n ULR=$intstren")=#
    fourpt_ed = m_ed["fourpt_momentum"]
    #plot_four_point(fourpt_ed; plot_title="ED $(lx)x$(ly) N=$n ULR=$intstren")
    #plot_four_point(fourpt_ed[:,3],2; plot_title="ED $(lx)x$(ly) N=$n ULR=$intstren")


    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)
    f_ttn = all_files_ttn[1]
    d_ttn,m_ttn = read_data(dataloc_ttn * "/" * f_ttn; output_level=0)

    #momoccs_ttn = m_ttn["mom_occs"]
    #plot(collect(0:lx-1),momoccs_ttn,"-p",c="r",label="TTN GS1")
    #momoccs_ttn_1 = m_ttn["mom_occs_1"]
    #plot(collect(0:lx-1),momoccs_ttn_1,"-p",c="r"; label="TTN GS1")
    #legend()

    fourpt_ttn = m_ttn["fourpt_momentum"]
    #plot_four_point(fourpt_ttn; plot_title="TTN GS1 $(lx)x$(ly) N=$n ULR=$intstren")
    #plot_four_point(fourpt_ttn[:,4],3; plot_title="TTN GS1 $(lx)x$(ly) N=$n ULR=$intstren")
    fourpt_ttn_1 = m_ttn["fourpt_momentum_1"]
    #plot_four_point(fourpt_ttn_1; plot_title="TTN GS1 $(lx)x$(ly) N=$n ULR=$intstren")
    #plot_four_point(fourpt_ttn_1[:,3],2; plot_title="TTN GS2 $(lx)x$(ly) N=$n ULR=$intstren")

    #=fourpt_gsmanifold_ttn = zeros(Float64,lx,lx)
    fourpt_gsmanifold_ttn .+= fourpt_ttn_1
    shifted_gs2 = zeros(Float64,lx,lx)
    for i in 1:lx-1
        for j in 1:lx-1
            shifted_gs2[i+1,j+1] = fourpt_ttn_1[i,j]
        end
    end
    fourpt_gsmanifold_ttn .+= shifted_gs2
    plot_four_point(fourpt_gsmanifold_ttn ./ 2; plot_title="TTN GS1+GS2 $(lx)x$(ly) N=$n ULR=$intstren")=#

end=#

#= save fock vector for TTN
if false
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))
    intstren = 0.0
    
    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)
    #=f_ttn = all_files_ttn[1]
    d_ttn,m_ttn = read_data(dataloc_ttn * "/wavefunc" * f_ttn; output_level=0)

    psi = d_ttn["ttn"]
    ttn_fockvector_1 = focking_vector(psi,lattice_params["full_basis"])
    datadict1 = Dict([("fock_vector",ttn_fockvector_1)])
    modify_data(datadict1,dataloc_ttn * "/" * f_ttn,"metadata")

    psi_1 = d_ttn["ttn_1"]
    ttn_fockvector_2 = focking_vector(psi_1,lattice_params["full_basis"])
    datadict2 = Dict([("fock_vector_1",ttn_fockvector_2)])
    modify_data(datadict2,dataloc_ttn * "/" * f_ttn,"metadata")=#



end=#

#= check overlap of TTN and ED states
if false
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))
    intstren = 300.0

    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)
    f_ttn = all_files_ttn[1]
    d_ttn,m_ttn = read_data(dataloc_ttn * "/" * f_ttn; output_level=0)

    f1_ttn = m_ttn["fock_vector"]
    f2_ttn = m_ttn["fock_vector_1"]

    dataloc_ed = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren)])
    all_files_ed = find_data_file(pdict_ed,"ed",dataloc_ed; file_type="jld2")
    display(all_files_ed)
    f_ed = all_files_ed[1]
    d_ed,m_ed = read_data(joinpath(dataloc_ed,f_ed))

    f1_ed = d_ed["state"][1]
    f2_ed = d_ed["state"][2]
    f3_ed = d_ed["state"][3]

    overlapmat = zeros(Float64,2,3)
    overlapmat[1,1] = abs2(adjoint(f1_ttn) * f1_ed)
    overlapmat[1,2] = abs2(adjoint(f1_ttn) * f2_ed)
    overlapmat[2,1] = abs2(adjoint(f2_ttn) * f1_ed)
    overlapmat[2,2] = abs2(adjoint(f2_ttn) * f2_ed)
    overlapmat[1,3] = abs2(adjoint(f1_ttn) * f3_ed)
    overlapmat[2,3] = abs2(adjoint(f2_ttn) * f3_ed)

    display(overlapmat)

end=#

#= 4pt momentum for TTN 8x4 from Fock vector
if false
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))
    intstren = 0.0

    dataloc_ttn = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict_ttn = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files_ttn = find_data_file(pdict_ttn,"ttn",dataloc_ttn)
    display(all_files_ttn)
    f_ttn = all_files_ttn[1]
    d_ttn,m_ttn = read_data(dataloc_ttn * "/" * f_ttn; output_level=0)

    fourpt_fock_1 = m_ttn["fourpt_momentum_fock"]
    fourpt_fock_2 = m_ttn["fourpt_momentum_fock_1"]

    psi1 = m_ttn["fock_vector"]
    psi2 = m_ttn["fock_vector_1"]

    m1 = [0.0,0/ly]
    m2 = [0.0,1/ly]

    lattice_params = get_lattice_params_from_ttn(m_ttn)

    fourpt_ttn_mixed1,fourpt_ttn_mixed2 = abs.(ft_fourpt([psi1,psi2],m1,m2,lattice_params; output_level=0))

    println("The unmixed Fock values are $fourpt_fock_1 and $fourpt_fock_2")
    println("From mixing we get $fourpt_ttn_mixed1 and $fourpt_ttn_mixed2")

    #plot_four_point(fourpt_fock_1; plot_title="TTN Fock GS1 $(lx)x$(ly) N=$n ULR=$intstren")
    #plot_four_point(fourpt_fock_2; plot_title="TTN Fock GS2 $(lx)x$(ly) N=$n ULR=$intstren")

    #=fourpt_ttn = m_ttn["fourpt_momentum"]
    plot_four_point(fourpt_ttn; plot_title="TTN GS1 $(lx)x$(ly) N=$n ULR=$intstren")

    dataloc_ed = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict_ed = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren)])
    all_files_ed = find_data_file(pdict_ed,"ed",dataloc_ed; file_type="jld2")
    display(all_files_ed)
    f_ed = all_files_ed[1]
    d_ed,m_ed = read_data(joinpath(dataloc_ed,f_ed))

    fourpt_ed = m_ed["fourpt_momentum"]
    plot_four_point(fourpt_ed; plot_title="ED $(lx)x$(ly) N=$n ULR=$intstren")=#
end=#


























"fin"
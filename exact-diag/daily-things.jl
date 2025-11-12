#####################################################
#=

This file contains any random functions written to do one-off tasks

Depends on:
    execute-ed.jl

=#
######################################################

include("execute-ed.jl")
include("plottings.jl")
include("../other-funcs/basic-2d-plottings.jl")

function datacollection_flatness(Lx::Int64,Ly::Int64,N::Int64; kwargs...)
    hanis::Float64 = get(kwargs,:hopping_anisotropy,1.0)
    if_hatsugai::Bool = get(kwargs,:if_hatsugai,true)

    intstrens = get(kwargs,:intstrens,nothing)
    if isnothing(intstrens)
        intstren_count::Int64 = get(kwargs,:intstren_count,11)
        intstren_start::Float64 = get(kwargs,:intstren_start,0.0)
        intstren_end::Float64 = get(kwargs,:intstren_end,2.0)
        intstrens = range(intstren_start,intstren_end,length=intstren_count)
    end

    tws_count::Int64 = get(kwargs,:tws_count,10)
    tws_start::Float64 = get(kwargs,:tws_start,0.0)
    tws_end::Float64 = get(kwargs,:tws_end,1.0)
    tws::Vector{Float64} = range(tws_start,tws_end,length=tws_count)

    println("Starting Flatness Data Collection for $(Lx)x$(Ly) N=$(N) from Intstrengths $(intstrens[1]) - $(intstrens[end]) and Twists $(tws_start) - $(tws_end)")
    sleep(1.0)

    nev = 10
    for (idx,intstren) in enumerate(intstrens)
        if if_hatsugai
            ref_multis,rm1,rm2 = get_reference_multiplets(Lx,Ly,N; interaction_strength=intstren,hopping_anisotropy=hanis,if_make_new=true)
        end
        for (idx2,tw1) in enumerate(tws)
            for (idx3,tw2) in enumerate(tws)
                params_dict = Dict([("output_level",1),("Lx",Lx),("Ly",Ly),("N",N),("tw1",tw1),("tw2",tw2),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",nev),("if_find_data",true),("if_save_data",true)])
                states,nrgs,rhos,filepath,if_found = run_normal_ed(params_dict; output_level=1)

                if if_hatsugai
                    if if_found
                        d,m = read_data_jld2(filepath; output_level=0)
                        if !haskey(m,"omega")
                            lambda1,lambda2,omega = get_hatsugaifull(states[1],states[2],ref_multis; if_save=true,filepath=filepath,ref_multis_filenames=[rm1,rm2])
                        end
                    else
                        lambda1,lambda2,omega = get_hatsugaifull(states[1],states[2],ref_multis; if_save=true,filepath=filepath,ref_multis_filenames=[rm1,rm2])
                    end
                end
            end
        end
    end
end


#= set ham in data file to nothing
if false
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    all_files = readdir(dataloc)
    filter!(x -> occursin("jld2",x),all_files)
    for f in all_files
        d,m = read_data(joinpath(dataloc,f); output_level=0)

        if !isnothing(m["H"])
            println("Modifying file: $f")
            modify_data(Dict([("H",nothing)]),joinpath(dataloc,f),"metadata"; output_level=0)
        end

    end
end=#

#= data collection for 4pt
if false
    lx,ly,n = 10,4,5
    intstren = 0.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true)])
    
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    display(all_files)

    for f in all_files
        d,m = read_data(dataloc * "/" * f; output_level=0)
        lattice_params = get_lattice_params_from_metadata(m)

        if haskey(m,"fourpt_momentum")
            println("Already has 4pt data")
            continue
        else
            fourpt_vals = four_point(d["state"][1],lattice_params)
            datadict = Dict([("fourpt_momentum",fourpt_vals)])
            modify_data(datadict,joinpath(dataloc,f),"metadata")
        end
        #two_fourpts = [ft_fourpt(d["state"][1],[0.0,n/ly],[0.0,n/ly],lattice_params) for n in 0:1]
        #datadict = Dict([("fourpt_momentum_diag",two_fourpts)])
        #modify_data(datadict,filepath,"metadata"; output_level=0)
    end
end=#

#= test 4pt diag max as order parameter for 6x3
if false
    for n in [3,4]
        lx,ly,n = Int(2*n),n,n
        intstrens = range(0.0,10.0,length=21)

        orderparams = zeros(Float64,length(intstrens))
        for (idx,intstren) in enumerate(intstrens)
            params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("nev",10),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("if_find_data",true),("if_save_data",true)])
            states,nrgs,rhos,filepath,if_found,latpara,hamiltpara = run_normal_ed(params_dict; output_level=1)

            fourpt_vals = four_point_diag(states[1],latpara; if_plot=false)

            orderparams[idx] = maximum(fourpt_vals)
        end
        plot(intstrens,orderparams,"-p",label="N=$n")
        xlabel("Interaction Strength")
        ylabel("Max 4pt Diagonal")
        title("k-DW Order Parameter")
    end

end=#

#= vary potential strength to look at occupancy visibility
if true
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/periodic-potential")
    intstren = 300.0
    lx,ly,n = 8,4,4
    anises = [1e-4,1e-3,1e-2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    cols = get_colors(length(anises))
    potstrens = round.(vcat(1e-3,1e-4,10 .^ range(-2,2,length=10)),digits=4)
    for (idx2,anis) in enumerate(anises)
        for (idx,potstren) in enumerate(potstrens)
            pdict_pp = Dict([("output_level",1),("dataloc",dataloc),("periodic_potential_strength",potstren),("tx",anis),("ty",1.0),("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",20),("if_find_data",true),("if_save_data",false)])
            states_pp,nrgs_pp,rhos_pp,filepath_pp,if_found_pp,lattice_params_pp,hamilt_params_pp = run_normal_ed(pdict_pp; output_level=1)

            #=occs = get_occupancy(states_pp[1],lattice_params_pp; if_plot = false)
            visibs[idx] = minimum(occs) / maximum(occs)
            
            if idx == 1
                scatter(potstren,visibs[idx],c=cols[idx2],label="tx= $anis")
            else
                scatter(potstren,visibs[idx],c=cols[idx2])
            end
            xlabel("Periodic Potential Strenght "*L"V_0")
            ylabel("Occs Flatness, Min(Occs) / Max(Occs)")
            xscale("log")
            title("Flatness vs Periodic Potential Strength, ULR = $intstren, $(lx)x$(ly) N=$n ")=#
        
        end
    end
end=#

#= little plot for 4x4 n=2 adiabatic condition tx and periodic potential strength
if false
    lx,ly,n = 4,4,2
    intstren = 300.0
    potstrens = vcat([1e-4,1e-3],round.(10 .^ range(-2,2,length=10),digits=4),100.1)
    anises = [1e-4,1e-3,1e-2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
    ftxs_max = zeros(Float64,length(anises)-1,length(potstrens)-1)
    fpps_max = zeros(Float64,length(anises)-1,length(potstrens)-1)
    ftxs_first = zeros(Float64,length(anises)-1,length(potstrens)-1)
    fpps_first = zeros(Float64,length(anises)-1,length(potstrens)-1)
    for (i2,v0) in enumerate(potstrens[1:end-1])
        #fig = figure()
        for (i1,tx) in enumerate(anises[1:end-1])
            params_dict = Dict([("output_level",1),("periodic_potential_strength",v0),("tx",tx),("ty",1.0),("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",40),("if_find_data",true),("if_save_data",false)])
            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)

            #plot_spectrum(anises,nrgs,i1,params_dict["nev"],"Physical Hopping",true; plot_title=" $(lx)x$(ly) V0=$v0")

            nrg_gaps = nrgs[2:end] .- nrgs[1]

            ftxs_0,ogtxs = tx_adiabatic_condition(states[1],states[2:end],nrg_gaps,lattice_params,hamilt_params)
            fpps_0,ogpps = periodic_potential_adiabatic_condition(states[1],states[2:end],nrg_gaps,lattice_params,hamilt_params)

            ftxs_max[i1,i2] = maximum(ftxs_0)
            fpps_max[i1,i2] = maximum(fpps_0)

            ftxs_first[i1,i2] = ftxs_0[2]
            fpps_first[i1,i2] = fpps_0[2]
        end
    end

    plot_adiabatic_condition(potstrens,anises,ftxs_max,fpps_max; plot_title="Max Value $(lx)x$(ly)")
    xlabel("Periodic Potential Strength, "*L"V_0")
    ylabel("Physical Hopping, "*L"t_x")
    xscale("log")

    plot_adiabatic_condition(potstrens,anises,ftxs_first,fpps_first; plot_title="First Excited State $(lx)x$(ly)")
    xlabel("Periodic Potential Strength, "*L"V_0")
    ylabel("Physical Hopping, "*L"t_x")
    xscale("log")


end=#

#= save adiabatic condition data for periodic potential
if false
    lx,ly,n = 8,4,4
    intstren = 300.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/periodic-potential")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    potstrens = vcat([1e-4,1e-3],round.(10 .^ range(-2,2,length=10),digits=4),100.1)
    anises = [1e-4,1e-3,1e-2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
    ftxs = zeros(Float64,length(anises)-1,length(potstrens)-1)
    fpps = zeros(Float64,length(anises)-1,length(potstrens)-1)
    for f in all_files
        d,m = read_data(joinpath(dataloc,f); output_level=0)

        lattice_params,hamilt_params = make_latticehamilt_params_from_metadata(m)
        
        #idx1 = findfirst(x -> x == m["tx"],anises)
        #idx2 = findfirst(x -> x == round(m["periodic_potential_strength"],digits=4),potstrens)

        #println("Found file with tx=$(m["tx"]), potstren=$(m["periodic_potential_strength"]) at indices idx1=$(idx1), idx2=$(idx2)")

        nrg_gaps = d["nrg"][2:end] .- d["nrg"][1]

        ftxs_0,ogtxs = tx_adiabatic_condition(d["state"][1],d["state"][2:end],nrg_gaps,lattice_params,hamilt_params)
        fpps_0,ogpps = periodic_potential_adiabatic_condition(d["state"][1],d["state"][2:end],nrg_gaps,lattice_params,hamilt_params)

        datadict = Dict([("ftx_0",ftxs_0),("fpp_0",fpps_0)])
        modify_data(datadict,joinpath(dataloc,f),"metadata"; output_level=0)

        #ftxs[idx1,idx2] = magval_tx
        #fpps[idx1,idx2] = magval_pp

    end
    #plot_adiabatic_condition(potstrens,anises,ftxs,fpps; plot_title="$(lx)x$(ly) N=$(n) ULR=$intstren")
    #xlabel("Periodic Potential Strength, "*L"V_0")
    #ylabel("Physical Hopping, "*L"t_x")
    #xscale("log")
end=#

#= plot adiabatic condition for periodic potential
if false
    lx,ly,n = 8,4,4
    intstren = 300.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/periodic-potential")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    potstrens = vcat([1e-4,1e-3],round.(10 .^ range(-2,2,length=10),digits=4),100.1)
    anises = [1e-4,1e-3,1e-2,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1]
    ftxs = zeros(Float64,length(anises)-1,length(potstrens)-1)
    fpps = zeros(Float64,length(anises)-1,length(potstrens)-1)
    for f in all_files
        d,m = read_data(joinpath(dataloc,f); output_level=0)
        
        idx1 = findfirst(x -> x == m["tx"],anises)
        idx2 = findfirst(x -> x == round(m["periodic_potential_strength"],digits=4),potstrens)

        println("Found file with tx=$(m["tx"]), potstren=$(m["periodic_potential_strength"]) at indices idx1=$(idx1), idx2=$(idx2)")

        ftx_val,ftx_ind = findmax(m["ftx_0"])
        fpp_val,fpp_ind = findmax(m["fpp_0"])

        display(m["ftx_0"])
        display(m["fpp_0"])

        ftxs[idx1,idx2] = ftx_val
        fpps[idx1,idx2] = fpp_val

    end
    plot_adiabatic_condition(potstrens,anises,ftxs,fpps; plot_title="$(lx)x$(ly) N=$(n) ULR=$intstren")
    xlabel("Periodic Potential Strength, "*L"V_0")
    ylabel("Physical Hopping, "*L"t_x")
    xscale("log")
end=#

#= look at higher ACs for dominant parameter(tx/ULR) and make paper plot
if false
    lx,ly,n = 8,4,4
    hanis = 1.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    display(all_files)

    f_ulrs = []
    allgaps = []
    ulrs = []

    for f in all_files

        fileparams = get_params_dict_from_filename(f)
        fileparams["interaction_strength"] > 10.0 && continue

        d,m = read_data(joinpath(dataloc,f); output_level=0)

        lattice_params,hamilt_params = make_latticehamilt_params_from_metadata(m)
        #=ops = Dict([("nev",50),("if_save_data",true),("if_find_data",true)])
        normparas = get_normal_params_from_lattham(lattice_params,hamilt_params,ops)
        states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(normparas; output_level=1)=#

        nrg_gaps = d["nrg"][3:end] .- d["nrg"][2]
        states = d["state"]

        #fuirs_1,ogfuirs1 = uir_adiabatic_condition(states[2],states[3:end],nrg_gaps,lattice_params,hamilt_params,nothing)
        fuirs_1,ogfuirs1 = uir_adiabatic_condition(states[1],states[1],0.0,lattice_params,hamilt_params,nothing)
        fuirs_2,ogfuirs2 = uir_adiabatic_condition(states[2],states[2],0.0,lattice_params,hamilt_params,nothing)

        println("OG F_UIR 1: $ogfuirs1")
        println("OG F_UIR 2: $ogfuirs2")

        #fuirs_0 = m["fuir_0"]

        #fuirs_avg = [0.5 .* (fuirs_0[i+1][1] .+ fuirs_1[i][1]) for i in 1:length(fuirs_1)]

        #datadict = Dict([("fuir_1",fuirs_1)])
        #modify_data(datadict,joinpath(dataloc,f),"metadata"; output_level=0)

        append!(allgaps,[[nrg_gaps]])
        append!(f_ulrs,[[fuirs_avg]])
        append!(ulrs,[m["U"][end]])

    end

    #=clog = [log10.(f_ulrs[i][1]) for i in 1:length(f_ulrs)]

    cols = ["#82AC9F","#C73E1D","#36213E"]

    target = cols[3]
    whitetohex = matplotlib.colors.LinearSegmentedColormap.from_list(
        "white_to_hex", ["#ffffff", target]
    )

    for i in 1:length(ulrs)
        xs = ulrs[i] .* ones(length(allgaps[i][1]))
        ys = allgaps[i][1]
        scatter(xs,ys,c=clog[i],cmap=whitetohex)
    end
    colorbar().set_label(L"log_{10} (F_{U_{IR}})")
    xlabel("Interaction Strength, "*L"U_{LR}")
    ylabel("Energy Gap")=#

    

end=#

#= look at higher adiabatic conditions
if false
    lx,ly,n = 8,4,4
    hanis = 1.0
    #ppstren = 1e-2
    intstren = 300.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")#/periodic-potential")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    display(all_files)
    txs = []
    pps = []
    ftxs = []
    fpps = []
    max_ftxs = []
    max_fpps = []
    for f in all_files
        d,m = read_data(joinpath(dataloc,f); output_level=0)

        #=lattice_params,hamilt_params = make_latticehamilt_params_from_metadata(m)
        ops = Dict([("nev",50),("if_save_data",false),("if_find_data",false)])
        normparas = get_normal_params_from_lattham(lattice_params,hamilt_params,ops)
        states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(normparas; output_level=1)

        nrg_gaps = nrgs[2:end] .- nrgs[1]
        
        ftxs_0,ogtxs = tx_adiabatic_condition(states[1],states[2:end],nrg_gaps,lattice_params,hamilt_params)
        fpps_0,ogpps = periodic_potential_adiabatic_condition(states[1],states[2:end],nrg_gaps,lattice_params,hamilt_params)
        #fuirs_0,ogfuirs = uir_adiabatic_condition(states[1],states[2:end],nrg_gaps,lattice_params,hamilt_params,nothing)

        mags = sqrt.(ftxs_0 .^2 .+ fpps_0 .^ 2)=#

        if haskey(m,"fpp_0") && haskey(m,"ftx_0")
            append!(txs,[m["tx"]])
            append!(pps,[m["periodic_potential_strength"]])
            max_ftx,maxidx_ftx = findmax(m["ftx_0"])
            max_fpp,maxidx_fpp = findmax(m["fpp_0"])
            append!(ftxs,[max_ftx])
            append!(fpps,[max_fpp])
            append!(max_ftxs,[maxidx_ftx])
            append!(max_fpps,[maxidx_fpp])
        else
            continue
        end
        
        #=scatter(1:length(ftxs_0),ftxs_0,label=L"V_0="*"$(round(m["periodic_potential_strength"],digits=4))")
        scatter(1:length(mags),mags,label=L"V_0="*"$(round(m["periodic_potential_strength"],digits=4))")
        xlabel("Excited State Index")
        ylabel("Adiabatic Condition, "*L"F_{tx}, F_{pp}")
        title("AC for Higher States: $(lx)x$(ly) ULR=$intstren, tx=$hanis")
        legend()
        yscale("log")=#
        
    end

    plotting_txs = sort(unique(txs))
    plotting_pps = sort(unique(pps))
    plotting_ftxs = zeros(Float64,length(plotting_txs),length(plotting_pps))
    plotting_fpps = zeros(Float64,length(plotting_txs),length(plotting_pps))
    plotting_maxftxs = zeros(Float64,length(plotting_txs),length(plotting_pps))
    plotting_maxfpps = zeros(Float64,length(plotting_txs),length(plotting_pps))
    plotting_txs = vcat(plotting_txs,[plotting_txs[end] + 0.01])
    plotting_pps = vcat(plotting_pps,[plotting_pps[end] + 0.01])
    for i in 1:length(mags)
        tx = findfirst(x -> x == txs[i],plotting_txs)
        pp = findfirst(x -> x == pps[i],plotting_pps)
        plotting_ftxs[tx,pp] = ftxs[i]
        plotting_fpps[tx,pp] = fpps[i]
        plotting_maxftxs[tx,pp] = max_ftxs[i]
        plotting_maxfpps[tx,pp] = max_fpps[i]
    end
    #=pcolormesh(plotting_pps,plotting_txs,log10.(plotting_mags))
    colorbar()
    ylabel("Hopping Anisotropy, "*L"t_x")
    xlabel("Periodic Potential Strength, "*L"V_0")=#
    #=plot_adiabatic_condition(plotting_pps,plotting_txs,plotting_fpps,plotting_ftxs; plot_title="$(lx)x$(ly) N=$(n) ULR=$intstren")
    xlabel("Periodic Potential Strength, "*L"V_0")
    ylabel("Hopping Anisotropy, "*L"t_x")
    xscale("log")=#
    
    fig = figure()
    pcolormesh(plotting_pps,plotting_txs,plotting_maxftxs)
    xlabel("Periodic Potential Strength, "*L"V_0")
    ylabel("Hopping Anisotropy, "*L"t_x")
    xscale("log")
    title("Most Relevant Excited State F_tx, $(lx)x$(ly) N=$(n) ULR=$intstren")
    colorbar()

    fig = figure()
    pcolormesh(plotting_pps,plotting_txs,plotting_maxfpps)
    xlabel("Periodic Potential Strength, "*L"V_0")
    ylabel("Hopping Anisotropy, "*L"t_x")
    xscale("log")
    title("Most Relevant Excited State F_pp, $(lx)x$(ly) N=$(n) ULR=$intstren")
    colorbar()



end=#

#= real/momentum space 4pt for all ulrs and txs
if false
    lx,ly,n = 8,4,4
    hanis = 1.0
    intstren = 300.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren),("hopping_anisotropy",hanis)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    display(all_files)
    
    d,m = read_data(joinpath(dataloc,all_files[1]); output_level=0)
    latparas = get_lattice_params_from_metadata(m)
    
    r4pt = realspace_fourpoint_full(d["state"][1],latparas; output_level=1)
    fig = figure()
    imshow(real.(sum(r4pt, dims=[3,4])[:,:,1,1]); origin="lower",vmin=0.0,vmax=0.5)
    xlabel("x")
    ylabel("x'")
    colorbar()
    title("Real Space 4pt for $(lx)x$(ly) N=$(n) tx=$hanis ULR=$intstren")

    m4pt = four_point(d["state"][1],latparas; if_plot=false)
    fig = figure()
    imshow(m4pt; origin="lower",vmin=0.0,vmax=0.5)
    xlabel("k_x")
    ylabel("k_x'")
    colorbar()
    title("Momentum Space 4pt for $(lx)x$(ly) N=$(n) tx=$hanis ULR=$intstren")
end=#


#= track overlap with periodic potential
if true
    anises = [1e-4,1e-3,1e-2,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
    intstren = 300.0
    lx,ly,n = 8,4,4
    all_overlaps1 = Matrix{Float64}(undef,20,length(anises))
    all_overlaps2 = Matrix{Float64}(undef,20,length(anises))
    for (idx,anis) in enumerate(anises)
        #pdict_found = Dict([("output_level",1),("tx",anis),("ty",1.0),("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",false)])
        #states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict_found; output_level=1)

        #plot_spectrum(anises,nrgs,idx,pdict_found["nev"],"Physical Hopping",true; plot_title=" $(lx)x$(ly) N=$n ULR=$intstren")

        pdict_pp = Dict([("output_level",1),("periodic_potential_strength",1e-2),("tx",anis),("ty",1.0),("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])
        states_pp,nrgs_pp,rhos_pp,filepath_pp,if_found_pp,lattice_params_pp,hamilt_params_pp = run_normal_ed(pdict_pp; output_level=1)

        plot_spectrum(anises,nrgs_pp,idx,pdict_pp["nev"],"Physical Hopping",true; plot_title=" with PP $(lx)x$(ly) N=$n ULR=$intstren")

        #for i in 1:20
        #    all_overlaps1[i,idx] = abs2(dot(states[1],states_pp[i]))
        #    all_overlaps2[i,idx] = abs2(dot(states[2],states_pp[i]))
        #end
        

        #=if idx == 1
            scatter(anis,overlaps["s11"][end],c="b",label="s11")
            scatter(anis,overlaps["s21"][end],c="g",label="s21")
            scatter(anis,overlaps["s12"][end],c="r",label="s12")
            scatter(anis,overlaps["s22"][end],c="m",label="s22")
            xlabel("Hopping Anisotropy, "*L"t_x")
            ylabel("Overlap with Periodic Potential")
            legend()
        else
            scatter(anis,overlaps["s11"][end],c="b")
            scatter(anis,overlaps["s21"][end],c="g")
            scatter(anis,overlaps["s12"][end],c="r")
            scatter(anis,overlaps["s22"][end],c="m")
        end=#
    end
    #=fig = figure()
    plot(anises,overlaps["s11"],"-o",label="State 11")
    plot(anises,overlaps["s21"],"-o",label="State 21")
    plot(anises,overlaps["s12"],"-o",label="State 12")
    plot(anises,overlaps["s22"],"-o",label="State 22")
    xlabel("Hopping Anisotropy, "*L"t_x")
    ylabel(L"\vert \langle \psi_{GSi} \vert \psi_{GSj}^{PP} \rangle \vert^2")
    title("Overlap with Periodic Potential GS1/2 for $(lx)x$(ly) N=$n ULR=$intstren")
    legend()=#
        
end=#

#= plot adiabatic condition for uir
if false
    lx,ly,n = 8,4,4
    hanis = 1.0
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis)])
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    #display(all_files)
    uirs = []
    intstrens = []
    for f in all_files
        d,m = read_data(joinpath(dataloc,f); output_level=0)

        if haskey(m,"fuir_12")
            fuir = m["fuir_12"]
            append!(uirs,[fuir])
            intstren = m["U"][end]
            append!(intstrens,[intstren])
        else
            continue
        end
        scatter(intstren,fuir,c="b")
        xlabel("Interaction Strength, "*L"U_{IR}")
        ylabel(L"F_{uir}")
        title("Adiabatic Condition for UIR, $(lx)x$(ly) N=$(n)")
    end
end=#

#= plot adiabatic condition
if false
    ly = 4
    for lx in [8,10]
    #lx,ly,n = 8,4,4
    n = Int(lx/2)
    txs = []
    intstren = 300.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    adiabs = []
    for f in all_files
        d,m = read_data(joinpath(dataloc,f); output_level=0)
        
        if haskey(m,"fuir_12") && haskey(m,"ftx_12")
            append!(txs,[m["tx"]])
            fuir = m["fuir_12"]
            ftx = m["ftx_12"]
            append!(adiabs,[sqrt(fuir^2 + ftx^2)])
        else
            continue
        end
    end
    #fig = figure()
    scatter(txs,adiabs,label="Lx=$lx")
    xlabel("Hopping Anisotropy, "*L"t_x")
    ylabel(L"\sqrt{F_{uir}^2 + F_{tx}^2}")
    title("Adiabatic Condition for ULR=$intstren")
    yscale("log")
    legend()
    end

end=#

#= testing for adiabatic condition
if false
    lx,ly,n = 8,4,4
    txs = [1e-3,1e-2,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0]
    #intstren = 300.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    intstrens = [300]#[0.0,10.0,50.0,100.0,150.0,200.0,250.0,300.0,350.0,400.0,450.0,500.0]
    fuirs = zeros(Float64,length(intstrens),length(txs))
    ftxs = zeros(Float64,length(intstrens),length(txs))
    for (idx2,tx) in enumerate(txs)
        for (idx,intstren) in enumerate(intstrens)
            pdict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("tx",tx),("ty",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",true),("if_save_data",true)])
            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)

            if if_found
                d,m = read_data(filepath; output_level=0)
                if haskey(m,"fuir_12") && haskey(m,"ftx_12")
                    fuirs[idx,idx2] = m["fuir_12"]
                    ftxs[idx,idx2] = m["ftx_12"]
                else
                    states = d["state"]
                    nrgs = d["nrg"]

                    fuir_2,ogval_uir2 = uir_adiabatic_condition(states[2],states[3],nrgs[3] - nrgs[2],lattice_params,hamilt_params,nothing)
                    fuirs[idx,idx2] = fuir_2

                    ftx_2,ogval_tx2 = tx_adiabatic_condition(states[2],states[3],nrgs[3] - nrgs[2],lattice_params,hamilt_params)
                    ftxs[idx,idx2] = ftx_2
                    
                    datadict = Dict([("fuir_12",fuir_2),("ftx_12",ftx_2)])
                    modify_data(datadict,filepath,"metadata"; output_level=0)
                end
                
            else
                fuir_2,ogval_uir2 = uir_adiabatic_condition(states[2],states[3],nrgs[3] - nrgs[2],lattice_params,hamilt_params,nothing)
                fuirs[idx,idx2] = fuir_2

                ftx_2,ogval_tx2 = tx_adiabatic_condition(states[2],states[3],nrgs[3] - nrgs[2],lattice_params,hamilt_params)
                ftxs[idx,idx2] = ftx_2
                
                datadict = Dict([("fuir_12",fuir_2),("ftx_12",ftx_2)])
                modify_data(datadict,filepath,"metadata"; output_level=0)
            end

            #pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",tx)])
            #all_files = find_data_file(pdict,"ed",dataloc; file_type="jld2")
            #display(all_files)
            #f = all_files[1]
            #d,m = read_data(joinpath(dataloc,f); output_level=0)
            #states = d["state"]
            #nrgs = d["nrg"]
            #lattice_params,hamilt_params = make_latticehamilt_params_from_metadata(m)

            


            #plot_spectrum(txs,nrgs,idx2,pdict["nev"],"Physical Hopping, tx",true; plot_title=" $(lx)x$(ly) N=$n ULR=$intstren")

            
        end
    end#
    #plot_adiabatic_condition(txs,intstrens,fuirs,ftxs; plot_title="$(lx)x$(ly) N=$(n)")
    #=rez = sqrt.(ftxs .^2 .+ fuirs .^ 2)
    fig = figure()
    imshow(log.(rez), origin="lower")
    colorbar()
    ylabel("Interaction Strength")
    xlabel("Hopping Anisotropy")

    fig = figure()
    imshow(log.(fuirs), origin="lower")
    title("Adiabatic Condition ULR")
    colorbar()
    ylabel("Interaction Strength")
    xlabel("Hopping Anisotropy")

    fig = figure()
    imshow(log.(ftxs), origin="lower")
    title("Adiabatic Condition Physical Hopping")
    colorbar()
    ylabel("Interaction Strength")
    xlabel("Hopping Anisotropy")=#
end=#

#= look at ranging Lx and fixed Ly
if false
    ly = 4
    n = 3
    lx = 2*n

    intstrens = range(0.0,10.0,length=21)
    laughlin = 0.0
    for (idx,intstren) in enumerate(intstrens)
        pdict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("if_reading",false),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
        states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)

        fourpt_diag = four_point_diag(states[1],lattice_params; if_plot=false)
        if intstren == 0.0
            global laughlin = maximum(fourpt_diag)
        end
        scatter(intstren,maximum(fourpt_diag) / laughlin,c="b")
        xlabel("Interaction Strength")
        ylabel("Max 4pt Diagonal")
    end

    lx = 8
    n = 4

    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("hopping_anisotropy",1.0)])
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    all_diags = []
    all_intstrens = []
    for f in all_files
        d,m = read_data(dataloc * "/" * f; output_level=0)
        lattice_params = get_lattice_params_from_metadata(m)

        if haskey(m,"fourpt_momentum_diag")
            fourpt_diag = m["fourpt_momentum_diag"]
        elseif haskey(m,"fourpt_momentum")
            fourpt_diag = diag(m["fourpt_momentum"])
        else
            fourpt_diag = four_point_diag(d["state"][1],lattice_params; if_plot=false)
        end

        #scatter(m["U"][end],maximum(fourpt_diag) / laughlin,c="r")
        append!(all_diags,[maximum(fourpt_diag)])
        append!(all_intstrens,[m["U"][end]])
    end

    laughlin_8 = all_diags[findfirst(x -> all_intstrens[x] == 0.0,1:length(all_intstrens))]
    scatter(all_intstrens,all_diags ./ laughlin_8,c="r",label="N=4")
    
end=#

#= plot energy gap ranging intstren and tx
if false
    lx,ly,n = 8,4,4
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge")
    intstrens = [100.0,200.0,300.0,400.0,500.0]
    anises = [1e-4,1e-3,1e-2,1e-1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    #xs_plot = zeros(Float64,length(intstrens)*length(anises))
    #ys_plot = zeros(Float64,length(intstrens)*length(anises))
    #gaps = zeros(Float64,length(intstrens)*length(anises))
    gaps = zeros(Float64,length(intstrens),length(anises))
    for (idx,intstren) in enumerate(intstrens)
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("interaction_strength",intstren)])
        all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
        for f in all_files
            d,m = read_data(joinpath(dataloc,f); output_level=0)

            anis_index = findfirst(x -> anises[x] == m["tx"], 1:length(anises))

            #println("Stren index is $(idx), Anis index is $(anis_index), and linear index is $(length(intstrens)*(idx-1) + anis_index)")

            #xs_plot[length(anises)*(idx-1) + anis_index] = m["tx"]
            #ys_plot[length(anises)*(idx-1) + anis_index] = intstren
            #gaps[length(anises)*(idx-1) + anis_index] = d["nrg"][3] - d["nrg"][1]

            gaps[idx,anis_index] = d["nrg"][3] - d["nrg"][1]
        end
    end

    #bin_count = 100
    #data_dict = bin_values(gaps,bin_count)
    #bv = [data_dict[val] for val in gaps] .* maximum(gaps) / bin_count
    #scatter(xs_plot, ys_plot, c=bv, cmap="viridis")

    pcolormesh(anises,intstrens,gaps; shading="auto",vmin=0.0)
    
    colorbar()
    xlabel("Physical Hopping, "*L"t_x")
    ylabel("IR Interaction Strength, "*L"U_{ir}")
    title("Energy Gap for $(lx)x$(ly) N=$(n)")
end=#

#= energy scaling fixed Ly
if false
    ly = 4
    intstren = 300.0
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/pinned-scaling")
    splittings = []
    gaps = []
    lxs = []
    for n in [3,4,5]
        lx = 2*n
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
        all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
        display(all_files)
        d,m = read_data(dataloc * "/" * all_files[1]; output_level=0)
        append!(lxs,[lx])
        append!(splittings,[d["nrg"][2] - d["nrg"][1]])
        append!(gaps,[d["nrg"][3] - d["nrg"][1]])
    end

    fig = figure()
    scatter(lxs,splittings,c="b")
    xlabel("Lx")
    ylabel("E1 - E0")
    title("Energy Splitting for Ly=$ly, ULR=$intstren")
    yscale("log")

    fig = figure()
    scatter(lxs,gaps,c="r")
    xlabel("Lx")
    ylabel("E2 - E0")
    title("Energy Gap for Ly=$ly, ULR=$intstren")
    ylim([-0.05,1.1*maximum(gaps)])
end=#


#= find if xi_crit depends on system size
if false
    intstren = 2.0

    cols = ["b","g","r","c"]
    for n in [3,4,5]
        lx = n
        ly = 2n
        dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/ulr-length")
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
        all_files = readdir(dataloc)
        filter!(x -> occursin("Lx-$lx",x),all_files)
        display(all_files)
        fig = figure()
        for (idx,f) in enumerate(all_files)
            d,m = read_data(joinpath(dataloc,f); output_level=0)
            xi = m["corr_length"]

            for i in 1:4
                scatter(xi,d["nrg"][i] - d["nrg"][1],c=cols[i])
            end
            
        end
        xlabel("Interaction Length")
        ylabel("E - E0")
        title("Spectrum Lx=$lx, Ly=$ly, N=$n")
    end
end=#

#= redo gamma/omega calcs for all files
if false
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    all_files = filter(x -> occursin("ed",x),readdir(dataloc))
    for (i,f) in enumerate(all_files)
        println(round(100*i/length(all_files),digits=3),"% done")
        bf = jldopen(dataloc * "/" * f,"r")
        if haskey(bf,"all_data")
            close(bf)
            d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        else
            close(bf)
            println("Error with file: $f")
            continue
        end

        if haskey(m,"omega")
            if !haskey(m,"redone_hatsugais")
                println("Working on file: $f")
                ref_multis::Vector{Vector{ComplexF64}} = [zeros(ComplexF64,2) for i in 1:4]
                ref_multis[1:2] = read_data_jld2(dataloc * "/" * m["rm1_name"]; output_level=0)[1]["state"][1:2]
                ref_multis[3:4] = read_data_jld2(dataloc * "/" * m["rm2_name"]; output_level=0)[1]["state"][1:2]
                gamma1,gamma2,omega = get_hatsugaifull(d["state"][1],d["state"][2],ref_multis; if_save=true,filepath=dataloc * "/" * f,ref_multis_filenames=[m["rm1_name"],m["rm2_name"]])
                modify_data_jld2(Dict([("redone_hatsugais",true)]),dataloc * "/" * f,"metadata"; output_level=0)
            end
        end
    end
end

# look at 6x5 n=3 at ULR = 0.0 and 100.0 and look at spectral flow of ground states
if false

    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    nev = 10
    cols = ["b","g","r","m","c"]
    if nev > length(cols)
        cols = repeat(cols,ceil(Int,nev/length(cols)))
    end

    for intstren in [0.0,100.0]
        fig = figure()
        title("Spectral Flow for ULR=$intstren at Theta_y=0.1")
        params_dict = Dict([("Lx",6),("Ly",5),("N",3),("interaction_strength",intstren),("twist_angle2",0.1)])
        all_files = find_data_file(params_dict,"ed",dataloc; output_level=0)
        for f in all_files
            d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
            for i in 1:3#length(d["nrg"])
                scatter(m["twist_angle"][1],d["nrg"][i],c=cols[i])
            end
        end
        xlabel("Theta_x / 2pi")
        ylabel("Energy")
    end
end

# see what sizes are available for rho1D = 1.0 from ED
if false
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    nev = 2
    lx,n = 6,3
    amplification_factor = 10.0

    cols = ["b","g","r","m","c"]
    if nev > length(cols)
        cols = repeat(cols,ceil(Int,nev/length(cols)))
    end

    intstrens = range(0.0,2.0,length=11)
    intstren = 0.0
    tws = range(0.0,1.0,length=21)
    #for (idx,intstren) in enumerate(intstrens)
        #=where_fqh = 1
        if intstren > 1.5
            where_fqh += 2
        end
        if intstren > 2
            where_fqh += 2
        end=#
    #these_nrgs = zeros(Float64,nev,length(tws))
    for ly in [3,5,6]
        fig = figure()
    for (idx,tw1) in enumerate(tws)
    #for (idx2,tw2) in enumerate(tws)
        
        tw2 = 0.0
        #tw1 = 0.0
            params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("nev",nev),("if_save_data",false),("if_find_data",false),("if_periodic_x",true),("if_periodic_y",true),("tw1",tw1),("tw2",tw2)])
            states,nrgs,rhos,filepath,if_found,latpara,hamiltpara = run_normal_ed(params_dict; output_level=1)
            #get_occupancy(states[1],latpara; plot_title="ULR=$intstren tw2=$tw2 E=$(round(nrgs[1],digits=3))")
            #get_occupancy(states[2],latpara; plot_title="ULR=$intstren tw2=$tw2 E=$(round(nrgs[2],digits=3))")
            #=for i in 1:length(nrgs)
                scatter3D(tw1,tw2,nrgs[i] - nrgs[1],c=cols[i])
            end
            xlabel(L"\theta_x / 2 \pi")
            ylabel(L"\theta_y / 2 \pi")
            title("Flux Direction = $(params_dict["flux_direction"])")=#
            #plot_spectrum(intstrens,nrgs,idx,params_dict["nev"],"Interaction Strength",true; plot_title=" $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"])")
            #plot_spectrum(tws,nrgs,idx2,params_dict["nev"],L"\theta_y / 2 \pi",false; plot_title=" $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"]) ULR=$intstren")
            plot_spectrum(tws,nrgs,idx,params_dict["nev"],L"\theta_x / 2 \pi",true; plot_title=" $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"]) ULR=$intstren")
    end

    
    fig = figure()
        for (idx2,tw2) in enumerate(tws)
            
            tw1 = 0.5
            #tw1 = 0.0
                params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("interaction_strength",intstren),("nev",nev),("if_save_data",false),("if_find_data",false),("if_periodic_x",true),("if_periodic_y",true),("tw1",tw1),("tw2",tw2)])
                states,nrgs,rhos,filepath,if_found,latpara,hamiltpara = run_normal_ed(params_dict; output_level=1)
                #get_occupancy(states[1],latpara; plot_title="ULR=$intstren tw2=$tw2 E=$(round(nrgs[1],digits=3))")
                #get_occupancy(states[2],latpara; plot_title="ULR=$intstren tw2=$tw2 E=$(round(nrgs[2],digits=3))")
                #=for i in 1:length(nrgs)
                    scatter3D(tw1,tw2,nrgs[i] - nrgs[1],c=cols[i])
                end
                xlabel(L"\theta_x / 2 \pi")
                ylabel(L"\theta_y / 2 \pi")
                title("Flux Direction = $(params_dict["flux_direction"])")=#
                #plot_spectrum(intstrens,nrgs,idx,params_dict["nev"],"Interaction Strength",true; plot_title=" $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"])")
                plot_spectrum(tws,nrgs,idx2,params_dict["nev"],L"\theta_y / 2 \pi",true; plot_title=" $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"]) ULR=$intstren")
        end
    end
end

# 3d plot of twist angles for 6x5 n=3
if false
    nev = 2
    cols = ["b","g","r","m","c"]
    if nev > length(cols)
        cols = repeat(cols,ceil(Int,nev/length(cols)))
    end
    tws = range(0.0,1.0,length=10)
    intstren = 0.0
    for (idx,tw1) in enumerate(tws)
        for (idx2,tw2) in enumerate(tws)
            params_dict = Dict([("Lx",6),("Ly",5),("N",3),("interaction_strength",intstren),("nev",nev),("if_save_data",false),("if_find_data",false),("if_periodic_x",true),("if_periodic_y",true),("tw1",tw1),("tw2",tw2)])
            states,nrgs,rhos,filepath,if_found,latpara,hamiltpara = run_normal_ed(params_dict; output_level=1)
            for i in 1:length(nrgs)
                scatter3D(tw1,tw2,nrgs[i],c=cols[i])
            end
        end
    end
end

# 4x8 n=4 look at twisting
if false
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    nev = 2
    intstren = 0.0
    #tw2 = 0.1
    params_dict = Dict([("Lx",4),("Ly",8),("N",4),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    all_gs1 = []
    all_gs2 = []
    tw1s = [0.0]
    tw2s = [0.0]
    nrgs = Dict([("1",[0.0]),("2",[0.0]),("3",[0.0])])
    omegas = [0.0*im]
    for f in all_files
        filename_params = get_params_dict_from_filename(f)
        if haskey(filename_params,"twist_angle2")
            #if filename_params["twist_angle2"] != tw2
            #    continue
            #end
            if filename_params["twist_angle2"] == 0.67 || filename_params["twist_angle2"] == 0.33
                continue
            end

            if !isinteger(filename_params["twist_angle2"]*10)
                continue
            end

            if !isinteger(filename_params["twist_angle1"]*10)
                continue
            end

            if filename_params["twist_angle1"] > 1.0
                continue
            end
        else
            continue
        end
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        println(f)
        #append!(all_gs1,[d["state"][1]])
        #append!(all_gs2,[d["state"][2]])
        append!(tw1s,[m["twist_angle"][1]])
        append!(tw2s,[m["twist_angle"][2]])
        for i in 1:nev
            append!(nrgs[string(i)],d["nrg"][i])
        end
        #append!(omegas,[m["omega"]])
    end

    #matrix_omegas = reshape(omegas,11,11)
    #plot_omega(tw2s[1:11],tw2s[1:11],matrix_omegas)

    #scatter(tw1s,nrgs["1"],label="1")
    #scatter(tw1s,nrgs["2"] .- nrgs["1"],label="2")
    #legend()


    #=overlap_11 = abs2.([dot(all_gs1[1],all_gs1[i]) for i in 1:length(all_gs1)])
    overlap_12 = abs2.([dot(all_gs1[1],all_gs2[i]) for i in 1:length(all_gs2)])
    overlap_22 = abs2.([dot(all_gs2[1],all_gs2[i]) for i in 1:length(all_gs2)])
    overlap_21 = abs2.([dot(all_gs2[1],all_gs1[i]) for i in 1:length(all_gs1)])

    fig = figure()
    plot(tws[1:end],overlap_11,"-p",label="11")
    plot(tws[1:end],overlap_22,"-p",label="22")
    legend()

    fig = figure()
    plot(tws[1:end],overlap_12,"-p",label="12")
    plot(tws[1:end],overlap_21,"-p",label="21")
    legend()=#
end

# compare twisting for ED/1DeffMPS 4x4 n=2
if false-
    nev = 3
    #intstrens = range(0.0,4.0,length=21)
    tws = range(-0.1,0.1,length=21)
    for (idx,tw1) in enumerate(tws)
        intstren = 100000.0
    #prev_next_fqh = [[],[]]
    #for (idx,intstren) in enumerate(intstrens)
        params_dict = Dict([("Lx",3),("Ly",7),("N",3),("nev",nev),("tw2",tw1),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("if_find_data",false),("if_save_data",false)])
        states,nrgs,rhos,filepath,if_found,latpara,hamiltpara = run_normal_ed(params_dict; output_level=1)
        plot_spectrum(tws,nrgs,idx,params_dict["nev"],L"\theta_y / 2 \pi",false; plot_title=" $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"]) ULR=$(params_dict["interaction_strength"])")
        #plot_spectrum(intstrens,nrgs,idx,params_dict["nev"],"Interaction Strength",true; plot_title=" $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"]) Theta_y=$(params_dict["tw2"])")
        #=if idx == 1
            prev_next_fqh[2] = states[2]
            append!(prev_next_fqh[1],[2.0])
        else
            overlaps_vals = [abs2(transpose(conj.(prev_next_fqh[2])) * states[i]) for i in 1:length(states)]
            append!(prev_next_fqh[1],[float(findfirst(x -> overlaps_vals[x] == maximum(overlaps_vals),1:length(overlaps_vals)))])
            prev_next_fqh[2] = states[Int(prev_next_fqh[1][end])]
        end
        scatter(intstren,nrgs[Int(prev_next_fqh[1][end])] - nrgs[1],c="k")=#
    end


end=#

#= sanity check that twistings are making different hams
if false
    lx,ly,n = 8,4,4
    intstren = 0.0
    
    tw11,tw21 = 0.2,0.3
    pdict1 = Dict([("Lx",lx),("Ly",ly),("N",n),("nev",nev),("filling",0.5),("tw2",tw21),("tw1",tw11),("lr","all"),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("if_find_data",false),("if_save_data",false)])
    lats1,hams1,runs1 = get_normal_model_params_ed(pdict1)
    lats1["full_basis"] = n_particle_basis(lats1; output_level=0,dataloc=runs1.basis_dataloc)
    fullham1 = buildHam(lats1,hams1; output_level=1)

    tw12,tw22 = 0.3,0.2
    pdict2 = Dict([("Lx",lx),("Ly",ly),("N",n),("nev",nev),("filling",0.5),("tw2",tw22),("tw1",tw12),("lr","all"),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("if_find_data",false),("if_save_data",false)])
    lats2,hams2,runs2 = get_normal_model_params_ed(pdict2)
    lats2["full_basis"] = n_particle_basis(lats2; output_level=0,dataloc=runs2.basis_dataloc)
    fullham2 = buildHam(lats2,hams2; output_level=1)
end=#

#= hatsugai data collection from including ed data
if false
    lx,ly,n = 6,5,3
    nev = 10
    #tws = range(-0.5,0.5,length=11)
    #tws2 = range(-0.75,-0.25,length=11)
    tws = range(0.0,1.0,length=11)
    tws2 = tws
    intstren = 100.0

    lambda1s::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
    lambda2s::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
    omegas::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
    ref_multis,rm1,rm2 = get_reference_multiplets(lx,ly,n; interaction_strength=intstren,if_save=false)
    for (idx,tw1) in enumerate(tws)
        for (idx2,tw2) in enumerate(tws2)
            params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("nev",nev),("filling",0.5),("tw2",tw2),("tw1",tw1),("lr","all"),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("if_find_data",false),("if_save_data",false)])
            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)

            lambda1s[idx,idx2],lambda2s[idx,idx2],omegas[idx,idx2] = get_hatsugaifull(states[1],states[2],ref_multis; if_save=false,filepath=filepath,ref_multis_filenames=[rm1,rm2])
        end
    end#
    #xlabel(L"\theta_x / 2 \pi")
    #ylabel(L"\theta_y / 2 \pi")
    #zlabel("Energy")
    #title("Energy Spectrum $(lx)x$(ly) N=$(n) ULR=$intstren")

    plot_omega(tws,tws2,omegas; plot_title=" $(lx)x$(ly) N=$(n) ULR=$intstren")
    plot_gamma(tws,tws2,lambda1s,1; plot_title=" $(lx)x$(ly) N=$(n) ULR=$intstren")
    plot_gamma(tws,tws2,lambda2s,2; plot_title=" $(lx)x$(ly) N=$(n) ULR=$intstren")

end=#

#=function get_closing_strength(intstrens,flatnesses::Vector{Float64})
    if length(intstrens) < 4
        which_y = 2
    else
        which_y = findfirst(i -> (flatnesses[i+1] - flatnesses[i]) / (flatnesses[i] - flatnesses[i-1]) < 1e-2,2:length(flatnesses)-1)
    end
    #scatter(intstrens[which_y],flatnesses[which_y],c="r")
    return (1.0 - flatnesses[findfirst(x -> x == 0.0,intstrens)]) / ((flatnesses[which_y] - flatnesses[findfirst(x -> x == 0.0,intstrens)]) / intstrens[which_y]),which_y
end

# finite size scaling of flatness
if false
    # ones we have (3,8,3),(4,8,4)
    # ED ones we want (5,8,5)
    # TTN ones we want (8,8,8)
    xs,ys = twist_flatness_ed(4,8,4)
    closing_stren,which_y = get_closing_strength(xs,ys)
    scatter(closing_stren,1.0,c="r",label="Guess")
end

# hatsugai data collection from existing ed data
if false
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    lx,ly,n = 8,3,3
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0)

    for (idx,f) in enumerate(all_files)
        filename_dict = get_params_dict_from_filename(f)
        if haskey(filename_dict,"twist_angle1") && (filename_dict["twist_angle1"] == 0.33 || filename_dict["twist_angle2"] == 0.67 || filename_dict["twist_angle1"] == 0.67 || filename_dict["twist_angle2"] == 0.33)
            continue
        end
        println(round(100*idx/length(all_files),digits=3),"% done")
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        if !haskey(m,"omega")
            println("Does not have Hatsugai data")
            ref_multis,rm1,rm2 = get_reference_multiplets(lx,ly,n; interaction_strength=filename_dict["interaction_strength"],if_make_new=true)
            lambda1,lambda2,omega = get_hatsugaifull(d["state"][1],d["state"][2],ref_multis; if_save=true,filepath=dataloc * "/" * f,ref_multis_filenames=[rm1,rm2])
        end
    end
end

# phase diagram using flatness ULR vs rho1D
if false
    # ones we want (8,5,5)
    configs = [(3,8,3),(8,3,3),(4,6,3),(8,4,4),(8,5,5)]
    #rez = make_phasediag_ulrrho1d_flatness(configs; max_intstren=500.0, if_plot=false)
    plot_phasediag_ulrrho1d_flatness(rez...)
end

# save denscorrs for all configs at all interaction strengths
if false
    configs = [(8,3,3),(4,6,3),(3,8,3),(8,4,4)]
    for (lx,ly,n) in configs
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
        dataloc = get_folder_location("cluster-data/exact-diag/torus")
        all_files = find_data_file(pdict,"ed",dataloc; output_level=0)
        
        filter!(x -> !occursin("twist_angle1",x),all_files)
        filter!(x -> !occursin("mk",x),all_files)

        intstrens = Float64[]
        results = Float64[]
        for f in all_files

            display(f)

            filepath = dataloc * "/" * f
            d,m = read_data_jld2(filepath; output_level=0)

            if !haskey(m,"dens_corr_mat")
                make_density_correlations(d["state"][1],get_lattice_params_from_metadata(m); if_save=true,filepath=filepath)
            end
            println("Finished $(lx)x$(ly) N=$(n) U=$(m["U"][end])")

        end
    end

end

# getting density-density correlation matrix at given angle for all phase diag configurations
function datacollection_dd(mom_angle::Float64,configs::Vector{Tuple{Int64,Int64,Int64}}=[(8,3,3),(4,6,3),(3,8,3),(8,4,4),(8,5,5)]; kwargs...)
    hanis::Float64 = get(kwargs,:hopping_anisotropy,1.0)
    if_redo::Bool = get(kwargs,:if_redo,false)
    if_plot::Bool = get(kwargs,:if_plot,false)

    for (lx,ly,n) in configs
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",hanis)])
        dataloc = get_folder_location("cluster-data/exact-diag/torus")
        all_files = find_data_file(pdict,"ed",dataloc; output_level=0)
        
        filter!(x -> !occursin("twist_angle1",x),all_files)
        filter!(x -> !occursin("mk",x),all_files)

        intstrens = Float64[]
        results = Float64[]
        for f in all_files

            filepath = dataloc * "/" * f
            d,m = read_data_jld2(filepath; output_level=0)

            push!(intstrens,m["U"][end])

            println("Working on $(lx)x$(ly) n=$n at Interaction Strength: $(m["U"][end])")

            latparas = get_lattice_params_from_metadata(m)

            if !haskey(m,"dens_corr_mat")
                denscorrs = make_density_correlations(d["state"][1],latparas; if_save=true,filepath=filepath)
            else
                denscorrs = m["dens_corr_mat"]
            end

        end
    end
    
end

# phase diag with ft-dd ratio
if false
    configs = [(8,4,4)]#[(8,3,3),(4,6,3),(3,8,3)]#,(8,4,4)]#,(8,5,5)]

    ulrs::Vector{Float64} = Float64[]
    ftdds::Vector{Float64} = Float64[]
    oneDrhos::Vector{Float64} = Float64[]

    #=for (idx,config) in enumerate(configs)
        lx,ly,n = config
        append!(flatnesses,[twist_flatness_1deff(lx,ly,n; if_plot_spectrum=false,plot_title=" $(lx)x$ly N=$n")])
        append!(ulrs,[400.0])
        append!(oneDrhos,[n/lx])
    end=#

    for (lx,ly,n) in configs
        local_strens,local_ftdd = get_ftdd_ratio(lx,ly,n)
        append!(ulrs,local_strens)
        append!(ftdds,local_ftdd ./ (lx*ly))#./ minimum(local_ftdd))
        append!(oneDrhos,ones(Float64,length(local_strens)) .* (n / lx))
    end

    bin_count = 100
    data_dict = bin_values(ftdds,bin_count)
    bv = [data_dict[val] for val in ftdds]
    min_nrgs2, max_nrgs2 = minimum(ftdds),maximum(ftdds)
    normalized_bv = [(val - minimum(bv)) / (maximum(bv) - minimum(bv)) * (max_nrgs2 - min_nrgs2) + min_nrgs2 for val in bv]

    fig = figure()
    scatter(oneDrhos, ulrs, c=normalized_bv, cmap="plasma")
    colorbar()
    xlabel(L"\rho_{1D}")
    ylabel("ULR")
    title("FT-DD Ratio Phase Diagram")
    #yscale("log")

end

# generic size playing with FT-DD rectangle/square
if false
    ks = collect(range(0*pi/2,2*pi/1,length=100))
    intstren = 0.0
    lx,ly,n = 8,4,4
    allowed_momenta_physical = filter(x-> x<2*pi+0.01,[2*pi*i / (lx/1) for i in 0:20])
    allowed_momenta_synthetic = filter(x-> x<2*pi+0.01,[2*pi*i / (ly/1) for i in 0:20])
    params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren)])
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    all_files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    filter!(x -> !occursin("twist_angle1",x),all_files)
    filter!(x -> !occursin("mk",x),all_files)

    if length(all_files) == 0
        println("No data found")
        params_dict["if_save_data"] = false
        params_dict["if_find_data"] = false
        states,nrgs,rhos,filepath,if_found,latparas,hamilt_params = run_normal_ed(params_dict; output_level=1)
        denscorrs = make_density_correlations(states[1],latparas)
    else
        display(all_files)

        f = all_files[1]
        filepath = dataloc * "/" * f
        d,m = read_data_jld2(filepath; output_level=0)

        latparas = get_lattice_params_from_metadata(m)
        denscorrs = m["dens_corr_mat"]
    end

    ftdd_results = zeros(Float64,length(ks),length(ks))
    for (idx,kx) in enumerate(ks)
        for (idx2,ky) in enumerate(ks)
            ftdd = ft_densitydensity_correlation([kx,ky],nothing,latparas; denscorrs=denscorrs,if_save=false,filepath=filepath)
            ftdd_results[idx,idx2] = abs(ftdd)
        end
    end

    norm_factor = integrate_2d_matrix(ftdd_results)
    ftdd_results = ftdd_results ./ norm_factor

    #=if lx == 4
        append!(plotting_ys,[[ftdd_results[1,:]]])
    elseif lx == 8
        append!(plotting_ys,[[ftdd_results[:,1]]])
    end=#

    #=allow_physmom_index = zeros(Int,length(allowed_momenta_physical))
    for (idx,ak) in enumerate(allowed_momenta_physical)
        allow_physmom_index[idx] = findfirst(x -> isapprox(x,ak;atol=1e-2),ks)
    end
    allow_synmom_index = zeros(Int,length(allowed_momenta_synthetic))
    for (idx,ak) in enumerate(allowed_momenta_synthetic)
        allow_synmom_index[idx] = findfirst(x -> isapprox(x,ak;atol=1e-2),ks)
    end=#

    fig = figure()
    imshow(ftdd_results,extent=(minimum(ks),maximum(ks),minimum(ks),maximum(ks)),origin="lower")
    xlabel(L"k_x")
    ylabel(L"k_y")
    colorbar()
    title("FT-DD for $(lx)x$(ly) N=$(n) ULR=$intstren")

    fig = figure()
    these_ys = ftdd_results[1,:]
    plot(ks ./ pi,these_ys)
    for k in allowed_momenta_synthetic
        plot([k / pi,k/pi],[minimum(these_ys),maximum(these_ys)],c="r")
    end
    xlabel(L"k_x")
    title("FT-DD at ky=0 for $(lx)x$(ly) N=$(n) ULR=$intstren")

    fig = figure()
    these_ys = ftdd_results[:,Int(round(length(ks)/2,digits=0))]
    plot(ks ./ pi,these_ys)
    for k in allowed_momenta_physical
        plot([k / pi,k/pi],[minimum(these_ys),maximum(these_ys)],c="r")
    end
    xlabel(L"k_y")
    title("FT-DD at kx=pi for $(lx)x$(ly) N=$(n) ULR=$intstren")

    get_distDDcorrs(denscorrs,"both"; plot_title=" $(lx)x$(ly) N=$(n) ULR=$intstren")

    #get_occupancy(d["state"][1],latparas; plot_title=" $(lx)x$(ly) N=$(n) ULR=$intstren")
    #end

    #=fig, ax = subplots()

    x1 = collect(ks)
    x2 = collect(ks)
    y1 = plotting_ys[1][1]
    y2 = plotting_ys[2][1]

    # Plot the first line
    ax.plot(x1, y1, c="b", label=L"\rho_{1D}"*"=1.0")
    ax.set_xlabel(L"k_x")
    ax.set_ylabel("FT-DD at ULR=400.0")
    ax.tick_params(axis="x", colors="b")

    # Create a twin x-axis
    ax2 = ax.twiny()

    # Plot the second line on the twin x-axis
    ax2.plot(x2, y2, color="r", label=L"\rho_{1D}"*"=0.5")
    ax2.set_xlabel(L"k_y")
    ax2.tick_params(axis="x", colors="r")

    # Add a legend
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(loc="upper left")
    ax2.legend(loc="upper right")

    # Adjust the layout
    tight_layout()=#
end

# phase transition by fourier transform of density profile at k=(pi,0) for given configs
if false
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

        #=if intstrens[end] == 0.0 || intstrens[end] == 1000.0 || (lx == 10 && isapprox(200.0,intstrens[end],atol=10.0))
            if_plot = true
        else
            if_plot = false
        end=#

        #=if if_plot
            if_plot = false
        else
            if_plot = true
        end=#

        latparas = get_lattice_params_from_metadata(m)
        occs = get_occupancy(d["state"][1],latparas; if_plot=if_plot,plot_title="$(lx)x$(ly) n=$n ULR=$(m["U"][end])")
        ftd = ft_density([pi,0],occs)
        append!(cdw_sfs,[abs(ftd)])

    end

    #fig = figure()
    scatter(intstrens,cdw_sfs,label="$(lx)x$(ly)")
    xlabel("Interaction Strength")
    ylabel("CDW SF at k=(pi,0)")
    title("CDW SF at "*L"\rho_{1D}=1/2")
    xscale("log")
    legend()
    end
end

# look at FT-D for ks range for other rho1D
if false
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    ks = range(0,2*pi,length=50)
    lx,ly,n = 10,5,5
    if_pinning = true
    #=pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0)

    filter!(x -> !occursin("twist_angle1",x),all_files)
    filter!(x -> !occursin("mk",x),all_files)
    filter!(x -> occursin("if_pinning",x),all_files)

    filter!(x -> occursin("interaction_strength-1000.0-",x),all_files)
    f = all_files[1]

        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)=#
        latparas = get_lattice_params_from_metadata(m)
        occs = get_occupancy(d["state"][1],latparas; if_plot=true,plot_title="$(lx)x$(ly) n=$n ULR=$(m["U"][end])")

        ftds = zeros(Float64,length(ks),length(ks))
        for (idx,kx) in enumerate(ks)
            for (idx2,ky) in enumerate(ks)
                ftd = ft_density([kx,ky],occs)
                ftds[idx,idx2] = abs(ftd)
            end
        end
        fig = figure()
        imshow(ftds,extent=(minimum(ks),maximum(ks),minimum(ks),maximum(ks)),origin="lower")
        xlabel(L"k_y")
        ylabel(L"k_x")
        colorbar()
        title("FT-Density for $(lx)x$(ly) N=$(n) ULR=$(m["U"][end])")
end

# flatness finite size scaling at Lx by 4 with rho_1D = 1/2
if false
    configs = [(4,4,2),(6,4,3),(8,4,4)]
    flatnesses = Dict([("4",[]),("6",[]),("8",[])])
    for (lx,ly,n) in configs

        if lx == 8
            tws = range(0.0,1.0,length=11)
        else
            tws = range(0.0,1.0,length=21)
        end

        for intstren in [0.0]
            #all_nrgs = Dict([("1",zeros(Float64,length(tws),length(tws))),("2",zeros(Float64,length(tws),length(tws))),("3",zeros(Float64,length(tws),length(tws)))])
            #fts = zeros(Float64,length(tws),length(tws))
            site_occs = Dict([("11",zeros(Float64,length(tws),length(tws))),("22",zeros(Float64,length(tws),length(tws))),("12",zeros(Float64,length(tws),length(tws))),("21",zeros(Float64,length(tws),length(tws)))])
            for (idx,tw1) in enumerate(tws)
                for (idx2,tw2) in enumerate(tws)
                    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("tw1",tw1),("tw2",tw2),("hopping_anisotropy",1.0),("if_pinning",true),("interaction_strength",intstren),("filling",0.5),("lr","all"),("if_find_data",true),("if_save_data",false),("nev",10)])
                    states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)
                    occs = get_occupancy(states[1],lattice_params; if_plot=false,plot_title="ULR=$(intstren) tw1=$tw1 tw2=$tw2")
                    site_occs["11"][idx,idx2] = occs[1,1]
                    site_occs["22"][idx,idx2] = occs[2,2]
                    site_occs["12"][idx,idx2] = occs[1,2]
                    site_occs["21"][idx,idx2] = occs[2,1]
                    #fts[idx,idx2] = real(ft_density([pi,0],occs))
                    #all_nrgs["1"][idx,idx2] = nrgs[1]
                    #all_nrgs["2"][idx,idx2] = nrgs[2]
                    #all_nrgs["3"][idx,idx2] = nrgs[3]
                end
            end

            for (k,v) in site_occs
                fig = figure()
                imshow(v,extent=(minimum(tws),maximum(tws),minimum(tws),maximum(tws)),origin="lower")
                xlabel(L"\theta_x / 2 \pi")
                ylabel(L"\theta_y / 2 \pi")
                colorbar()
                title("Occupancy for $(lx)x$(ly) N=$(n) ULR=$intstren k=$k")
            end
            #=fig = figure()
            imshow(fts,extent=(minimum(tws),maximum(tws),minimum(tws),maximum(tws)),origin="lower")
            xlabel(L"\theta_x / 2 \pi")
            ylabel(L"\theta_y / 2 \pi")
            colorbar()
            title("FT-Density for $(lx)x$(ly) N=$(n) ULR=$intstren")=#

            #flatness = maximum((all_nrgs["2"] .- all_nrgs["1"]) ./ (all_nrgs["3"] .- all_nrgs["1"]))
            #append!(flatnesses[string(lx)],[flatness])
        end
    end

    #=fig = figure()
    for (lx,flatness) in flatnesses
        if lx == 4
            plot([lx],flatness[1],c="b",label="ULR=0.0")
            plot([lx],flatness[2],c="r",label="ULR=1000.0")
        else
            plot([lx],flatness[1],c="b")
            plot([lx],flatness[2],c="r")
        end
    end
    xlabel("Lx")
    ylabel("Flatness")
    title("Flatness Finite Size Scaling")
    legend()=#

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
end

# looking at laughlin finite scaling at ULR=0.0
if false
    configs = [(4,4,2),(6,4,3),(8,4,4),(10,4,5)]
    for config in configs
        lx,ly,n = config
        pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("if_pinning",false),("hopping_anisotropy",1.0),("interaction_strength",0.0),("filling",0.5),("lr","all"),("if_find_data",true),("if_save_data",false),("nev",10)])
        states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)

        get_occupancy(states[1],lattice_params; if_plot=true,plot_title="$(lx)x$(ly) n=$n ULR=0.0")
        get_occupancy(states[2],lattice_params; if_plot=true,plot_title="E1 - E0 = $(round(nrgs[2] - nrgs[1],digits=5))")
    end
end

# do density-density correlations like Alberto
if false
    lx,ly,n = 8,4,4
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0)

    filter!(x -> !occursin("twist_angle1",x),all_files)
    filter!(x -> !occursin("mk",x),all_files)
    filter!(x -> !occursin("pinning",x),all_files)

    filter!(x -> occursin(".0-",x),all_files)

    #= for when has saved non-Alberto version of density-density
    for f in all_files
        d,m = read_data_jld2(f,dataloc; output_level=0)
        if haskey(m,"dens_corr_mat")
            println(f)
            dd = m["dens_corr_mat"]
        end
    end=#

    us = []
    fts = []
    #for f in all_files
    f = all_files[1]
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        #=if !(m["U"][end] in [0.0,1.0,20.0,100.0])
            continue
        end=#
        latparas = get_lattice_params_from_metadata(m)
        #occs = get_occupancy(d["state"][1],latparas; if_plot=false,plot_title="$(lx)x$(ly) n=$n ULR=$(m["U"][end])")
        #=dds = fourpoint_alberto(d["state"][1],latparas; occs=occs, plot_title="$(lx)x$(ly) n=$n ULR=$(m["U"][end])",if_plot=false)
        mdata = Dict([("densitydensity",dds),("occs",occs)])
        modify_data_jld2(mdata,dataloc * "/" * f,"metadata"; output_level=1)
        ft_dds_stripe = abs(ft_densitydensity([pi,0],dds))
        append!(us,m["U"][end])
        append!(fts,ft_dds_stripe)=#
        #=moms = range(0.0,2*pi,length=100)
        ft_dds = zeros(Float64,length(moms),length(moms))
        for (idx,kx) in enumerate(moms)
            for (idx2,ky) in enumerate(moms)
                ft_dds[idx2,idx] = abs(ft_densitydensity([kx,ky],dds))
            end
        end
        #display(ft_dds)
        fig = figure()
        imshow(ft_dds,extent=(minimum(moms),maximum(moms),minimum(moms),maximum(moms)),origin="lower")
        xlabel(L"k_{phys}")
        ylabel(L"k_{synth}")
        colorbar()
        title("FT-Density-Density for $(lx)x$(ly) N=$(n) ULR=$(m["U"][end])")=#

        dds = m["densitydensity"]
        #centersite = [Int64(ceil(lx/2)),Int64(ceil(ly/2))]
        #pairdist = pairdistribution(dds,occs; if_plot=true,plot_title="$(lx)x$(ly) n=$n ULR=$(m["U"][end])")
        #=fig = figure()
        plot(1:lx,pairdist[centersite[2],:],"-p",label=m["U"][end])
        xlabel("Physical Site")

        #plot(1:ly,pairdist[:,centersite[1]],"-p",label=m["U"][end])
        #xlabel("Synthetic Site")

        ylabel("Pairdist")
        title("Pairdist Slice")
        legend()=#

    #end
    #=fig = figure()
    scatter(us,fts)
    xlabel("Interaction Strength")
    ylabel("FT-DD at k=("*L"\pi"*",0)")
    title("FT-DD at k=("*L"\pi"*",0) for $(lx)x$(ly) N=$(n)")=#
end

# checking 10x5 N=5 unpinned data
if false
    lx,ly,n = 10,5,5
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0)

    display(all_files)
    for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        latparas = get_lattice_params_from_metadata(m)
        occs = m["occs"]
        dds = fourpoint_alberto(d["state"][1],latparas; occs=occs, plot_title="$(lx)x$(ly) n=$n ULR=$(m["U"][end])",if_plot=true)
        modify_data_jld2(Dict([("densitydensity",dds)]),filepath,"metadata"; output_level=1)
    end
end

# spatial entanglement bipartition general
if false
    cols = ["b","r","g"]
    for (idx,intstren) in enumerate([0.0])#,1.0,100.0])
    #for (idx,intstren) in enumerate(range(0.0,20.0,length=10))
        if true
            lx,ly,n = 8,4,4
            #intstren = 0.0
            pdict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",2),("if_find_data",true),("if_save_data",false)])
            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=1)
        end

        psi = states[1]
        
        perims = [6,8,10,12,14,16]
        all_subsysts_phys = [[1,2],[1,2,3],[1,2,3,7,8,9],[1,2,3,4,7,8,9,10],[1,2,3,4,5,7,8,9,10,11],[1,2,3,4,5,6,7,8,9,10,11,12]]
        all_subsysts_vert = [[1,7],[1,7,13],[1,2,7,8,13,14],[1,2,3,7,8,9,13,14,15],[1,2,3,4,7,8,9,10,13,14,15,16],[1,2,3,4,5,7,8,9,10,11,13,14,15,16,17]]
        for (i,subsysA) in enumerate(all_subsysts_phys)
            subsysB = find_subsystem_B(subsysA,lx,ly)
            rho_A, unique_A_configs = compute_reduced_density_matrix(psi, lattice_params["full_basis"], subsysA, subsysB)
            ee = entanglement_entropy(rho_A)
            if i == 1
                scatter(perims[i],ee,c="b",label="Physical")
            else
                scatter(perims[i],ee,c="b")
            end
        end
        for (i,subsysA) in enumerate(all_subsysts_vert)
            subsysB = find_subsystem_B(subsysA,lx,ly)
            rho_A, unique_A_configs = compute_reduced_density_matrix(psi, lattice_params["full_basis"], subsysA, subsysB)
            ee = entanglement_entropy(rho_A)
            if i == 1
                scatter(perims[i],ee,c="r",label="Synthetic")
            else
                scatter(perims[i],ee,c="r")
            end
        end
        xlabel("Perimeter")
        ylabel("Entanglement Entropy")
        title("Entanglement Entropy for Oriented Bipartitions 6x3 N=3 ED")
        legend()
        xlim(0,18)
        ylim(-1,3)

        #=all_subsysts = [[1,2,3],[1,7,13],[1,2,7,8]]
        for (i,subsysA) in enumerate(all_subsysts)
            subsysB = find_subsystem_B(subsysA,lx,ly)
            rho_A, unique_A_configs = compute_reduced_density_matrix(psi, lattice_params["full_basis"], subsysA, subsysB)
            ee = entanglement_entropy(rho_A)
            if idx == 1
                scatter(intstren,ee,c=cols[i],label="$subsysA")
            else
                scatter(intstren,ee,c=cols[i])
            end
        end
        xlabel("Interaction Strength")
        ylabel("Entanglement Entropy")
        legend()=#
    end

end

# find bipartition entanglement entropy for rho_1D = 1.0 hopefully zero. It is not
if false
    lx,ly,n = 4,8,4
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0)

    keepints = [0.0,0.52632,2.0,100.0,400.0]
    cols = ["b","r","g"]
    for f in all_files

        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

        if !(m["U"][end] in keepints)
            continue
        end

        latparas = get_lattice_params_from_metadata(m)
        psi = d["state"][1]

        perims = [6,8,12,16]
        all_subsysts = [[1,5],[1,2,5,6],[1,2,5,6,9,10,13,14],[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]]

        if m["U"][end] == 0.0
            col = "b"
        elseif m["U"][end] == 0.52632
            col = "r"
        elseif m["U"][end] == 2.0
            col = "g"
        elseif m["U"][end] == 100.0
            col = "k"
        elseif m["U"][end] == 400.0
            col = "m"
        end


        for (i,subsysA) in enumerate(all_subsysts)
            subsysB = find_subsystem_B(subsysA,lx,ly)
            rho_A, unique_A_configs = compute_reduced_density_matrix(psi, latparas["full_basis"], subsysA, subsysB)
            ee = entanglement_entropy(rho_A)
            if i == 1
                scatter(perims[i],ee,c=col,label="$(m["U"][end])")
            else
                scatter(perims[i],ee,c=col)
            end
        end
    end
    legend()
    xlabel("Perimeter")
    ylabel("Entanglement Entropy")

end

# find TEE of superfluid/charge ordered states 
if false
    lx,ly,n = 8,4,4
    intstren = 0.0

    params_dict = Dict([("output_level",1),("if_check_fluxes",false),("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("lr","all"),("alpha",0.0),("nev",10),("if_find_data",false),("if_save_data",false)])
    states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)
    psi = states[1]

    occs = get_occupancy(states[1],lattice_params; plot_title=" $(lx)x$ly N=$n Superfluid")

    #all_subsysts = [[1,2],[1,2,7,8],[1,2,3,7,8,9],[1,2,3,4,7,8,9,10]]
    #perims = [6,8,10,12]
    #all_subsysts = [[1],[1,2],[1,2,5,6],[1,2,3,5,6,7],[1,2,3,4,5,6,7,8]]
    #perims = [4,6,8,10,12]
    all_subsysts = [[1,2],[1,2,9,10],[1,2,3,9,10,11],[1,2,3,4,9,10,11,12]]
    perims = [6,8,10,12]
    ees = []

    fig = figure()
    for (i,subsysA) in enumerate(all_subsysts)
        subsysB = find_subsystem_B(subsysA,lx,ly)
        rho_A, unique_A_configs = compute_reduced_density_matrix(psi, lattice_params["full_basis"], subsysA, subsysB)
        ee = entanglement_entropy(rho_A)
        append!(ees,[ee])
        scatter(perims[i],ee,c="b")
    end
    
    xlabel("Perimeter")
    ylabel("Entanglement Entropy")
    xlim([0,13])
    ylim([-1,1.1*maximum(ees)])

    linfit = linear_fit(perims,ees)
    plot(range(0,13,length=2),linfit[1] .+ range(0,13,length=2) .* linfit[2],c="r")

    #title("TEE for Charge Ordered State, Anis=$(params_dict["hopping_anisotropy"])")
    title("TEE for Superfluid, "*L"\gamma"*"=$(round(-linfit[1],digits=3))")


end

# FQH optimal control comparison
if false
    lx,ly,n = 6,6,3

    #=alphas = range(0.15,0.35,length=20)
    bds = zeros(Float64,length(alphas))
    for (idx,alph) in enumerate(alphas)
        params_dict = Dict([("output_level",1),("if_check_fluxes",false),("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",false),("if_periodic_y",false),("hopping_anisotropy",1.0),("interaction_strength",0.0),("lr",0),("alpha",alph),("nev",10),("if_find_data",false),("if_save_data",false)])
        states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)
        psi = states[1]

        scatter(alph,nrgs[1] - nrgs[1],c="b")
        scatter(alph,nrgs[2] - nrgs[1],c="g")
        scatter(alph,nrgs[3] - nrgs[1],c="r")

        #=occs = get_occupancy(psi,lattice_params; if_plot=false)
        bulkdens = sum(occs[3:4,3:4]) / 4
        bds[idx] = bulkdens
        scatter(alph,bulkdens,c="b")=#
    end
    xlabel(L"\alpha")
    #ylabel("Bulk Density")
    ylabel("Energy Gap")

    #linfit = linear_fit(alphas,bds)
    #plot(range(0.2,0.25,length=2),linfit[1] .+ range(0.2,0.25,length=2) .* linfit[2],c="r",label="$(round(linfit[2],digits=3))")
    =#

    tws = range(0.0,1.0,length=11)
    for (idx,tw1) in enumerate(tws)
        for (idx2,tw2) in enumerate(tws)
            params_dict = Dict([("output_level",1),("if_check_fluxes",false),("tw1",tw1),("tw2",tw2),("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",false),("if_periodic_y",false),("hopping_anisotropy",1.0),("interaction_strength",0.0),("lr",0),("alpha",0.22),("nev",10),("if_find_data",false),("if_save_data",false)])
            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)
            psi = states[1]

            scatter3D(tw1,tw2,nrgs[2] - nrgs[1],c="b")
            scatter3D(tw1,tw2,nrgs[3] - nrgs[1],c="r")
            xlabel(L"\theta_x")
            ylabel(L"\theta_y")
        end
    end

end=#

#= Felix Palm, looking at interaction length and fourpt momentum flatness as order parameter
if true
    lx,ly,n = 8,4,4
    dataloc = get_folder_location("cluster-data/exact-diag/torus/new-gauge/ulr-length")
    pdict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ed",dataloc; output_level=0,file_type="jld2")
    
    xis = range(0.0,4.0,length=11)
    intstrens = range(0.0,4.0,length=11)[2:end]
    flatnesses = ones(Float64,length(intstrens),length(xis))
    gaps = ones(Float64,length(intstrens),length(xis))
    splittings = ones(Float64,length(intstrens),length(xis))
    for f in all_files

        d,m = read_data(joinpath(dataloc,f); output_level=0)

        xi = m["corr_length"]

        xi_index = findfirst(x -> xis[x] == xi,1:11)
        intstren_index = findfirst(x -> intstrens[x] == m["U"][1],1:11)

        #=if !haskey(m,"fourpt_momentum_1")
            continue
        end

        fourpt1 = m["fourpt_momentum"]
        fourpt2 = m["fourpt_momentum_1"]

        fourpt_mixed = 0.5 * (fourpt1 .+ fourpt2)

        subset_fourpt = vcat([diag(fourpt_mixed,i) for i in 3:lx-3]...)
        flatness = minimum(subset_fourpt) / maximum(subset_fourpt)

        

        flatnesses[intstren_index,xi_index] = flatness=#


        gaps[intstren_index,xi_index] = d["nrg"][3] - d["nrg"][2]
        splittings[intstren_index,xi_index] = d["nrg"][2] - d["nrg"][1]

    end

    #= normalize from laughlin tao-thouless
    for i in 1:length(intstrens)
        flatnesses[i,:] ./= flatnesses[i,1]
    end

    imshow(flatnesses,extent=(minimum(xis),maximum(xis),minimum(intstrens),maximum(intstrens)),origin="lower",aspect="auto",vmin=0.0,vmax=1.0)
    xlabel("Interaction Length")
    ylabel("Interaction Strength")
    title("Normalized Fourpt Flatness for 8x4 N=4")
    colorbar()=#

    fig = figure()
    imshow(gaps,extent=(minimum(xis),maximum(xis),minimum(intstrens),maximum(intstrens)),origin="lower",aspect="auto",vmin=0.0)
    xlabel("Interaction Length")
    ylabel("Interaction Strength")
    title("Energy Gap")
    colorbar()

    fig = figure()
    imshow(log.(splittings),extent=(minimum(xis),maximum(xis),minimum(intstrens),maximum(intstrens)),origin="lower",aspect="auto")
    xlabel("Interaction Length")
    ylabel("Interaction Strength")
    title("Groundstate Splitting")
    colorbar().set_label(L"log_{10} (\Delta_{01})")
end=#

if true
    lx,ly,n = 8,4,4
    intstrens = range(0.0,5.0,length=11)
    for (idx,intstren) in enumerate(intstrens)
        params_dict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",false),("if_periodic_y",false),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])
        states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)

        if idx == 1
            occs = get_occupancy(states[1],lattice_params; if_plot=true,plot_title="$(lx)x$(ly) N=$n ULR=$(intstren)")
            fig = figure()
        end

        plot_spectrum(intstrens,nrgs,idx,params_dict["nev"],"Interaction Strength",true; plot_title="")

        if idx == length(intstrens)
            fig = figure()
            occs = get_occupancy(states[1],lattice_params; if_plot=true,plot_title="$(lx)x$(ly) N=$n ULR=$(intstren)")
        end
    end
end































"fin"
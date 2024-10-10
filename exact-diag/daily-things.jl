#####################################################
#=

This file contains any random functions written to do one-off tasks

Depends on:
    execute-ed.jl

=#
######################################################

include("execute-ed.jl")


# redo gamma/omega calcs for all files
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

# calculate hatsugai data for 6x5 n=3 ULR=100.0
if false
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    params_dict = Dict([("Lx",6),("Ly",5),("N",3),("interaction_strength",100.0)])
    all_files = find_data_file(params_dict,"ed",dataloc; output_level=0)
    ref_multiplets,rm1_name,rm2_name = get_reference_multiplets(6,5,3; interaction_strength=100.0)
    for (i,f) in enumerate(all_files)
        if f == rm1_name || f == rm2_name
            continue
        end
        println(round(100*i/length(all_files),digits=3),"% done")
        println("Working on file: $f")
        filepath = dataloc * "/" * f
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        get_hatsugaifull(d["state"][1],d["state"][2],ref_multiplets; if_save=true,filepath=filepath,ref_multis_filenames=[rm1_name,rm2_name])
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
    nev = 13
    amplification_factor = 10.0

    cols = ["b","g","r","m","c"]
    if nev > length(cols)
        cols = repeat(cols,ceil(Int,nev/length(cols)))
    end

    intstrens = range(0.0,100.0,length=30)
    tws = range(0.0,1.0,length=50)
    for (idx,intstren) in enumerate(intstrens)
        #=where_fqh = 1
        if intstren > 1.5
            where_fqh += 2
        end
        if intstren > 2
            where_fqh += 2
        end=#
    #these_nrgs = zeros(Float64,nev,length(tws))
    #for (idx,tw1) in enumerate(tws)
    #for (idx2,tw2) in enumerate(tws)
        #intstren = 1.75
        
        tw2 = 0.5
        tw1 = 0.0
            params_dict = Dict([("Lx",3),("Ly",7),("N",3),("interaction_strength",intstren),("nev",nev),("if_save_data",false),("if_find_data",false),("if_periodic_x",true),("if_periodic_y",true),("tw1",tw1),("tw2",tw2)])
            states,nrgs,rhos,filepath,if_found,latpara,hamiltpara = run_normal_ed(params_dict; output_level=1)
            #get_occupancy(states[1],latpara; plot_title="ULR=$intstren tw2=$tw2 E=$(round(nrgs[1],digits=3))")
            #get_occupancy(states[2],latpara; plot_title="ULR=$intstren tw2=$tw2 E=$(round(nrgs[2],digits=3))")
            #=for i in 1:length(nrgs)
                scatter3D(tw1,intstren,nrgs[i] - nrgs[1],c=cols[i])
            end
            xlabel("Theta_x / 2pi")
            ylabel("ULR")=#
            plot_spectrum(intstrens,nrgs,idx,params_dict["nev"],"Interaction Strength"; plot_title=" $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"])")
            #plot_spectrum(tws,nrgs,idx,params_dict["nev"],"Theta_x / 2pi",false; plot_title=" $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"]) ULR=$intstren")
        end
        #title("Spectrum for Lx=3 Ly=7 N=3 tw2=$tw2")
    #end
end

if true
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




































"fin"
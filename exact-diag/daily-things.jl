#####################################################
#=

This file contains any random functions written to do one-off tasks

Depends on:
    execute-ed.jl

=#
######################################################

include("execute-ed.jl")
#include("plottings.jl")
#include("../other-funcs/basic-2d-plottings.jl")

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

# try charge pumping
if false
    lx,ly,n = 6,6,3
    nev = 10

    intstren = 0.0
    tws = [0.0,0.5]#range(0.0,1.0,length=21)

    pumped_charge = 0.0
    cppuls = zeros(Float64,length(tws))
    for (idx,tw1) in enumerate(tws)
        params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("nev",nev),("tw1",tw1),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("if_find_data",false),("if_save_data",false)])
        states,nrgs,rhos,filepath,if_found,latpara,hamiltpara = run_normal_ed(params_dict; output_level=1)

        occs = get_occupancy(states[1],latpara; if_plot=false)
        cppuls[idx] = sum([sum(occs[i,:] .* i) for i in 1:ly])
        #=if idx == 1
            global pumped_charge += cppul
        else
            global pumped_charge -= cppul
        end=#
        #plot_spectrum(tws,nrgs,idx,params_dict["nev"],L"\theta_x / 2 \pi",true; plot_title=" $(params_dict["Lx"])x$(params_dict["Ly"]) N=$(params_dict["N"]) ULR=$intstren")
    end
    #fig = figure()
    #scatter(tws,cppuls,c="b")
    #xlabel(L"\theta_x / 2 \pi")
    #ylabel("Charge Polarization per Unit Length")
    println("Pumped charge Q = $(round(cppuls[end] - cppuls[1],digits=3))")

end

# compare twisting for ED/1DeffMPS 4x4 n=2
if false
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


end

# hatsugai for 3x7 n=3 data collection
if false
    lx,ly,n = 3,7,3
    nev = 3
    #tws = range(-0.5,0.5,length=11)
    #tws2 = range(-0.75,-0.25,length=11)
    tws = range(0.0,1.0,length=10)
    tws2 = tws
    intstrens = range(0.0,2.0,length=11)

    for intstren in intstrens
        #lambda1s::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
        #lambda2s::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
        #omegas::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
        #ref_multis,rm1,rm2 = get_reference_multiplets(lx,ly,n; interaction_strength=intstren)
        for (idx,tw1) in enumerate(tws)
            for (idx2,tw2) in enumerate(tws2)
                params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("nev",nev),("tw2",tw2),("tw1",tw1),("interaction_strength",intstren),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("if_find_data",true),("if_save_data",true)])
                states,nrgs,rhos,filepath,if_found,latpara,hamiltpara = run_normal_ed(params_dict; output_level=1)

                #scatter3D(tw1,tw2,nrgs[1],c="b")
                #scatter3D(tw1,tw2,nrgs[2],c="g")
                #scatter3D(tw1,tw2,nrgs[3],c="r")

                #lambda1s[idx,idx2],lambda2s[idx,idx2],omegas[idx,idx2] = get_hatsugaifull(states[1],states[2],ref_multis; if_save=true,filepath=filepath,ref_multis_filenames=[rm1,rm2])
            end
        end
        #xlabel(L"\theta_x / 2 \pi")
        #ylabel(L"\theta_y / 2 \pi")
        #zlabel("Energy")
        #title("Energy Spectrum $(lx)x$(ly) N=$(n) ULR=$intstren")

        #plot_omega(tws,tws2,omegas; plot_title=" $(lx)x$(ly) N=$(n) ULR=$intstren")
        #plot_gamma(tws,tws2,lambda1s,1; plot_title=" $(lx)x$(ly) N=$(n) ULR=$intstren")
        #plot_gamma(tws,tws2,lambda2s,2; plot_title=" $(lx)x$(ly) N=$(n) ULR=$intstren")
    end

end

function get_closing_strength(intstrens,flatnesses::Vector{Float64})
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

# hatsugai data collection for (8,3,3)
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
if true
    ks = range(0.0,2*pi,length=100)
    intstren = 100.0
    lx,ly,n = 3,8,3
    params_dict = Dict([("Lx",lx),("Ly",ly),("N",n),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren)])
    dataloc = get_folder_location("cluster-data/exact-diag/torus")
    all_files = find_data_file(params_dict,"ed",dataloc; output_level=0)

    filter!(x -> !occursin("twist_angle1",x),all_files)
    filter!(x -> !occursin("mk",x),all_files)

    f = all_files[1]
    filepath = dataloc * "/" * f
    d,m = read_data_jld2(filepath; output_level=0)

    latparas = get_lattice_params_from_metadata(m)
    denscorrs = m["dens_corr_mat"]

    ftdds = zeros(Float64,length(ks),length(ks))
    for (idx,kx) in enumerate(ks)
        for (idx2,ky) in enumerate(ks)
            ftdd = ft_densitydensity_correlation([kx,ky],nothing,latparas; denscorrs=denscorrs,if_save=false,filepath=filepath)
            ftdds[idx,idx2] = abs(ftdd)
        end
    end
    fig = figure()
    imshow(reverse(ftdds,dims=1),extent=(0,2*pi,0,2*pi),origin="lower")
    xlabel(L"k_x")
    ylabel(L"k_y")
    colorbar()
    title("FT-DD for $(lx)x$(ly) N=$(n) ULR=$intstren")
end












































"fin"
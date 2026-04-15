#####################################################
#=

This file contains daily things for synth-dims TTNs

Depends on:
    synth-dims/long-range-ttn.jl
    review-practice-codes/observables.jl
    review-practice-codes/plottings.jl
    synth-dims/hatsugai-mbcn.jl

=#
######################################################

include("../other-funcs/include-other-files.jl")
#include_other_files(["synth-dims/long-range-ttn.jl","review-practice-codes/observables.jl","synth-dims/hatsugai-mbcn.jl","other-funcs/basic-2d-observables.jl"])
#include_other_files(["review-practice-codes/plottings.jl","other-funcs/basic-2d-plottings.jl"])
#include_other_files(["synth-dims/oneD-effective-LR.jl","synth-dims/plottings-oneD.jl"])

function datacollection_flatness_1deff(Lx::Int64,Ly::Int64,N::Int64; kwargs...)
    hanis::Float64 = get(kwargs,:hopping_anisotropy,1.0)
    if_hatsugai::Bool = get(kwargs,:if_hatsugai,true)
    cutoff::Float64 = get(kwargs,:cutoff,1e-8)
    nrgtol::Float64 = get(kwargs,:nrgtol,1e-6)
    mdim::Int64 = get(kwargs,:mdim,200)


    tws_count::Int64 = get(kwargs,:tws_count,10)
    tws_start::Float64 = get(kwargs,:tws_start,0.0)
    tws_end::Float64 = get(kwargs,:tws_end,1.0)
    tws::Vector{Float64} = range(tws_start,tws_end,length=tws_count)

    println("Starting Flatness Data Collection for 1D effective $(Lx)x$(Ly) N=$(N) from Twists $(tws_start) - $(tws_end)")
    sleep(1.0)

    if if_hatsugai
        ref_multis,rm1,rm2 = get_reference_multiplets(Lx,Ly,N; hopping_anisotropy=hanis,if_make_new=true)
    end
    for (idx2,tw1) in enumerate(tws)
        for (idx3,tw2) in enumerate(tws)

            params_dict = Dict([("Lphys",Lx),("Lsynth",Ly),("particles",N),("tw2",tw2),("tw1",tw1),("if_remapping",false),("es_count",2),("nrgtol",nrgtol),("cutoff",cutoff),("mdim",mdim),("if_periodic_phys",true),("if_periodic_synth",true),("filling",0.5),("if_find_data",true),("if_save_data",true)])
            psis,rhos,nrgs,model_paras,if_found = run_normal_1deffmps(params_dict)
            filepath = model_paras[:location]*"/"*model_paras[:name]*".jld2"

            if isnothing(psis[1])
                continue
            end

            if if_hatsugai
                if if_found
                    d,m = read_data_jld2(filepath; output_level=0)
                    if !haskey(m,"omega")
                        lambda1,lambda2,omega = get_hatsugaifull(psis[1],psis[2],ref_multis; if_save=true,filepath=filepath,ref_multis_filenames=[rm1,rm2])
                    end
                else
                    lambda1,lambda2,omega = get_hatsugaifull(psis[1],psis[2],ref_multis; if_save=true,filepath=filepath,ref_multis_filenames=[rm1,rm2])
                end
            end
        end
    end
end

#= testing factory run on TTN with new gauge with seed ttn
if false
    lx,ly,n = 16,8,8
    layers = Int(log(2,lx*ly))
    intstren = 300.0

    which_one = 1
    dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    previous_pdict = Dict([("layers",layers),("onsite_strength",intstren),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files_previous = find_data_file(previous_pdict,"wavefuncttn",dataloc)
    display(all_files_previous)
    f_previous = all_files_previous[which_one]
    d_previous,m_previous = read_data(joinpath(dataloc,f_previous); output_level=0)
    previous_ttn = d_previous["ttn"]

    params_dict = Dict([("hopping_anisotropy",1.0),("expander_fraction",1e-5),("seed_ttn",previous_ttn),("if_redo",true),("particles",n),("layers",layers),("mdim",400),("if_save_data",true),("filling",0.5),("onsite_strength",intstren),("lr","all"),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_results = run_synth_dims_generic(params_dict)

end=#

#= rerun 16x8 laughlin to converge further
if false
    BLAS.set_num_threads(5)

    lx,ly,n = 12,6,6
    #layers = Int(log(2,lx*ly))
    intstren = 300.0
    dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge/pinned-scaling")
    pdict = Dict([("Lx",lx),("Ly",ly),("onsite_strength",intstren),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    display(all_files)

    newparas = Dict([("nrgtol",1e-7),("mdim",500)])
    rez = reconverge_ttn(joinpath(dataloc,all_files[1]); new_parameters=newparas)
end=#

#= data collection of 4pt MPO
if false
    lx,ly,n = 16,4,8
    #layers = Int(log(2,lx*ly))
    #intstren = 0.0

    dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict = Dict([("particles",n),("Lx",lx),("Ly",ly),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    display(all_files)
    #f = all_files[1]
    for f in all_files
        d,m = read_data(dataloc * "/wavefunc" * f; output_level=0)

        #if "fourpt_momentum" in keys(m)
        #    continue
        #end

        psi1 = d["ttn"]
        psi2 = d["ttn_1"]

        #fourpt_mpo = four_point([psi1,psi2]; output_level=0)
        fourpt_mpo = four_point(psi1; output_level=1)
        fourpt_mpo_2 = four_point(psi2; output_level=1)

        datadict = Dict([("fourpt_momentum",fourpt_mpo),("fourpt_momentum_1",fourpt_mpo_2)])
        #datadict = Dict([("fourpt_momentum",fourpt_mpo)])
        modify_data(datadict,dataloc * "/" * f,"metadata"; output_level=0)
    end
end=#

#= calculate entanglement spectrum for new data
if false
    BLAS.set_num_threads(5)
    dataloc = get_folder_location("cluster-data/synth-dims/torus/")
    pdict = Dict([("layers",5),("particles",4),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc) 
    display(all_files)
    
    for f in all_files
        d,m = read_data(dataloc * "/wavefunc" * f; output_level=0)
        
        fulldata = zeros(ComplexF64,pdict["layers"]-1,maxlinkdim(d["ttn"]))
        for k in 0:pdict["layers"] - 2
            entspec = spatial_entanglement_spectrum(d["ttn"]; layers_down=k)
            fulldata[k+1,:] = vcat(entspec,zeros(ComplexF64,maxlinkdim(d["ttn"]) - length(entspec)))
        end
        datadict = Dict([("entanglement_spectrum",fulldata)])
        modify_data(datadict,dataloc * "/" * f,"metadata"; output_level=0)
    end


end=#

# plot entanglement entropy to look for kDW phase transition
#= no reason to suspect transition visible from here
if false
    dataloc = get_folder_location("cluster-data/synth-dims/torus/")
    pdict = Dict([("layers",5),("particles",4),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc) 
    display(all_files)
    
    entropies = Dict([("1",[]),("2",[]),("3",[]),("4",[])])
    intstrens = []
    for f in all_files
        d,m = read_data(joinpath(dataloc,f); output_level=0)
        entanglement_spectrum = real.(m["entanglement_spectrum"])
        append!(intstrens,[m["onsite_strength"]])
        for i in 1:4
            local_ent_spec = filter(x -> x != 0.0,entanglement_spectrum[i,:])
            println("For ULR $(m["onsite_strength"]) and layer $(i) the smallest eigenvalue is ",minimum(local_ent_spec)," and has values ",length(local_ent_spec))
            #display(local_ent_spec)
            append!(entropies[string(i)],[entanglement_entropy(local_ent_spec)])
        end
    end
    for i in 1:4
        scatter(intstrens,entropies[string(i)],label="$(i)")
    end
    xlabel("Interaction Strength")
    ylabel("Entanglement Entropy")
    legend()
end=#

#= check if expander is so required 8x4 and 16x8
if false
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))
    intstren = 0.0


    params_dict = Dict([("hopping_anisotropy",1.0),("es_count",1),("particles",n),("layers",layers),("mdim",200),("if_save_data",true),("filling",0.5),("onsite_strength",intstren),("lr","all"),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_results = run_synth_dims_generic(params_dict)

    #=dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    filter!(x -> !occursin("gpu",x),all_files)
    display(all_files)
    f = all_files[1]
    d,m = read_data(dataloc * "/wavefunc" * f; output_level=0)
    psi_withexp = d["ttn"]
    #psi2_withexp = d["ttn_1"]=#

    #=pdict_noexp = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files_noexp = find_data_file(pdict_noexp,"ttn",dataloc)
    filter!(x -> occursin("gpu",x),all_files_noexp)
    display(all_files_noexp)
    f_noexp = all_files_noexp[1]
    d_noexp,m_noexp = read_data(dataloc * "/" * f_noexp; output_level=0)
    psi_noexp = d_noexp["ttn"]
    #psi2_noexp = d_noexp["ttn_1"]=#

    # check overlap of states



end=#

#= plot 4pt 16x8
if false
     
    lx,ly,n = 16,8,8
    layers = Int(log(2,lx*ly))
    #intstren = 0.0

    dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict = Dict([("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)

    f = all_files[5]
    d,m = read_data(joinpath(dataloc,f))

    fourpt_vals = m["fourpt_momentum"]
    #plot_four_point(fourpt_vals; plot_title="GS1 16x8 N=8 ULR=$(m["onsite_strength"])")

    #fourpt_vals2 = m["fourpt_momentum_1"]
    fourpt_vals2 = zeros(Float64,size(fourpt_vals))
    for i in 1:size(fourpt_vals,1)-1
        for j in 1:size(fourpt_vals,2)-1
            fourpt_vals2[i,j] = fourpt_vals[i+1,j+1]
        end
    end
    fourpt_vals2[end,:] = vcat(fourpt_vals[1,:][2:end],fourpt_vals[1,:][1])
    fourpt_vals2[:,end] = vcat(fourpt_vals[:,1][2:end],fourpt_vals[:,1][1])
    #plot_four_point(fourpt_vals2; plot_title="GS2 16x8 N=8 ULR=$(m["onsite_strength"])")

    mixed_fourpt = 0.5 .* (fourpt_vals .+ fourpt_vals2)
    #plot_four_point(mixed_fourpt; plot_title="GS1+GS2 16x8 N=8 ULR=$(m["onsite_strength"])")
    

end=#

#= compare 4pt for 8x4 and 16x8 for finite size scaling
if false
    
    big_lx,big_ly,big_n = 16,8,8
    big_layers = Int(log(2,lx*ly))
    intstren = 0.0

    dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    big_pdict = Dict([("layers",big_layers),("onsite_strength",intstren),("particles",big_n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    big_all_files = find_data_file(big_pdict,"ttn",dataloc)

    big_f = big_all_files[1]    
    big_d,big_m = read_data(joinpath(dataloc,big_f))

    big_fourpt_vals = 0.5 .* (big_m["fourpt_momentum"] .+ big_m["fourpt_momentum_1"])
    plot_four_point(big_fourpt_vals; plot_title="GS1 16x8 N=8 ULR=$(big_m["onsite_strength"])")

    small_lx,small_ly,small_n = 8,4,4
    small_layers = Int(log(2,small_lx*small_ly))
    small_pdict = Dict([("layers",small_layers),("onsite_strength",intstren),("particles",small_n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    small_all_files = find_data_file(small_pdict,"ttn",dataloc)
    small_f = small_all_files[1]
    small_d,small_m = read_data(joinpath(dataloc,small_f))
    small_fourpt_vals = 0.5 .* (small_m["fourpt_momentum"] .+ small_m["fourpt_momentum_1"])
    plot_four_point(small_fourpt_vals; plot_title="GS1 8x4 N=4 ULR=$(small_m["onsite_strength"])")

end=#

#= look at finite size scaling of commensurate filling interaction strength spectrum
if false
    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    pdict = Dict([("layers",6),("particles",8),("hopping_anisotropy",1.0),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    all_files = vcat(all_files[1],all_files[3:end])

    chosen_intstrens = collect(range(0.0,2.0,length=11))[1:5]

    nrgs = Dict([("1",[]),("2",[]),("3",[])])
    intstrens = []
    for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        if_done,all_checks = check_nrg_convergence(m)
        #=if !(m["onsite_strength"] in chosen_intstrens)
            continue
        end=#
        append!(intstrens,[m["onsite_strength"]])
        for i in 1:3
            keyname = i == 1 ? "observer" : "observer_$(i-1)"
            if keyname in keys(m) && all_checks[string(i-1)]
                append!(nrgs[string(i)],[m[keyname].nrg[end]])
            else
                append!(nrgs[string(i)],[0.0])
            end
        end
        #get_occupancy(d["densmat"]; plot_title="Intstren = $(m["onsite_strength"])",fix_colorbar=true)
    end

    fig = figure()
    for i in 1:3
        scatter(intstrens,nrgs[string(i)] .- nrgs["1"],label="E$(i-1)")
    end
    xlabel("Interaction Strength")
    ylabel("Energy Difference")
    legend()
    ylim([-0.1,0.7])
end

# look at CDW structure factor
if false
    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    pdict = Dict([("layers",6),("particles",8),("hopping_anisotropy",1.0),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    all_files = vcat(all_files[1],all_files[3:end])

    d,m = read_data_jld2(dataloc * "/" * all_files[1]; output_level=0)

    denscorrs = make_density_correlations(d["ttn"])
    range_cdwsf_angles(100,denscorrs; if_plot=true,plot_title="V = $(m["onsite_strength"])")
end

# 1Deff hatsugai
if false
    lx,ly,n = 4,4,2

    #tws2 = range(0.1,0.3,length=11)
    #tws = range(0.7,0.9,length=11)
    tw1 = 0.33
    tw2 = 0.67

    #lambda1s::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
    #lambda2s::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
    #omegas::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
    #ref_multis,rm1,rm2 = get_reference_multiplets(lx,ly,n)
    #for (idx,tw1) in enumerate(tws)
        #for (idx2,tw2) in enumerate(tws2)
            params_dict = Dict([("Lphys",lx),("Lsynth",ly),("particles",n),("tw2",tw2),("tw1",tw1),("if_remapping",false),("es_count",2),("nrgtol",1e-6),("mdim",200),("if_periodic_phys",true),("if_periodic_synth",true),("filling",0.5),("if_find_data",false),("if_save_data",false)])
            psis,rhos,nrgs,model_paras,if_found_data = run_normal_1deffmps(params_dict)

            params_dict["tw1"] = tw2
            params_dict["tw2"] = tw1
            psis2,rhos2,nrgs2,model_paras2,if_found_data2 = run_normal_1deffmps(params_dict)
            
            #scatter3D(tw1,tw2,nrgs[1],c="b")
            #scatter3D(tw1,tw2,nrgs[2],c="g")
            #scatter3D(tw1,tw2,nrgs[3],c="r")

            #=if if_found_data
                d,m = read_data_jld2(model_paras[:location]*"/"*model_paras[:name]*".jld2"; output_level=0)
                lambda1s[idx,idx2] = m["gamma1"]
                lambda2s[idx,idx2] = m["gamma2"]
                omegas[idx,idx2] = m["omega"]
            else
                lambda1s[idx,idx2], lambda2s[idx,idx2], omegas[idx,idx2] = get_hatsugaifull(psis[1],psis[2],ref_multis; if_save=true,filepath=model_paras[:location]*"/"*model_paras[:name],ref_multis_filenames=[rm1,rm2])
            end
        end
    end=#
    #xlabel(L"\theta_x / 2\pi")
    #ylabel(L"\theta_y / 2\pi")
    #zlabel("Energy")

    dataloc = get_folder_location("cluster-data/synth-dims/twists")
    pdict = Dict([("Lphys",lx),("Lsynth",ly),("nbosons",n),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_files = find_data_file(pdict,"mps",dataloc)

    #tw1s = sort(unique(vcat(range(0.1,0.3,length=11),range(0.7,0.9,length=11),range(0.0,1.0,length=21))))
    #tw2s = sort(unique(vcat(range(0.7,0.9,length=11),range(0.1,0.3,length=11),range(0.0,1.0,length=21))))
    tw1s = range(0.0,1.0,length=21)
    tw2s = range(0.0,1.0,length=21)



    lambda1s::Matrix{ComplexF64} = zeros(ComplexF64,length(tw1s),length(tw2s))
    lambda2s::Matrix{ComplexF64} = zeros(ComplexF64,length(tw1s),length(tw2s))
    #omegas::Matrix{ComplexF64} = zeros(ComplexF64,length(tw1s),length(tw2s))
    #nrgs::Dict{String,Matrix{Float64}} = Dict([("1",zeros(Float64,length(tw1s),length(tw2s))),("2",zeros(Float64,length(tw1s),length(tw2s))),("3",zeros(Float64,length(tw1s),length(tw2s)))])
    for f in all_files
        filename_dict = get_params_from_filename(f)
        d,m = read_data_jld2(dataloc * "/wavefunc" * f; output_level=0)
        try
            tw1,tw2 = filename_dict["tw1"],filename_dict["tw2"]
        catch
            continue
        end
        if tw1 in [0.33,0.67] || tw2 in [0.33,0.67]
            continue
        end
        if !(tw1 in tw1s) || !(tw2 in tw2s)
            continue
        end
        println("Twist1 is ",tw1," and Twist2 is ",tw2)
        idx = findfirst(x->tw1s[x]==tw1,1:length(tw1s))
        idx2 = findfirst(x->tw2s[x]==tw2,1:length(tw2s))

        lambda1s[idx,idx2] = get_gamma(d["mps"],d["mps_1"],psis[1:2])
        lambda2s[idx,idx2] = get_gamma(d["mps"],d["mps_1"],psis2[1:2])

        #lambda1s[idx,idx2] = m["gamma1"]
        #lambda2s[idx,idx2] = m["gamma2"]
        #omegas[idx,idx2] = m["omega"]
        #=nrgs["1"][idx,idx2] = m["observer"].energies[end]
        for i in 2:3
            nrgs[string(i)][idx,idx2] = m["observer_$(i-1)"].energies[end]
        end=#
    end

    plot_omega(tw1s,tw2s,omegas; if_perfect_grid=true)
    plot_gamma(tw1s,tw2s,lambda1s,1)
    plot_gamma(tw1s,tw2s,lambda2s,2)
    
    #=cs = ["b","g","r"]
    fig = figure()
    for i in 1:3
        for ii in 1:length(tw1s)
            for jj in 1:length(tw2s)
                scatter3D(tw1s[ii],tw2s[jj],nrgs[string(i)][ii,jj],c=cs[i])
            end
        end
    end=#
end

# making ft_dd data
if false
    configs = [(8,3,3),(4,6,3),(3,8,3)]#,(8,4,4),(8,5,5)]
    for (lx,ly,n) in configs
    #lx,ly,n = 8,3,3
        pdict = Dict([("Lphys",lx),("Lsynth",ly),("nbosons",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
        dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
        all_files = find_data_file(pdict,"mps",dataloc; output_level=0)
        
        filter!(x -> !occursin("twist_angle1",x),all_files)
        filter!(x -> !occursin("mk",x),all_files)

        display(all_files)

        filepath = dataloc * "/" * all_files[1]

        d,m = read_data_jld2(filepath; output_level=0)

        local_nrgs = get_local_nrgs(m)
        if local_nrgs[1] != sort(local_nrgs)[1]
            println("Levels are not sorted")
            which_level_gs = findfirst(x -> minimum(local_nrgs) == local_nrgs[x],1:length(local_nrgs))
            wavefunc_string = "mps_$(which_level_gs-1)"
        else
            which_level_gs = 1
            wavefunc_string = "mps"
        end

        d_wavefunc = read_data_jld2(dataloc * "/wavefunc" * all_files[1]; output_level=0)

        results = zeros(Float64,2)
        results[1] = abs(ft_densitydensity_correlation(pi/2,d_wavefunc[wavefunc_string]; if_save=true,filepath=filepath))
        results[2] = abs(ft_densitydensity_correlation(0.0,d_wavefunc[wavefunc_string]; if_save=true,filepath=filepath))

        println("Final FT-DD value for $(lx)x$(ly) n=$n is ",results)
    end
end

# save denscorrs for all configs
if false
    configs = [(8,3,3),(4,6,3),(3,8,3),(8,4,4),(8,5,5)]
    for (lx,ly,n) in configs
        pdict = Dict([("Lphys",lx),("Lsynth",ly),("nbosons",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
        dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
        all_files = find_data_file(pdict,"mps",dataloc; output_level=0)
        
        filter!(x -> !occursin("twist_angle1",x),all_files)
        filter!(x -> !occursin("mk",x),all_files)

        f = all_files[1]

        filepath_metadata = dataloc * "/" * f
        filepath_wavefunc = dataloc * "/wavefunc" * f
        d,m = read_data_jld2(filepath_metadata; output_level=0)

        if !haskey(m,"dens_corr_mat")
            d_wavefunc = read_data_jld2(filepath_wavefunc; output_level=0)
            make_density_correlations(d_wavefunc[get_correct_gs_mpsstring(m)]; if_save=true,filepath=filepath_metadata)
        end
        println("Finished $(lx)x$(ly) N=$(n)")

    end
end

# save ft-dd at given angle for all configs
if false
    configs = [(8,3,3),(4,6,3),(3,8,3),(8,4,4),(8,5,5)]
    for (lx,ly,n) in configs
        pdict = Dict([("Lphys",lx),("Lsynth",ly),("nbosons",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
        dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
        all_files = find_data_file(pdict,"mps",dataloc; output_level=0)
        
        filter!(x -> !occursin("twist_angle1",x),all_files)
        filter!(x -> !occursin("mk",x),all_files)

        f = all_files[1]

        filepath_metadata = dataloc * "/" * f
        d,m = read_data_jld2(filepath_metadata; output_level=0)

        #if !haskey(m,"dens_corr_mat")
            denscorrs = m["dens_corr_mat"]
            ftdd_val_0p5 = ft_densitydensity_correlation(pi/2,nothing; denscorrs=denscorrs,if_save=true,filepath=filepath_metadata)
            ftdd_val_0p0 = ft_densitydensity_correlation(0.0,nothing; denscorrs=denscorrs,if_save=true,filepath=filepath_metadata)
        #end
        println("Finished $(lx)x$(ly) N=$(n)")

        scatter(n/lx,abs(ftdd_val_0p5),c="b")
        scatter(n/lx,abs(ftdd_val_0p0),c="r")
        xlabel(L"\rho_{1D}")
        ylabel("FT-DD Value")
        title("FT-DD at 0pi and pi/2 ranging rho_1D")

    end

end

# look at full angle range FT-DD for given configuration
if false
    lx,ly,n = 3,8,3
    pdict = Dict([("Lphys",lx),("Lsynth",ly),("nbosons",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])

    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    all_files = find_data_file(pdict,"mps",dataloc; output_level=0)
    
    filter!(x -> !occursin("twist_angle1",x),all_files)
    filter!(x -> !occursin("mk",x),all_files)

    f = all_files[1]

    filepath_metadata = dataloc * "/" * f
    d,m = read_data_jld2(filepath_metadata; output_level=0)

    denscorrs = m["dens_corr_mat"]

    momenta = range(0.0,2*pi,length=100)
    for k in momenta
        ftdd_val = ft_densitydensity_correlation(k,nothing; denscorrs=denscorrs)
        scatter(k / pi,abs(ftdd_val),c="b")
    end
    xlabel("Momenta / pi")
    ylabel("FT-DD Value")
    title("FT-DD for $(lx)x$(ly) N=$(n)")

end

# remove all seed_ttn data and set to nothing
if false
    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    all_files = readdir(dataloc)
    filter!(x -> occursin("ttn",x),all_files)
    println("Found $(length(all_files)) files")
    seed_count = 0
    for f in all_files
        bf = jldopen(dataloc * "/" * f,"a+")
        m = bf["metadata"]
        if haskey(m,"seed_ttn") && !isnothing(m["seed_ttn"])
            global seed_count += 1
            close(bf)
            modify_data_jld2(Dict([("seed_ttn",nothing)]),dataloc * "/" * f,"metadata"; output_level=0)
            println("Removed from $f")
        else
            close(bf)
        end
    end
    println("Of the $(length(all_files)) files, $(seed_count) had seed_ttn data")
end

# density-density stuff for 16x8
if false
    dataloc = get_folder_location("cluster-data/synth-dims/torus")
    pdict = Dict([("layers",7),("particles",8),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    display(all_files)

    #ulrs = []
    #ft_dds = []
    for f in all_files
    #f = all_files[1]
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

        if m["onsite_strength"] in [5.0,50.0,1.0]
            continue
        end


        if !haskey(m,"densitydensity")
            continue
        end

        #= checking convergence of energies
        all_convs = [false,false,false]
        for i in 1:3
            observerkey = i == 1 ? "observer" : "observer_$(i-1)"
            if observerkey in keys(m)
                all_convs[i] = abs(m[observerkey].nrg[end] - m[observerkey].nrg[end-1]) < m["nrgtol"]
            end
        end
        println("For $(m["onsite_strength"]) all convs are ",all_convs)=#

        occs = get_occupancy(d["densmat"]; plot_title="ULR = $(m["onsite_strength"])",if_plot=false,if_synth_rectangle=true)
        dds = m["densitydensity"]
        #ft_dds_stripe = abs(ft_densitydensity([pi,0],dds))


        
        pairdist = pairdistribution(dds,occs; if_plot=false, plot_title="ULR = $(m["onsite_strength"])",vmax=1.5)
        #=centersite = [Int64(ceil(size(pairdist,2)/2)),Int64(ceil(size(pairdist,1)/2))]
        #fig = figure()
        #plot(1:size(pairdist,2),pairdist[centersite[2],:],"-p",label=m["onsite_strength"])
        #xlabel("Physical Site")
        plot(1:size(pairdist,1),pairdist[:,centersite[1]],"-p",label=m["onsite_strength"])
        xlabel("Synthetic Site")
        ylabel("Pair Dist")
        title("Pair Distribution Synthetic Slice 16x8 N=8")
        legend()=#

        rez = pairdist_ellipticalness(pairdist)

        #=scatter(rez[2],rez[3],c="b")
        xlabel("Center X")
        ylabel("Center Y")
        xlim([-1.0,1.0])
        ylim([-1.0,1.0])=#


        #=println("For ULR = $(m["onsite_strength"]) the x_var is $(round(rez[4],digits=4)) and y_var is $(round(rez[5],digits=4))")
        scatter(m["onsite_strength"],rez[1],c="b")
        xlabel("Interaction Strength")
        ylabel("Variance Ratio, "*L"\Delta x / \Delta y")
        title("Ellipticalness of Pair Distribution 16x8 N=8")
        xscale("log")=#

        #plot_fourpointcorrelator(dds; if_plot=true,plot_title="ULR=$(m["onsite_strength"])")
        #=ks = range(-pi,pi,length=100)
        allmoms::Matrix{Float64} = zeros(Float64,length(ks),length(ks))
        for (idx1,kx) in enumerate(ks)
            for (idx2,ky) in enumerate(ks)
                allmoms[idx2,idx1] = abs(ft_densitydensity([kx,ky],dds))
            end
        end
        fig = figure()
        imshow(allmoms,extent=(minimum(ks),maximum(ks),minimum(ks),maximum(ks)),origin="lower")
        xlabel(L"k_{phys}")
        ylabel(L"k_{synth}")
        colorbar()
        title("FT-Density-Density for 16x8 N=8 ULR=$(m["onsite_strength"])")=#
        #append!(ulrs,[m["onsite_strength"]])
        #append!(ft_dds,[ft_dds_stripe])


    end
    #=fig = figure()
    scatter(ulrs,ft_dds)
    xlabel("Interaction Strength")
    #ylabel("FT-DD at k=("*L"\pi"*",0)")
    title("FT-DD at k=("*L"\pi"*",0)")=#
end

# trying DD for 12x6
if false
    dataloc = get_folder_location("cluster-data/synth-dims/torus")
    pdict = Dict([("layers",7),("particles",6),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)

    f = all_files[1]
    d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
    d_wavefunc,m_wavefunc = read_data_jld2(dataloc * "/wavefunc" * f; output_level=0)

    occs = get_occupancy(d["densmat"]; plot_title="Intstren = $(m["onsite_strength"])",if_plot=true,if_synth_rectangle=true)
    dds = fourpoint_alberto(d_wavefunc["ttn"]; if_plot=true,plot_title="Intstren = $(m["onsite_strength"])")

end

# doing scaling of bipartition size
if false
    dataloc = get_folder_location("cluster-data/synth-dims/torus")
    pdict = Dict([("layers",7),("particles",8),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc) 
    display(all_files)
    
    #=
    f = all_files[1]
    d = read_data_jld2(dataloc * "/" * f; output_level=0)

    ee = zeros(Float64,pdict["layers"] - 1)
    perims = [16,12,8,6]
    for k in 0:pdict["layers"] - 2
        entspec = spatial_entanglement_spectrum(d["ttn"]; layers_down=k)
        ee[k+1] = entanglement_entropy(entspec)
        println("Local Entanglement Entropy at Layer $(pdict["layers"]-k) with bonddim $(length(entspec)) is ",ee[k+1])
    end
    scatter(perims,ee,c="b")
    xlabel("Perimeter")
    ylabel("Entanglement Entropy")=#

    all_ees = []
    perims = []
    intstrens = [0.0,1.0,100.0]
    markers = ["o","s","^"]
    global i = 0
    for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        if !(m["onsite_strength"] in intstrens)
            continue
        end
        global i += 1
        ee = m["ee_scaling"]
        ee_vec = Float64[]
        perims_vec = Int64[]
        for (k,v) in ee
            append!(ee_vec,v)
            append!(perims_vec,calculate_perimeter(k))
        end
        append!(all_ees,[ee_vec])
        global perims = perims_vec
        scatter(perims_vec[2:end],ee_vec[2:end],marker=markers[i],label=m["onsite_strength"])
    end
    xlabel("Perimeter")
    ylabel("Entanglement Entropy")
    legend()
    xlim([0,26])
    ylim([-1,5])

    first_fit = linear_fit(perims[2:end],all_ees[1,:][1][2:end])
    plot(range(0,35,length=2),first_fit[1] .+ (first_fit[2] .* range(0,35,length=2)))
    title("Entanglement Entropy Scaling 16x8 N=8, Top EE = $(round(-first_fit[1],digits=3))")


    #=for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        entspec = m["entanglement_spectrum"]
        ee = entanglement_entropy(entspec)
        scatter(m["onsite_strength"],ee,c="b")
    end=#
end

# abcd calculation of TEE from cluster data
if false
    dataloc = get_folder_location("cluster-data/synth-dims/torus")
    pdict = Dict([("layers",7),("particles",8),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc) 
    display(all_files)

    cols = Dict([("0.0","b"),("1.0","g"),("100.0","r")])
    for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

        if !(m["onsite_strength"] in [0.0,1.0,100.0])
            continue
        end

        as = [0,0,0,2,0,4]
        for i in [4,6]
            if i == 4
                scatter(as[i],m["tee_$i"],c=cols[string(m["onsite_strength"])],label=m["onsite_strength"])
            else
                scatter(as[i],m["tee_$i"],c=cols[string(m["onsite_strength"])])
            end
        end
        legend()
        xlabel("Bipartition Width, a")
        ylabel("Topological Entanglement Entropy")

    end
end

# time scaling for long range interaction strength
if false
    dataloc = get_folder_location("cluster-data/synth-dims/torus")
    pdict = Dict([("layers",7),("particles",8),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc) 
    display(all_files)    

    for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

        intstren = []
        times = []
        if !haskey(m,"runtime")
            continue
        else
            append!(intstren,m["onsite_strength"])
            append!(times,m["runtime"])
        end
        scatter(intstren,times)
        xlabel("Interaction Strength")
        ylabel("Runtime")
        yscale("log")

    end

end
=#

#= spatial entanglement spectrum for 16x8
if false
    lx,ly,n = 16,8,8
    dataloc = get_folder_location("cluster-data/synth-dims/torus")
    pdict = Dict([("layers",7),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc) 

    for (idx,f) in enumerate(all_files)
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)

        #plot_spectrum(intstrens,m["entanglement_spectrum"][1:20],idx,20,"Interaction Strength",false)
        
        if !haskey(m,"entanglement_spectrum")
            d_wavefunc,m_wavefunc = read_data(dataloc * "/wavefunc" * f; output_level=0)
            entspec = spatial_entanglement_spectrum(d_wavefunc["ttn"]; if_save=true,filepath=dataloc * "/" * f)
        end

    end

end=#

#= test Noah's GPU implementation for accuracy comparing to CPU
if false
    lx,ly,n = 4,4,2
    intstren = 0.0
    #=layers = Int(log(2,lx*ly))
    dataloc = get_folder_location("cluster-data/synth-dims/torus/new-gauge")
    pdict = Dict([("layers",layers),("particles",n),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    filter!(x -> occursin("if_gpu-true",x),all_files)
    f = all_files[1]
    d_gpu,m_gpu = read_data(joinpath(dataloc,f); output_level=0)
    d_wavefunc_gpu,m_wavefunc_gpu = read_data(joinpath(dataloc,"wavefunc"*f); output_level=0)=#


    params_dict = Dict([("if_gpu",false),("num_sweeps",5),("nrgtol",1e-6),("lr","all"),("hopping_anisotropy",1.0),("Lx",lx),("Ly",ly),("es_count",1),("expander_fraction",0.5),("particles",n),("mdim",100),("if_save_data",false),("filling",0.5),("if_find_data",false),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true)])
    all_states, hamilt, all_obs, all_densmats, all_runtimes = run_synth_dims_generic(params_dict)

    
end=#

#= test NSight profiling
using CUDA
if false
    lx,ly,n = 8,4,4
    intstren = 10.0

    params_dict = Dict([("if_gpu",true),("nrgtol",5e-5),("lr","all"),("hopping_anisotropy",1.0),("Lx",lx),("Ly",ly),("es_count",0),("expander_fraction",1e-5),("particles",n),("mdim",300),("if_save_data",false),("filling",0.5),("if_find_data",false),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true)])
    CUDA.@allowscalar all_states, hamilt, all_obs, all_densmats, all_runtimes = run_synth_dims_generic(params_dict)

end=#

function get_inter_coeff(s1,s2,t_strength,phi,edge_length_x,edge_length_y; kwargs...) 
	hopping_anisotropy = get(kwargs, :hopping_anisotropy, 1.0)
	#t_strength_phys = t_strength * hopping_anisotropy
	flux_direction = get(kwargs,:flux_direction,"phys")
	

	t_strength_phys,t_strength_synth = 1.0,1.0#get_hopping_strengths(t_strength,hopping_anisotropy)
	
	if get(kwargs, :no_magF, false)
		phi = 0.0
	end
	
	if s1[1] == s2[1] # Synthetic Dimension Hopping

		stren = -t_strength_synth
		flux_direction == "synth" ? stren *= exp(im*2*pi*(phi*(s1[1]-1))) : nothing
		return stren
	elseif s1[2] == s2[2] # Physical Dimension Hopping

		stren = -t_strength_phys
		flux_direction == "phys" ? stren *= exp(im*2*pi*(phi*(s1[2]-1))) : nothing
		return stren
	else
		return 0.0
	end

end

# testing expander for Trento
if true
    lx,ly,n = 4,4,8
    intstren = 0.0
    num_layers = Int(log(2,lx*ly))
    
    params_dict = Dict([("if_gpu",false),("if_check_fluxes",false),("outputlevel",1),("nrgtol",5e-5),("lr","all"),("hopping_anisotropy",1.0),("Lx",lx),("Ly",ly),("es_count",0),("expander_fraction",0.9),("particles",n),("mdim",40),("if_save_data",false),("filling",0.5),("if_find_data",false),("onsite_strength",intstren),("if_periodic_phys",true),("if_periodic_synth",true)])

    net = TTN.BinaryRectangularNetwork(num_layers, TTN.ITensorNode, "Boson";conserve_qns=true,dim=2)
    println("Made Network")


    restricted_size = (lx,ly)
    phys_edge_length = lx
    virt_edge_length = ly
    which_axis = 1
    t_strength = 1.0
    phi = 2 * n / (lx*ly)
    if_per = [params_dict["if_periodic_phys"],params_dict["if_periodic_synth"]]
    twist_angle = [0.0,0.0]
    hopping = TTN.OpSum()
    for s_phys in 1:restricted_size[1]
        for s_synth in 1:restricted_size[2]
            starting_site = [s_phys,s_synth]
            twist = 0
            for which_axis in [1,2]
                    ending_site = starting_site .+ ((which_axis == 1,which_axis == 2))

                    # enforce boundary conditions
                    if ending_site[which_axis] > restricted_size[which_axis]
                        if if_per[which_axis]
                            ending_site[which_axis] = 1
                            twist = 1
                        else
                            continue
                        end
                    end

                    if ending_site[which_axis] < 1
                        if if_per[which_axis]
                            ending_site[which_axis] = restricted_size[which_axis]
                            twist = 2
                        else
                            continue
                        end
                    end

                    coeff = get_inter_coeff(starting_site,ending_site,t_strength,phi,phys_edge_length,virt_edge_length)
                    twist == 1 ? coeff *= exp(im*twist_angle[which_axis]*2*pi) : nothing
                    twist == 2 ? coeff *= exp(-im*twist_angle[which_axis]*2*pi) : nothing
                    global hopping += (coeff,"Adag",Tuple(starting_site),"A",Tuple(ending_site))
                    global hopping += (conj(coeff),"Adag",Tuple(ending_site),"A",Tuple(starting_site))
                    twist = 0
            end
        end
    end
    ham_tpo = TTN.TPO(hopping,physical_lattice(net))
    println("Made Hamiltonian")

    #states = fill_states(n,lx*ly,1)
	#old_ttn = TTN.ProductTreeTensorNetwork(net,states)
    old_ttn = TTN.initialize_ttn(TTN.ProductTreeTensorNetwork(net,fill("0", lx*ly)),20,n; part_type="Boson")
    println("Made Initial TTN")

    observer = NRGVarObserver(params_dict["nrgtol"])
    expander = DefaultExpander(params_dict["expander_fraction"])
    println("Running DMRG")
    sp = TTN.dmrg(old_ttn,ham_tpo; expander=expander, number_of_sweeps=100, maxdims=100, noise=0.0, outputlevel=1, observer=observer, cutoff=1e-8, eigsolve_krylovdim=15, eigsolve_verbosity=0, use_gpu=false)


end



#= depreciation calculation for moving out of old apartment
# depreciation parameters
factor_per_year = 0.9
number_of_years = 2 + 7/12  # moved into your old room June 2023, moving out Jan 2026

# initial values of items
kallax = 60
clothing_stand = 5 # if by clothing stand you mean the wooden boxes, this is nailed-together scrap wood
drawers = 2*15     # definitely not 30 each
desk = 10          # nowhere near new when I received it, definitely cheap to start
chair = 10         # also not new, had many holes, rips, and stains
total_value = kallax + clothing_stand + drawers + desk + chair
println("Total initial value is ",total_value)

# calculated depreciation
final_value = total_value * factor_per_year^(number_of_years)
depreciation_amount = total_value - final_value
println("Depreciation amount after $(number_of_years) years is ",depreciation_amount)

# moving and disposal costs
van_rental_cost = 60 / 4               # all these items makeup one quarter of a van load, van rental is €60
moving_items_cost = 15 * 0.25     # cost of physically moving all the items at €15 / hour for 15 minutes
total_moving_cost = van_rental_cost + moving_items_cost
println("Total moving cost is ",total_moving_cost)

# final payment due
println("Final Payment due is ",depreciation_amount - total_moving_cost)=#

























"fin"

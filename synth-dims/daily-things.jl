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
include_other_files(["synth-dims/long-range-ttn.jl","review-practice-codes/observables.jl","review-practice-codes/plottings.jl","synth-dims/hatsugai-mbcn.jl","other-funcs/basic-2d-plottings.jl"])
include_other_files(["synth-dims/oneD-effective-LR.jl","synth-dims/plottings-oneD.jl"])

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

# look at finite size scaling of commensurate filling interaction strength spectrum
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






























"fin"
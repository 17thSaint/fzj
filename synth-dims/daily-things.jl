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
include_other_files(["synth-dims/long-range-ttn.jl","review-practice-codes/observables.jl","review-practice-codes/plottings.jl","synth-dims/hatsugai-mbcn.jl"])
include_other_files(["synth-dims/oneD-effective-LR.jl","synth-dims/plottings-oneD.jl"])


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

    #=tws2 = range(0.1,0.3,length=11)
    tws = range(0.7,0.9,length=11)

    lambda1s::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
    lambda2s::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
    omegas::Matrix{ComplexF64} = zeros(ComplexF64,length(tws),length(tws2))
    ref_multis,rm1,rm2 = get_reference_multiplets(lx,ly,n)
    for (idx,tw1) in enumerate(tws)
        for (idx2,tw2) in enumerate(tws2)
            params_dict = Dict([("Lphys",lx),("Lsynth",ly),("particles",n),("tw2",tw2),("tw1",tw1),("if_remapping",false),("es_count",2),("nrgtol",1e-6),("mdim",200),("if_periodic_phys",true),("if_periodic_synth",true),("filling",0.5),("if_find_data",true),("if_save_data",true)])
            psis,rhos,nrgs,model_paras,if_found_data = run_normal_1deffmps(params_dict)
            
            #scatter3D(tw1,tw2,nrgs[1],c="b")
            #scatter3D(tw1,tw2,nrgs[2],c="g")
            #scatter3D(tw1,tw2,nrgs[3],c="r")

            if if_found_data
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

    tw1s = sort(unique(vcat(range(0.1,0.3,length=11),range(0.7,0.9,length=11),range(0.0,1.0,length=21))))
    tw2s = sort(unique(vcat(range(0.7,0.9,length=11),range(0.1,0.3,length=11),range(0.0,1.0,length=21))))

    lambda1s::Matrix{ComplexF64} = zeros(ComplexF64,length(tw1s),length(tw2s))
    lambda2s::Matrix{ComplexF64} = zeros(ComplexF64,length(tw1s),length(tw2s))
    omegas::Matrix{ComplexF64} = zeros(ComplexF64,length(tw1s),length(tw2s))
    for f in all_files
        d,m = read_data_jld2(dataloc * "/" * f; output_level=0)
        tw1,tw2 = m["twist_angle"]
        if tw1 in [0.33,0.67] || tw2 in [0.33,0.67]
            continue
        end
        idx = findfirst(x->tw1s[x]==tw1,1:length(tw1s))
        idx2 = findfirst(x->tw2s[x]==tw2,1:length(tw2s))
        lambda1s[idx,idx2] = m["gamma1"]
        lambda2s[idx,idx2] = m["gamma2"]
        omegas[idx,idx2] = m["omega"]
    end

    plot_omega(tw1s,tw2s,omegas; if_perfect_grid=false)
    plot_gamma(tw1s,tw2s,lambda1s,1)
    plot_gamma(tw1s,tw2s,lambda2s,2)
end






























"fin"
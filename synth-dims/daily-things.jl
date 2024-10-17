#####################################################
#=

This file contains daily things for synth-dims TTNs

Depends on:
    other-funcs/basic-2d-stuff.jl
    review-practice-codes/ttn.jl
    review-practice-codes/observables.jl
    review-practice-codes/plottings.jl

=#
######################################################

include("../other-funcs/include-other-files.jl")
include_other_files(["review-practice-codes/ttn.jl","synth-dims/long-range-ttn.jl","review-practice-codes/observables.jl","review-practice-codes/plottings.jl"])



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
        if !(m["onsite_strength"] in chosen_intstrens)
            continue
        end
        append!(intstrens,[m["onsite_strength"]])
        for i in 1:3
            keyname = i == 1 ? "observer" : "observer_$i"
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






























"fin"
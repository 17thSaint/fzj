#####################################################
#=

This file contains daily things for fPEPS, ### header + if true/false blocks run from the REPL
(pwd inside fzj/, for include_other_files)

Depends on:
    other-funcs/include-other-files.jl (include_other_files)
    fpeps/load-fpeps.jl (all of fpeps)
    Random

=#
######################################################

include("../other-funcs/include-other-files.jl")
include_other_files(["fpeps/load-fpeps.jl"])

using Random,PyPlot,Statistics


### test CDW order parameters (S(pi), correlation ratio, Binder) from samples against exact enumeration,
#= random 4x4 N=2 PEPS (120 configurations in the sector)
if false
    Random.seed!(1)
    lx,ly,n = 4,4,2
    bond_dim = 3
    nsamples = 4000

    hilbert = ITensors.siteinds("S=1/2",lx,ly)
    peps = PEPS(ComplexF64, hilbert; bond_dim=bond_dim)
    max_counts = particle_sector_counts(lx, ly, n)

    occs_ex, weights_ex = exact_fpeps_occupations(peps, n)
    column_corr_ex, evencolumn_dist_ex = column_cdw_statistics(occs_ex, weights_ex)

    @time occs, weights = sample_fpeps_occupations(peps, nsamples; max_counts)
    column_corr, evencolumn_dist = column_cdw_statistics(occs, weights)
    println("$(size(occs,3)) samples on $(Threads.nthreads()) threads, particle numbers $(unique(sum(occs,dims=(1,2))))")
    println("max |column_corr - exact| = $(maximum(abs.(column_corr .- column_corr_ex))), max |evencolumn_dist - exact| = $(maximum(abs.(evencolumn_dist .- evencolumn_dist_ex)))")

    s_pi_ex, corr_ratio_ex, binder_ex = cdw_order_parameters(column_corr_ex, evencolumn_dist_ex)
    ops = fpeps_cdw_order_parameters(occs, weights)
    println("S(pi):       exact $(round(s_pi_ex, digits=4)), sampled $(round(ops.s_pi, digits=4)) ± $(round(ops.s_pi_err, digits=4))")
    println("corr ratio:  exact $(round(corr_ratio_ex, digits=4)), sampled $(round(ops.corr_ratio, digits=4)) ± $(round(ops.corr_ratio_err, digits=4))")
    println("Binder:      exact $(round(binder_ex, digits=4)), sampled $(round(ops.binder, digits=4)) ± $(round(ops.binder_err, digits=4))")
end=#

### CDW order parameters of saved fPEPS runs, statistics saved into the metadata of their fpeps-*.jld2 (writes to cluster-data),
#= read back with fpeps_cdw_order_parameters(metadata) or read_column_cdw_statistics(metadata)
if false
    dataloc = get_folder_location("cluster-data/synth-dims/fpeps")
    nsamples = 4000
    for f in filter(f -> startswith(f,"fpeps-"), readdir(dataloc))
        ops = fpeps_cdw_order_parameters(joinpath(dataloc,f), nsamples; if_save_data=true, output_level=2)
    end
end=#


# check convergence of 8x8 cluster data
if true
    lx,ly,n = 8,8,4
    dataloc = get_folder_location("cluster-data/synth-dims/fpeps/dd-ints")
    pdict = Dict([("Lx",lx),("Ly",ly),("particles",n)])
    all_files = find_data_file(pdict, "fpeps", dataloc; file_type="jld2")
    
    for f in all_files

        d,m = read_data(joinpath(dataloc, f); output_level=0)
        
        nrg_history = d["history"]
        average_count = 20
        running_average = [mean(nrg_history.energy[max(1,i-average_count):i]) for i in 1:length(nrg_history.energy)]
        fig = figure()
        plot(nrg_history.niter, running_average .- minimum(nrg_history.energy), "-o")
        yscale("log")
        ylabel("Energy - min(Energy)")
        xlabel("Iteration Number")
        title("fPEPS convergence for MagSpac = $(m["magnetic_spacing"])")
    end
end









"fin"

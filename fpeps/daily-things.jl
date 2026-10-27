#####################################################
#=

This file contains daily things for fPEPS

Depends on:
    fpeps/build-ham.jl
    fpeps/observables.jl
    other-funcs/data-storage-funcs.jl

=#
######################################################

include("../other-funcs/include-other-files.jl")
include_other_files(["fpeps/build-ham.jl","fpeps/observables.jl","other-funcs/data-storage-funcs.jl"])

using QuantumNaturalfPEPS, ITensors, ITensorMPS, JLD2, Random


### test CDW order parameters (S(pi), correlation ratio, Binder) from samples against exact enumeration,
### random 4x4 N=2 PEPS (120 configurations in the sector)
if true
    Random.seed!(1)
    lx,ly,n = 4,4,2
    bond_dim = 3
    nsamples = 4000

    hilbert = ITensors.siteinds("S=1/2",lx,ly)
    peps = PEPS(ComplexF64, hilbert; bond_dim=bond_dim)
    max_counts = [n, lx*ly - n]

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
end




"fin"

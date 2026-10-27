#####################################################
#=

This file contains observables for fPEPS wavefunctions, estimated from samples of |ψ|²

Depends on:
	other-funcs/basic-2d-stuff.jl (LinearAlgebra)
	other-funcs/basic-2d-observables.jl (column_cdw_statistics, cdw_order_parameters)
	other-funcs/data-storage-funcs.jl (read_data, add_wavefunc_to_filepath, remove_wavefunc_from_filename)
	QuantumNaturalfPEPS

Samples S come from get_sample (S[x,y] = 0 is Up = occupied, so n = 1 - S) and are reweighted by
|ψ(S)|²/p_c(S), since the row-by-row sampling with a finite contract_dim is only approximately |ψ|².
The CDW order parameters are shared with ED (other-funcs/basic-2d-observables.jl), only the column
occupations come from samples instead of the ED basis

=#
######################################################


include("../other-funcs/basic-2d-stuff.jl")
include("../other-funcs/basic-2d-observables.jl")

using QuantumNaturalfPEPS


# occs[x,y,i] = n_{x,y} of sample i and the importance weights (mean 1) for nsamples samples of the peps.
# max_counts = [N, Lx*Ly - N] samples in the N particle sector, as in execute-fpeps.jl
function sample_fpeps_occupations(peps, nsamples::Int; max_counts=nothing, threaded::Bool=true)
    Lx,Ly = size(peps)
    QuantumNaturalfPEPS.update_double_layer_envs!(peps)

    nr_threads = threaded ? Threads.nthreads() : 1
    k = ceil(Int, nsamples / nr_threads)
    occs = zeros(Int8, Lx, Ly, k*nr_threads)
    logψs = zeros(ComplexF64, k*nr_threads)
    logpcs = zeros(Float64, k*nr_threads)

    Threads.@threads for t in 1:nr_threads
        for i in (1+(t-1)*k):(t*k)
            S, logpc, env_top = QuantumNaturalfPEPS.get_sample(peps; max_counts)
            logψs[i], _ = QuantumNaturalfPEPS.get_logψ_and_envs(peps, S, env_top)
            logpcs[i] = real(logpc)
            occs[:,:,i] .= 1 .- S
        end
    end

    weights = QuantumNaturalfPEPS.compute_importance_weights(logψs, logpcs)
    n_eff = sum(weights)^2 / sum(abs2, weights)
    # the default combinatorial sector proposal already gives n_eff ~ 0.4 nsamples for a random 4x4 N=2 PEPS
    n_eff < 0.1 * length(weights) && @warn "effective sample size $(round(n_eff, digits=1)) of $(length(weights)), sampling is far from |ψ|²"

    return occs, weights
end

# all configurations of the N particle sector with weights |ψ(S)|² from exact contraction, in the format of
# sample_fpeps_occupations. Only for small systems (benchmarking the sampling)
function exact_fpeps_occupations(peps, N::Int)
    Lx,Ly = size(peps)
    configs = [c for c in 0:(2^(Lx*Ly)-1) if count_ones(c) == N]

    occs = zeros(Int8, Lx, Ly, length(configs))
    logψs = zeros(ComplexF64, length(configs))
    for (i,c) in enumerate(configs)
        occs[:,:,i] .= reshape(digits(Int8, c, base=2, pad=Lx*Ly), Lx, Ly)
        logψs[i] = QuantumNaturalfPEPS.logψ_exact(peps, Int64.(1 .- occs[:,:,i]))
    end
    log_probs = 2 .* real.(logψs)

    return occs, exp.(log_probs .- maximum(log_probs))
end

# column occupations ncols[x,i] = sum_y n_{x,y} of sample i, for column_cdw_statistics in other-funcs/basic-2d-observables.jl
column_occupations(occs::Array{Int8,3}) = dropdims(sum(occs,dims=2),dims=2)

# column statistics of weighted samples, see column_cdw_statistics in other-funcs/basic-2d-observables.jl
column_cdw_statistics(occs::Array{Int8,3},weights::Vector{Float64}) = column_cdw_statistics(column_occupations(occs),weights)

# S(pi), correlation ratio and Binder cumulant (cdw_order_parameters) with jackknife errors over nbins bins of samples.
# The samples are independent (no Markov chain), so the bins only set the cost of the jackknife. The column statistics
# are computed once per bin, the leave-one-bin-out ones are the weighted totals minus that bin
function fpeps_cdw_order_parameters(occs::Array{Int8,3}, weights::Vector{Float64}; nbins::Int=20)
    ncols = column_occupations(occs)
    nbins = min(nbins, size(ncols,2))
    bins = [b:nbins:size(ncols,2) for b in 1:nbins]

    bin_weights = [sum(view(weights,bin)) for bin in bins]
    bin_stats = [column_cdw_statistics(ncols[:,bin], weights[bin]) for bin in bins]
    total_weight = sum(bin_weights)
    total_corr = sum(bin_weights[b] .* bin_stats[b][1] for b in 1:nbins)
    total_dist = sum(bin_weights[b] .* bin_stats[b][2] for b in 1:nbins)

    values = collect(cdw_order_parameters(total_corr ./ total_weight, total_dist ./ total_weight))
    jack = [collect(cdw_order_parameters((total_corr .- bin_weights[b] .* bin_stats[b][1]) ./ (total_weight - bin_weights[b]),
                                         (total_dist .- bin_weights[b] .* bin_stats[b][2]) ./ (total_weight - bin_weights[b])))
            for b in 1:nbins]
    jack_mean = sum(jack) / nbins
    errors = sqrt.((nbins - 1) / nbins .* sum(j -> (j .- jack_mean).^2, jack))

    s_pi, corr_ratio, binder = values
    return (; s_pi, corr_ratio, binder, s_pi_err=errors[1], corr_ratio_err=errors[2], binder_err=errors[3])
end

# CDW order parameters of a saved fPEPS (fpeps-*.jld2 or wavefuncfpeps-*.jld2), sampled in the particle sector of its metadata
function fpeps_cdw_order_parameters(filepath::String, nsamples::Int; nbins::Int=20, threaded::Bool=true)
    filepath = joinpath(dirname(filepath), remove_wavefunc_from_filename(basename(filepath)))
    _, metadata = read_data(filepath)
    wavefunc_data, _ = read_data(add_wavefunc_to_filepath(filepath))
    lx,ly,n = metadata["Lx"],metadata["Ly"],metadata["particles"]

    occs, weights = sample_fpeps_occupations(wavefunc_data["peps"], nsamples; max_counts=[n, lx*ly - n], threaded)
    return fpeps_cdw_order_parameters(occs, weights; nbins)
end

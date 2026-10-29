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
# max_counts = [N, Lx*Ly - N] samples in the N particle sector, as in execute-fpeps.jl.
# output_level 1 prints start/end and the effective sample size, 2 also the progress every 10% of the samples
function sample_fpeps_occupations(peps, nsamples::Int; max_counts=nothing, threaded::Bool=true, output_level::Int=0)
    Lx,Ly = size(peps)
    start_time = time()
    QuantumNaturalfPEPS.update_double_layer_envs!(peps)
    output_level > 0 && println("Double-layer environments of the $(Lx)x$(Ly) PEPS in $(round(time()-start_time,digits=1)) s")

    nr_threads = threaded ? Threads.nthreads() : 1
    k = ceil(Int, nsamples / nr_threads)
    occs = zeros(Int8, Lx, Ly, k*nr_threads)
    logψs = zeros(ComplexF64, k*nr_threads)
    logpcs = zeros(Float64, k*nr_threads)

    # samples done over all threads, for the progress output
    nr_done = Threads.Atomic{Int}(0)
    report_every = max(1, length(logpcs) ÷ 10)
    output_level > 0 && println("Sampling $(length(logpcs)) configurations on $nr_threads threads")
    flush(stdout)
    start_time = time()

    Threads.@threads for t in 1:nr_threads
        for i in (1+(t-1)*k):(t*k)
            S, logpc, env_top = QuantumNaturalfPEPS.get_sample(peps; max_counts)
            logψs[i], _ = QuantumNaturalfPEPS.get_logψ_and_envs(peps, S, env_top)
            logpcs[i] = real(logpc)
            occs[:,:,i] .= 1 .- S

            if output_level > 1
                done = Threads.atomic_add!(nr_done, 1) + 1
                if done % report_every == 0
                    elapsed = time() - start_time
                    println("  $done/$(length(logpcs)) samples ($(round(Int,100*done/length(logpcs))) %), $(round(elapsed,digits=1)) s, ~$(round(elapsed*(length(logpcs)-done)/done,digits=1)) s left")
                    flush(stdout)
                end
            end
        end
    end

    weights = QuantumNaturalfPEPS.compute_importance_weights(logψs, logpcs)
    n_eff = sum(weights)^2 / sum(abs2, weights)
    # the default combinatorial sector proposal already gives n_eff ~ 0.4 nsamples for a random 4x4 N=2 PEPS
    n_eff < 0.1 * length(weights) && @warn "effective sample size $(round(n_eff, digits=1)) of $(length(weights)), sampling is far from |ψ|²"
    output_level > 0 && println("Sampled $(length(weights)) in $(round(time()-start_time,digits=1)) s, effective sample size $(round(n_eff,digits=1)) ($(round(Int,100*n_eff/length(weights))) %)")

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

# column statistics per bin of samples (sample i in bin mod1(i,nbins)) and the bin weights, for the jackknife errors
function binned_column_cdw_statistics(occs::Array{Int8,3}, weights::Vector{Float64}; nbins::Int=20)
    ncols = column_occupations(occs)
    nbins = min(nbins, size(ncols,2))
    bins = [b:nbins:size(ncols,2) for b in 1:nbins]
    bin_stats = [column_cdw_statistics(ncols[:,bin], weights[bin]) for bin in bins]
    return first.(bin_stats), last.(bin_stats), [sum(view(weights,bin)) for bin in bins]
end

# weighted sums over bins of the column statistics, the totals are these over sum(bin_weights)
bin_sums(corr_bins, dist_bins, bin_weights) = sum(bin_weights .* corr_bins), sum(bin_weights .* dist_bins)

# S(pi), correlation ratio and Binder cumulant (cdw_order_parameters) of all bins with jackknife errors. The samples are
# independent (no Markov chain), so the bins only set the cost of the jackknife; the leave-one-bin-out statistics are
# the weighted sums minus that bin
function jackknife_cdw_order_parameters(corr_bins::Vector{Matrix{Float64}}, dist_bins::Vector{Vector{Float64}}, bin_weights::Vector{Float64})
    nbins = length(bin_weights)
    total_weight = sum(bin_weights)
    total_corr, total_dist = bin_sums(corr_bins, dist_bins, bin_weights)

    values = collect(cdw_order_parameters(total_corr ./ total_weight, total_dist ./ total_weight))
    jack = [collect(cdw_order_parameters((total_corr .- bin_weights[b] .* corr_bins[b]) ./ (total_weight - bin_weights[b]),
                                         (total_dist .- bin_weights[b] .* dist_bins[b]) ./ (total_weight - bin_weights[b])))
            for b in 1:nbins]
    jack_mean = sum(jack) / nbins
    errors = sqrt.((nbins - 1) / nbins .* sum(j -> (j .- jack_mean).^2, jack))

    s_pi, corr_ratio, binder = values
    return (; s_pi, corr_ratio, binder, s_pi_err=errors[1], corr_ratio_err=errors[2], binder_err=errors[3])
end

fpeps_cdw_order_parameters(occs::Array{Int8,3}, weights::Vector{Float64}; nbins::Int=20) =
    jackknife_cdw_order_parameters(binned_column_cdw_statistics(occs, weights; nbins)...)

# metadata entries of the GS (k = 1): the ED keys "column_corr1", "evencolumn_dist1" over all samples (column_cdw_statistics_dict)
# plus the per-bin statistics for the jackknife, the number of samples and the optimisation step niter they were drawn at
function fpeps_cdw_statistics_dict(corr_bins::Vector{Matrix{Float64}}, dist_bins::Vector{Vector{Float64}}, bin_weights::Vector{Float64}, nsamples::Int, niter)
    total_corr, total_dist = bin_sums(corr_bins, dist_bins, bin_weights)
    stats = column_cdw_statistics_dict(total_corr ./ sum(bin_weights), total_dist ./ sum(bin_weights), 1)
    return merge(stats, Dict{String,Any}([("column_corr_bins1",corr_bins),("evencolumn_dist_bins1",dist_bins),("bin_weights1",bin_weights),
        ("nsamples1",nsamples),("niter_sampled1",niter)]))
end

# CDW order parameters with jackknife errors back from metadata saved by fpeps_cdw_order_parameters(filepath, ...; if_save_data=true)
function fpeps_cdw_order_parameters(metadata::Dict)
    haskey(metadata,"column_corr_bins1") || error("no binned column statistics in the metadata")
    return jackknife_cdw_order_parameters(Matrix{Float64}.(metadata["column_corr_bins1"]), Vector{Float64}.(metadata["evencolumn_dist_bins1"]),
        Vector{Float64}(metadata["bin_weights1"]))
end

# CDW order parameters of a saved fPEPS (fpeps-*.jld2 or wavefuncfpeps-*.jld2), sampled in the particle sector of its metadata.
# if_save_data adds the statistics (fpeps_cdw_statistics_dict) to the metadata of the fpeps-*.jld2 file, overwriting older ones.
# output_level 1 prints the run, sampling start/end and the result, 2 also the sampling progress (see sample_fpeps_occupations)
function fpeps_cdw_order_parameters(filepath::String, nsamples::Int; nbins::Int=20, threaded::Bool=true, if_save_data::Bool=false, output_level::Int=0)
    filepath = joinpath(dirname(filepath), remove_wavefunc_from_filename(basename(filepath)))
    _, metadata = read_data(filepath; output_level)
    wavefunc_data, _ = read_data(add_wavefunc_to_filepath(filepath); output_level)
    lx,ly,n = metadata["Lx"],metadata["Ly"],metadata["particles"]
    output_level > 0 && println("CDW order parameters of $(basename(filepath)): $(lx)x$(ly) N=$n, D = $(metadata["bond_dim"]), niter = $(metadata["niter"])")

    occs, weights = sample_fpeps_occupations(wavefunc_data["peps"], nsamples; max_counts=[n, lx*ly - n], threaded, output_level)
    binned_stats = binned_column_cdw_statistics(occs, weights; nbins)
    if if_save_data
        modify_data(fpeps_cdw_statistics_dict(binned_stats..., size(occs,3), metadata["niter"]), filepath, "metadata")
        output_level > 0 && println("Saved the column statistics to the metadata of $(basename(filepath))")
    end

    ops = jackknife_cdw_order_parameters(binned_stats...)
    if output_level > 0
        println("S(pi) = $(round(ops.s_pi,digits=4)) ± $(round(ops.s_pi_err,digits=4)), corr ratio = $(round(ops.corr_ratio,digits=4)) ± $(round(ops.corr_ratio_err,digits=4)), Binder = $(round(ops.binder,digits=4)) ± $(round(ops.binder_err,digits=4))")
        flush(stdout)
    end
    return ops
end

#####################################################
#=

This file contains plots for the dipolar interaction length scaling paper.

Depends on:
    synth-dims/long-range-ttn.jl
    review-practice-codes/observables.jl
    review-practice-codes/plottings.jl

    exact-diag/execute-ed.jl
    exact-diag/observables.jl
    exact-diag/plottings.jl

=#
######################################################

include("../other-funcs/include-other-files.jl")
include_other_files(["synth-dims/long-range-ttn.jl","review-practice-codes/observables.jl","other-funcs/basic-2d-observables.jl","exact-diag/execute-ed.jl","exact-diag/observables.jl"])
include_other_files(["other-funcs/basic-2d-plottings.jl","review-practice-codes/plottings.jl","exact-diag/plottings.jl"])

# make plot of shifted zeeman minima along with interaction profile
function plot_shifted_zeeman_minima_and_interaction_profile()
    num_minima = 11
    omega = 20.0
    shift_val = 1.0
    cols = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf","#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"]
    xs = range(-1.5, 12.5, length=200)
    central_minimum = [0.5 * omega * x^2 for x in xs]

    fig,axs = subplots(2, 1, figsize=(6, 8))
    axs[1].plot(xs,central_minimum,c=cols[1],label="m=0")
    for i in 1:num_minima-1
        shifted_minimum = [0.5 * omega * (x - i*shift_val)^2 for x in xs]
        axs[1].plot(xs,shifted_minimum,c=cols[i+1],label="m=$(i)")
    end
    axs[1].set_xlabel("Position (arb. units)")
    axs[1].set_ylabel("Potential Energy (arb. units)")
    axs[1].legend()
    axs[1].set_ylim(-0.5, 3)
    axs[1].set_xlim(-0.8, 11)


    # make plot of interaction profile
    grads_vals = [0.5, 1.0, 2.0]
    cols_int = ["k","r","g"]
    grads_x_starts = [0.5, 0.2, 0.1]
    for i in 1:length(grads_vals)
        grads = grads_vals[i]
        xs_int = vcat([grads_x_starts[i]], [x for x in 1:11])
        interaction_profile = [300 / (grads * r)^3 for r in xs_int]
        axs[2].plot(xs_int,interaction_profile,"-p",c=cols_int[i],label="s=$(grads)")
    end
    axs[2].set_xlabel("Zeeman Sublevel Separation, "*L"\vert m - m' \vert")
    axs[2].set_ylabel("Interaction Strength (arb. units)")
    axs[2].set_xlim(-0.8, 11)
    axs[2].legend()
    axs[2].set_yscale("log")
    tight_layout()
end





























"fin"
#####################################################
#=

This file contains the figures of merit used to optimize, with QuOCS dCRAB, the magnetic_spacing
ramp of the dd interaction profile used in tevo-daily-things.jl:

    dd gs manifold (a = a_start) --(magnetic_spacing ramp a(t))--> dd gs manifold (a = a_end)

with either a(t) itself as the control (compute_fidelity_magspac_ramp, config_magspacRamp.py) or
the magnetic gradient B(t) that drives it (compute_fom_magspac_gradient, config_magspacGrad.py).

Unlike the other control-functions files the evolution does not go through run_timeevo.
Two reasons, both specific to ramping the spacing:
  - speed: timeham rebuilds H with buildHam at every half step because a(t) changes every
    sample; here H is linear in the coupling vector, H(a) = H0 + sum_k U_k(a) M_k, so the
    matrices are built once in setup;
  - stability: the pulse is allowed down to a = 0.1, where U1 = intstren/a^3 is 125x its
    a = 0.5 value. RK4 at a fixed dt is only stable for the t = 0 couplings, so the evolution
    is a unitary split-operator scheme instead, unconditionally stable.
The interaction matrices M_k are diagonal in the Fock basis, so H(t) = H0 + D(t) with H0
(hopping) constant and D(t) = sum_k U_k(a(t)) diag(M_k). Each substep of length d is the Strang
splitting exp(-i D d/2) exp(-i H0 d) exp(-i D d/2), with a(t) linear between the half-step
samples; exp(-i H0 d) is computed once per substep size, so a substep is one 120x120 matmul
plus diagonal phases. The splitting error grows with the spread of D, so each interval is cut
into m = ceil(h * spread(D) / theta) substeps: one everywhere a >= 0.5, more where the pulse
dips. Checked (theta = 2) against a 64-substep reference and against run_timeevo on the ZVD
pulse: agreement to ~2e-6 in the fidelity, ~10 ms/eval (vs ~2 s with exact exponentials).

Depends on:
    execute-ed.jl (run_normal_ed, buildHam, long_range_scaling)
    control-functions.jl (optimization_ed_params, groundstate_manifold_fidelity)
    time-evolution.jl (magnetic_gradient_response)

=#
######################################################

using LinearAlgebra

# the dd coupling vector at spacing a, exactly as timeham builds it (including
# long_range_scaling's 5-digit rounding), so the fast evolution sees the same U as run_timeevo.
# applyHam drops every coupling at or below interaction_cutoff, so those are zeroed here too
# (at large a the rounded tail vanishes)
function magspac_couplings(a::Real,n_u::Int,Ly::Int,intstren::Float64,cutoff::Float64)
    u = long_range_scaling(n_u-1,Ly,intstren; scaling="dd",magnetic_spacing=Float64(a))
    return ifelse.(abs.(u) .> cutoff, u, 0.0)
end

# pulse-independent part of the figure of merit, run once per optimization: the H = H0 +
# sum_k U_k M_k decomposition and the endpoint groundstate manifolds
function setup_magspac_ramp(parameters_dictionary)

    speccount::Int = Int(parameters_dictionary["speccount"])
    intstren::Float64 = Float64(parameters_dictionary["intstren"])
    a_start::Float64 = Float64(parameters_dictionary["a_start"])
    a_end::Float64 = Float64(parameters_dictionary["a_end"])

    pdict = optimization_ed_params(parameters_dictionary,"scaling_type"=>"dd","magnetic_spacing"=>a_start,"interaction_strength"=>intstren)
    _,_,_,_,_,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=0)
    Ly::Int = lattice_params["Ly"]

    # H is linear in the coupling vector, so a handful of buildHam calls give the whole family
    # exactly; checked against the starting H buildHam made. Unit vectors do not work as the
    # probes: applyHam takes the interaction range from the number of couplings above the
    # cutoff, so e_k would read as range 0. Perturbing an all-ones vector keeps the full range
    n_u = length(hamilt_params["U"])
    cutoff::Float64 = hamilt_params["interaction_cutoff"]
    hp = copy(hamilt_params)
    hp["U"] = ones(n_u)
    H_ones = Matrix(buildHam(lattice_params,hp; output_level=0))
    Ms = map(1:n_u) do k
        hp["U"] = [i == k ? 2.0 : 1.0 for i in 1:n_u]
        Matrix(buildHam(lattice_params,hp; output_level=0)) - H_ones
    end
    H0 = H_ones - sum(Ms)
    ham(a) = H0 + sum(magspac_couplings(a,n_u,Ly,intstren,cutoff) .* Ms)
    @assert ham(a_start) ≈ Matrix(hamilt_params["H"]) "H is not linear in the coupling vector U; the H0 + sum_k U_k M_k decomposition is invalid"
    @assert all(isdiag,Ms) "interaction matrices are not diagonal; the split-operator evolution assumes they are"

    # both endpoints have a speccount-fold degenerate groundstate, which a single Lanczos draw
    # sometimes returns only part of; the space is small, so diagonalize densely
    starting_states, target_states = map((a_start,a_end)) do a
        eig = eigen(Hermitian(ham(a)))
        nrgs = eig.values
        @assert nrgs[speccount] - nrgs[1] < 1e-8 && nrgs[speccount+1] - nrgs[speccount] > 1e-3 "groundstate manifold at a = $(a) is not an isolated $(speccount)-fold degenerate level: $(nrgs[1:speccount+1])"
        [Vector{ComplexF64}(eig.vectors[:,i]) for i in 1:speccount]
    end

    # D(a) = diag_matrix * U(a): column k is the diagonal of M_k
    diag_matrix = reduce(hcat,[real(diag(M)) for M in Ms])

    # indexable as setup[1], setup[2] for the starting / target manifolds; expH0_cache holds
    # exp(-i H0 d) per substep size d, filled on demand by evolve_magspac_ramp
    return (starting_states=starting_states,target_states=target_states,H0=Hermitian(H0),diag_matrix=diag_matrix,
            Ly=Ly,intstren=intstren,cutoff=cutoff,expH0_cache=Dict{Float64,Matrix{ComplexF64}}())
end

# evolve the starting manifold through the spacing samples a (spacing h apart in time) and
# return the final states, by the adaptive Strang splitting described in the header
function evolve_magspac_ramp(spacing_values::AbstractVector{<:Real},h::Float64,setup; theta::Float64=2.0)
    n_u = size(setup.diag_matrix,2)
    diag_energies(a) = setup.diag_matrix * magspac_couplings(a,n_u,setup.Ly,setup.intstren,setup.cutoff)

    psi = reduce(hcat,setup.starting_states)
    d_prev = diag_energies(spacing_values[1])
    for i in 1:length(spacing_values)-1
        d_next = diag_energies(spacing_values[i+1])
        spread = max(maximum(d_prev) - minimum(d_prev), maximum(d_next) - minimum(d_next))
        m = max(1, ceil(Int, h*spread/theta))
        dsub = h/m
        expH0 = get!(() -> exp(-im*dsub*setup.H0), setup.expH0_cache, dsub)
        for j in 1:m
            w = (j - 0.5)/m
            half_phase = cis.(-dsub/2 .* ((1 - w) .* d_prev .+ w .* d_next))
            psi = half_phase .* (expH0 * (half_phase .* psi))
        end
        d_prev = d_next
    end

    return [psi[:,i] for i in 1:size(psi,2)]
end

# a(t) as the control: fidelity of the spacing samples QuOCS hands over (possibly complex-typed)
# on the dt/2 grid
function compute_fidelity_magspac_ramp(pulses,parameters_dictionary,setup)
    spacing_values = Float64.(real.(collect(pulses[1])))
    final_states = evolve_magspac_ramp(spacing_values,Float64(parameters_dictionary["dt"])/2,setup)
    return groundstate_manifold_fidelity(final_states,setup.target_states)
end

# B(t) as the control: fidelity of the spacing ramp the gradient drives from rest at a_start,
#     a(t) = a_start + (1/w) int_0^t B(t') sin(w(t - t')) dt'   (magnetic_gradient_response),
# minus the residual ringing. The optimization pins B(T) to the holding gradient
# w^2 (a_end - a_start), so after T the mode oscillates about a_end with amplitude
#     R = sqrt((a(T) - a_end)^2 + (adot(T)/w)^2),
# which is zero only if the lattice arrives at a_end at rest and stays there. Unweighted, both
# terms are dimensionless / in units of spacing and of the same size as the fidelity gains.
# Pulses driving a below a_floor are not evolved (U ~ 1/a^3 diverges at 0, and the split-operator
# substep count with it); they score below any evolved pulse, graded by the depth of the dip.
# Returns (FoM, fidelity, ringing, min a); fidelity is NaN for a rejected pulse
function magspac_gradient_fom(gradient_values::AbstractVector{<:Real},parameters_dictionary,setup)
    h = Float64(parameters_dictionary["dt"]) / 2
    w = Float64(parameters_dictionary["trap_frequency"])
    a_start = Float64(parameters_dictionary["a_start"])
    a_end = Float64(parameters_dictionary["a_end"])
    a_floor = Float64(parameters_dictionary["a_floor"])

    x, xdot = magnetic_gradient_response(gradient_values,h,w)
    a = a_start .+ x
    ringing = hypot(a[end] - a_end, xdot[end] / w)
    min_a = minimum(a)
    min_a < a_floor && return (-1.0 - (a_floor - min_a), NaN, ringing, min_a)

    fidelity = groundstate_manifold_fidelity(evolve_magspac_ramp(a,h,setup),setup.target_states)
    return (fidelity - ringing, fidelity, ringing, min_a)
end

compute_fom_magspac_gradient(pulses,parameters_dictionary,setup) = magspac_gradient_fom(Float64.(real.(collect(pulses[1]))),parameters_dictionary,setup)

# The ZVD shaper for the trapped mode: impulses of relative size [1/4, 1/2, 1/4] at 0, pi/w, 2pi/w,
# which cancel the residual oscillation and its first derivative in w. Convolved with a step it is
# the gradient staircase B = Bf [1/4, 3/4, 1], Bf = w^2 (a_end - a_start), which settles the spacing
# at a_end from 2pi/w on.
zvd_impulses(trap_frequency::Float64) = ([0.25, 0.5, 0.25], [0.0, 1.0, 2.0] .* (pi/trap_frequency))

# the ZVD spacing ramp, a(t) = a_start + da sum_{t_i <= t} A_i (1 - cos w(t - t_i)): the closed-form
# response to the staircase, sampled on the half-step grid (initial guess when a(t) is the control)
function zvd_spacing_pulse(a_start::Float64,a_end::Float64,trap_frequency::Float64,ramptime::Float64,dt::Float64)
    amps, timps = zvd_impulses(trap_frequency)
    @assert timps[end] <= ramptime*(1 + 1e-12) "ZVD pulse settles at $(timps[end]), after the end of the ramp $(ramptime)"
    return [a_start + (a_end - a_start)*sum((A*(1 - cos(trap_frequency*(t - ti))) for (A,ti) in zip(amps,timps) if ti <= t); init=0.0) for t in halfstep_times(ramptime,dt)]
end

# the ZVD gradient staircase on the half-step grid (initial guess when B(t) is the control).
# An interior switch landing exactly on a sample gets the average of the levels either side, which
# makes the trapezoid integrate the step exactly (see the ZVD three-step block of
# tevo-daily-things.jl); a switch between samples is smeared across one half-step instead. A switch
# on the last sample takes the new level, so a ramp ending at 2pi/w ends on the holding gradient Bf
function zvd_gradient_pulse(a_start::Float64,a_end::Float64,trap_frequency::Float64,ramptime::Float64,dt::Float64)
    amps, timps = zvd_impulses(trap_frequency)
    @assert timps[end] <= ramptime*(1 + 1e-12) "ZVD pulse settles at $(timps[end]), after the end of the ramp $(ramptime)"
    bf = trap_frequency^2 * (a_end - a_start)
    h = dt/2
    times = halfstep_times(ramptime,dt)
    level(t) = bf * sum((A for (A,ti) in zip(amps,timps) if ti <= t + 1e-9h); init=0.0)
    return map(enumerate(times)) do (i,t)
        on_interior_switch = 1 < i < length(times) && any(ti -> abs(t - ti) < 1e-9h, timps)
        on_interior_switch ? (level(t - h/2) + level(t))/2 : level(t)
    end
end

"fin"

#####################################################
#=

This file contains the figure-of-merit function used to optimize, with QuOCS dCRAB,
the magnetic_spacing ramp of the dd interaction profile used in tevo-daily-things.jl:

    dd gs manifold (a = a_start) --(magnetic_spacing ramp a(t))--> dd gs manifold (a = a_end)

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
pulse: agreement to ~2e-6 in the fidelity, ~20 ms/eval (vs ~2 s with exact exponentials).

Depends on:
    execute-ed.jl (run_normal_ed, buildHam, long_range_scaling)
    control-functions.jl (groundstate_manifold_fidelity)

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

    pdict = Dict{String,Any}(
        "output_level"=>0,"Lx"=>Int(parameters_dictionary["Lx"]),"Ly"=>Int(parameters_dictionary["Ly"]),"N"=>Int(parameters_dictionary["N"]),
        "lr"=>"all","if_periodic_x"=>true,"if_periodic_y"=>true,"hopping_anisotropy"=>1.0,"filling"=>0.5,
        "scaling_type"=>"dd","magnetic_spacing"=>a_start,"interaction_strength"=>intstren,
        "nev"=>speccount,"if_find_data"=>false,"if_save_data"=>false,
    )
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

    # both endpoints have a speccount-fold degenerate groundstate, which a single Lanczos draw
    # sometimes returns only part of; the space is small, so diagonalize densely
    endpoint_states = map((a_start,a_end)) do a
        eig = eigen(Hermitian(ham(a)))
        nrgs = eig.values
        @assert nrgs[speccount] - nrgs[1] < 1e-8 && nrgs[speccount+1] - nrgs[speccount] > 1e-3 "groundstate manifold at a = $(a) is not an isolated $(speccount)-fold degenerate level: $(nrgs[1:speccount+1])"
        [Vector{ComplexF64}(eig.vectors[:,i]) for i in 1:speccount]
    end

    @assert all(isdiag,Ms) "interaction matrices are not diagonal; the split-operator evolution assumes they are"
    diag_ms = [real(diag(M)) for M in Ms]
    # exp(-i H0 d) per substep size d, filled on demand by evolve_magspac_ramp
    expH0_cache = Dict{Float64,Matrix{ComplexF64}}()

    return (endpoint_states[1],endpoint_states[2],H0,diag_ms,Ly,intstren,cutoff,expH0_cache)
end

# evolve the starting manifold through the spacing samples a (spacing h apart in time) and
# return the final states, by the adaptive Strang splitting described in the header
function evolve_magspac_ramp(spacing_values::AbstractVector{<:Real},h::Float64,setup; theta::Float64=2.0)
    starting_states,_,H0,diag_ms,Ly,intstren,cutoff,expH0_cache = setup
    n_u = length(diag_ms)
    diag_energies(a) = sum(magspac_couplings(a,n_u,Ly,intstren,cutoff) .* diag_ms)

    psi = reduce(hcat,starting_states)
    d_prev = diag_energies(spacing_values[1])
    for i in 1:length(spacing_values)-1
        d_next = diag_energies(spacing_values[i+1])
        spread = max(maximum(d_prev) - minimum(d_prev), maximum(d_next) - minimum(d_next))
        m = max(1, ceil(Int, h*spread/theta))
        dsub = h/m
        expH0 = get!(() -> exp(-im*dsub*Hermitian(H0)), expH0_cache, dsub)
        for j in 1:m
            w = (j - 0.5)/m
            half_phase = cis.(-dsub/2 .* ((1 - w) .* d_prev .+ w .* d_next))
            psi = half_phase .* (expH0 * (half_phase .* psi))
        end
        d_prev = d_next
    end

    return [psi[:,i] for i in 1:size(psi,2)]
end

function compute_fidelity_magspac_ramp(pulses,parameters_dictionary,setup)
    # QuOCS hands the pulse over as a (possibly complex-typed) array on the half-step grid
    spacing_values = Float64.(real.(collect(pulses[1])))
    h = Float64(parameters_dictionary["dt"]) / 2
    final_states = evolve_magspac_ramp(spacing_values,h,setup)
    return real(groundstate_manifold_fidelity(final_states,setup[2]))
end

compute_fidelity_magspac_ramp(pulses,parameters_dictionary) = compute_fidelity_magspac_ramp(pulses,parameters_dictionary,setup_magspac_ramp(parameters_dictionary))

# the ZVD spacing ramp used as the initial guess: B(t) = Bf*[1/4,3/4,1] switching at 0, pi/w,
# 2pi/w, and a(t) its closed-form response a0 + da*sum_{t_i <= t} A_i (1 - cos w(t - t_i)),
# sampled on the half-step grid pulse_ramp expects (ceil(ramptime/(dt/2)) + 1 samples)
function zvd_spacing_pulse(a_start::Float64,a_end::Float64,trap_frequency::Float64,ramptime::Float64,dt::Float64)
    amps = [0.25, 0.5, 0.25]
    timps = [0.0, 1.0, 2.0] .* (pi/trap_frequency)
    @assert timps[end] <= ramptime "ZVD pulse settles at $(timps[end]), after the end of the ramp $(ramptime)"
    times = (0:Int(ceil(ramptime/(dt/2)))) .* (dt/2)
    return [a_start + (a_end - a_start)*sum((A*(1 - cos(trap_frequency*(t - ti))) for (A,ti) in zip(amps,timps) if ti <= t); init=0.0) for t in times]
end

"fin"

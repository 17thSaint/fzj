#####################################################
#=

This file contains the setup function used to optimize, with QuOCS' AD (automatic
differentiation) algorithm, the same interaction strength ramp that
intstren-ramp-control-functions.jl optimizes with dCRAB:

    ULR manifold (U = intstren_start) --(interaction_strength ramp)--> FCI manifold (U = intstren_end)

QuOCS' ADAlgorithm differentiates straight through the figure of merit with jax.grad,
so the FoM itself must be pure JAX (see config_intstrenRamp_AD.py). Only the
pulse-independent constants are computed here in Julia, once at setup:

    - the starting/target ground-state manifolds from ED
    - the decomposition H(u) = H_hop + u * H_int of the ramped Hamiltonian

The decomposition comes from the saved undressed matrices in reading-hamiltonian.jl:
H_hop is the dressed hopping (interaction-independent) and H_int is the interaction
dressed at the largest strength the ramp can reach (the pulse's upper limit) and
divided back out, which is exact because dressInteraction is linear in the overall
strength (U = u * profile) for every scaling type in long_range_scaling. Dressing at
the maximum rather than at unit strength matters for decaying profiles: a far-distance
coupling of an exp/rydberg/dd profile can sit below interaction_cutoff at u = 1 (and
would be silently dropped from H_int) while being above it at the strongly interacting
end of the ramp. No buildHam call is needed, so setup stays cheap at large Hilbert
space dimensions. The hopping is returned in COO form and the (diagonal) interaction
as its diagonal, ready for a sparse/dense JAX propagation.

The split is verified through the eigenstate residuals ||(H_hop + u*H_int)psi - E*psi||
of the ED manifolds at both ramp endpoints, which checks the dressed matrices against
whatever Hamiltonian run_normal_ed actually diagonalized.

Depends on:
    control-functions.jl (for the execute-ed.jl include chain, which already brings
        in reading-hamiltonian.jl's getHopping / getInteraction via two-dimensions.jl)
    intstren-ramp-control-functions.jl (the dCRAB FoM this mirrors)

=#
######################################################

function setup_intstren_ramp_ad(parameters_dictionary)

    Lx::Int = Int(parameters_dictionary["Lx"])
    Ly::Int = Int(parameters_dictionary["Ly"])
    N::Int = Int(parameters_dictionary["N"])
    speccount::Int = Int(parameters_dictionary["speccount"])
    intstren_start::Float64 = Float64(parameters_dictionary["intstren_start"])
    intstren_end::Float64 = Float64(parameters_dictionary["intstren_end"])

    common_params = Dict{String,Any}(
        "output_level"=>0,"Lx"=>Lx,"Ly"=>Ly,"N"=>N,"lr"=>parameters_dictionary["lr"],
        "if_periodic_x"=>parameters_dictionary["if_periodic_x"],"if_periodic_y"=>parameters_dictionary["if_periodic_y"],
        "hopping_anisotropy"=>1.0,"filling"=>0.5,"nev"=>speccount,"if_find_data"=>false,"if_save_data"=>false,
    )

    # forward the interaction profile settings so non-flat (exp / rydberg / dd / gaussian)
    # scalings reach get_normal_model_params_ed instead of silently falling back to flat
    for k in ("scaling_type","corr_length","sigma","blockade_radius","magnetic_spacing")
        haskey(parameters_dictionary,k) && (common_params[k] = parameters_dictionary[k])
    end

    # strongly interacting ULR starting manifold
    pdict_starting = merge(common_params,Dict("interaction_strength"=>intstren_start))
    states_starting,nrgs_starting,_,_,_,lattice_params,hamilt_params = run_normal_ed(pdict_starting; output_level=0)

    # target FCI ground-state manifold
    pdict_ending = merge(common_params,Dict("interaction_strength"=>intstren_end))
    states_ending,nrgs_ending,_,_,_,_,_ = run_normal_ed(pdict_ending; output_level=0)

    # these terms would not scale with the interaction strength, so they cannot be
    # baked into the H_hop + u * H_int split (pinning would even end up inside H_int
    # via dressInteraction and wrongly scale with the ramp)
    get(hamilt_params,"if_pinning",false) && error("if_pinning is not supported by the AD setup")
    get(hamilt_params,"disorder_strength",0.0) != 0.0 && error("disorder is not supported by the AD setup")
    get(hamilt_params,"periodic_potential_strength",0.0) != 0.0 && error("periodic potential is not supported by the AD setup")

    # interaction-independent part: the dressed hopping from the saved undressed matrix
    ham_hop = getHopping(lattice_params,hamilt_params; output_level=0)

    # interaction part: dress at the largest strength the ramp can reach and divide it
    # back out (linear in the overall strength, so u * H_int reproduces the dressing at
    # strength u); dressing at the maximum keeps decaying-profile couplings that would
    # fall below interaction_cutoff at small strengths
    intstren_max::Float64 = Float64(get(parameters_dictionary,"intstren_max",max(intstren_start,intstren_end)))
    intstren_max > 0.0 || error("intstren_max must be positive")
    hp_max = copy(hamilt_params)
    lr_dist::Int = get(hamilt_params,"lr_dist",length(hamilt_params["U"])-1)
    hp_max["interaction_strength"] = intstren_max
    hp_max["U"] = long_range_scaling(lr_dist,Ly,intstren_max;
        scaling=hp_max["scaling_type"],corr_length=hp_max["corr_length"],
        sigma=hp_max["sigma"],blockade_radius=hp_max["blockade_radius"],
        magnetic_spacing=hp_max["magnetic_spacing"])
    ham_int = getInteraction(lattice_params,hp_max; output_level=0) ./ intstren_max

    isdiag(ham_int) || error("interaction matrix is expected to be diagonal")
    hint_diag = Vector{ComplexF64}(diag(ham_int))

    # verify the split against the Hamiltonians run_normal_ed diagonalized, through
    # the eigenstate residuals at both ramp endpoints
    max_residual = 0.0
    for (u,states,nrgs) in ((intstren_start,states_starting,nrgs_starting),(intstren_end,states_ending,nrgs_ending))
        for i in 1:speccount
            psi = Vector{ComplexF64}(states[i])
            residual = norm(ham_hop * psi + (u .* hint_diag) .* psi - nrgs[i] .* psi)
            max_residual = max(max_residual,residual)
        end
    end

    starting_states = reduce(hcat,[Vector{ComplexF64}(states_starting[i]) for i in 1:speccount])
    target_states = reduce(hcat,[Vector{ComplexF64}(states_ending[i]) for i in 1:speccount])

    # hopping in COO form (1-based indices) for the Python side to assemble
    hop_rows,hop_cols,hop_vals = findnz(ham_hop)

    return hop_rows,hop_cols,hop_vals,hint_diag,size(ham_hop,1),starting_states,target_states,max_residual
end

"fin"

#####################################################
#=

This file contains tests of the time evolution code

Depends on:
    execute-ed.jl

=#
######################################################

include("execute-ed.jl")
include("time-evolution.jl")
include("control-functions.jl")
include("../other-funcs/basic-2d-plottings.jl")
include("plottings.jl")

using Test


#= testing the proper overlap with the manifold as a fidelity measure
if false
    lx,ly,n = 4,4,2
    intstren = 10.0
    num_samples = 100

    start_xi = 50.0
    pdict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("scaling_type","exp"),("corr_length",start_xi),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])

    all_states_gs1 = []
    all_states_gs2 = []

    for i in 1:num_samples
        states1,nrgs1,rhos1,filepath1,if_found1,lattice_params1,hamilt_params1 = run_normal_ed(pdict; output_level=0)
        push!(all_states_gs1,states1[1])
        push!(all_states_gs2,states1[2])
    end

    all_fidelities = zeros(Float64,num_samples,num_samples)
    for i1 in 1:num_samples
        for i2 in i1+1:num_samples
            fid_mat = zeros(ComplexF64,2,2)
            fid_mat[1,1] = adjoint(all_states_gs1[i1]) * all_states_gs1[i2]
            fid_mat[1,2] = adjoint(all_states_gs1[i1]) * all_states_gs2[i2]
            fid_mat[2,1] = adjoint(all_states_gs2[i1]) * all_states_gs1[i2]
            fid_mat[2,2] = adjoint(all_states_gs2[i1]) * all_states_gs2[i2]
            all_fidelities[i1,i2] = 0.5 * tr(adjoint(fid_mat)*fid_mat)
        end
    end

    fig = figure()
    imshow(all_fidelities,origin="lower",vmin=0,vmax=1)
    colorbar()
    xlabel("Sample index")
    ylabel("Sample index")
    title("Fidelity between different samples of the same state")
end=#

#= 4x4 N=2 ULR check
if false
    lx,ly,n = 4,4,2
    
    intstrens = range(0.0,10.0,length=20)

    for intstren in intstrens
        params_dict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])

        states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)

        if intstren == 0.0
            scatter(intstren,nrgs[1] - nrgs[1],c="b",label="E0")
            scatter(intstren,nrgs[2] - nrgs[1],c="g",label="E1")
            scatter(intstren,nrgs[3] - nrgs[1],c="k",label="E2")
            for i in 4:length(nrgs)
                scatter(intstren,nrgs[i] - nrgs[1],c="k")
            end
        else
            scatter(intstren,nrgs[1] - nrgs[1],c="b")
            scatter(intstren,nrgs[2] - nrgs[1],c="g")
            for i in 3:length(nrgs)
                scatter(intstren,nrgs[i] - nrgs[1],c="k")
            end
        end

        xlabel("Interaction strength")
        ylabel("Energy - E0")
        title("$(lx)x$(ly) N=$(n) Spectrum")
        legend()
    end
end=#

#= 4x4 N=2 twisting check
if false
    lx,ly,n = 4,4,2
    
    twist_angles = range(0.0,1.0,length=11)

    gaps02 = zeros(Float64,length(twist_angles),length(twist_angles))
    gaps12 = zeros(Float64,length(twist_angles),length(twist_angles))
    for (i1,tw1) in enumerate(twist_angles)
        for (i2,tw2) in enumerate(twist_angles)
            params_dict = Dict([("output_level",1),("Lx",lx),("Ly",ly),("N",n),("tw1",tw1),("tw2",tw2),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",0.0),("filling",0.5),("nev",20),("if_find_data",false),("if_save_data",false)])

            states,nrgs,rhos,filepath,if_found,lattice_params,hamilt_params = run_normal_ed(params_dict; output_level=1)

            gaps02[i1,i2] = nrgs[3] - nrgs[1]
            gaps12[i1,i2] = nrgs[3] - nrgs[2]
        end
    end

    fig = figure()
    imshow(gaps02,origin="lower",extent=(0,1,0,1),vmin=0)
    colorbar()
    xlabel(L"\theta_x/ 2\pi")
    ylabel(L"\theta_y/ 2\pi")
    title("$(lx)x$(ly) N=$(n) Gap E2-E0")

    fig = figure()
    imshow(gaps12,origin="lower",extent=(0,1,0,1),vmin=0)
    colorbar()
    xlabel(L"\theta_x/ 2\pi")
    ylabel(L"\theta_y/ 2\pi")
    title("$(lx)x$(ly) N=$(n) Gap E2-E1")

end=#

#= 4x4 N=2 dt benchmark
if false

    lx,ly,n = 6,3,3
    
    intstren = 0.0

    #tmax_global = 10.0
    #dt_global = 0.001
    dataloc = get_folder_location("cluster-data/exact-diag/time-evo/dt-benchmark")

    start_tx = 1e-8

    params_dict_starting = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("tx",start_tx),("ty",1.0),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",5),("if_find_data",false),("if_save_data",false)])

    states_i,nrgs_i,rhos_i,filepath_i,if_found_i,lattice_params_i,hamilt_params_i = run_normal_ed(params_dict_starting; output_level=0)
    println("Found starting state")

    psi0_1 = states_i[1]
    psi0_2 = states_i[2]

    end_tx = 1.0

    params_dict_ending = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("tx",end_tx),("ty",1.0),("lr","all"),("if_periodic_x",true),("if_periodic_y",true),("hopping_anisotropy",1.0),("interaction_strength",intstren),("filling",0.5),("nev",5),("if_find_data",false),("if_save_data",false)])

    states_f,nrgs_f,rhos_f,filepath_f,if_found_f,lattice_params_f,hamilt_params_f = run_normal_ed(params_dict_ending; output_level=0)
    println("Found ending state")

    psif_1 = states_f[1]
    psif_2 = states_f[2]

    inst_fid = groundstate_manifold_fidelity([psif_1,psif_2],[psi0_1,psi0_2])
    println("Instantaneous fidelity at start: $(inst_fid)")

    speccount = 2
    time_running_args = (nev=speccount,output_level=1,if_instant_gs=false,if_save_data=false,dataloc=dataloc,)

    tmax_global = 0.5
    ramptime = 0.1
    dts = 10 .^ range(-3,-1,length=21)
    for dt_global in dts
        tevo_params = Dict([ ("tx",(linear_ramp,start_tx,end_tx,ramptime)),("dt",dt_global),("tmax",tmax_global) ])
        tevo_gs,tevo_dict,intspec,saving_args = run_timeevo([psi0_1,psi0_2],tevo_params,lattice_params_i,hamilt_params_i; time_running_args...)

        final_gs_manifold = [tevo_gs[1][:,end-1],tevo_gs[2][:,end-1]]
        final_fidelity = groundstate_manifold_fidelity(states_f[1:speccount],final_gs_manifold)
        println("Final fidelity for dt=$(dt_global): $(final_fidelity)")
        scatter(dt_global,final_fidelity,c="b")
    end
    xscale("log")
    xlabel("Time step dt")
    ylabel("Fidelity with target manifold")
    title("Fidelity vs dt $(lx)x$(ly) N=$(n) U=$(intstren) ramp tx ramptime=$(ramptime)")

end=#

# Sums the occupancy on all sites not in allowed_coords.
# get_occupancy returns an Ly×Lx matrix where entry [y,x] is the expected
# particle count at coordinate (x,y), so occ[y,x] for each (x,y) in allowed.
function occ_outside_allowed(state, lp, allowed_coords)
    occ = get_occupancy(state, lp; if_plot=false)
    inside = sum(occ[y, x] for (x, y) in allowed_coords)
    return sum(occ) - inside
end

# Shared parameter template for all position_state tests.
# PBC in both directions is required so the magnetic flux alpha = N/(0.5*Lx*Ly)
# stays below the hard limit of 0.4 checked by check_fluxes. Lattice sizes are
# chosen so that 2N is divisible by Lx or Ly (flux quantization condition on the torus).
# No interactions and no saving to disk; output silenced.
base_pdict(lx, ly, n, nev) = Dict(
    "Lx"=>lx, "Ly"=>ly, "N"=>n, "nev"=>nev,
    "if_periodic_x"=>true, "if_periodic_y"=>true,
    "interaction_strength"=>0.0, "filling"=>0.5,
    "if_save_data"=>false, "if_find_data"=>false, "output_level"=>0
)

# --- two particles restricted to 2x2 corner of 4x4 lattice ---
# The 4 allowed sites {(1,1),(2,1),(1,2),(2,2)} span a 2x2 block.
# With barrier_strength=1e6 on all other sites, the N low-energy eigenstates
# should live almost entirely in the restricted binomial(4,2)=6 dimensional subspace.
# We request nev=4 states (all well below the barrier) and check that each has
# negligible occupancy (<1e-3) outside the allowed block.
if false
    @testset "position_state: N=2 confined to 2x2 corner on 4x4" begin
        pdict = base_pdict(4, 4, 2, 4)
        allowed = [(1,1), (2,1), (1,2), (2,2)]
        states, nrgs, lp, _ = position_state(allowed, pdict; output_level=0)

        for state in states
            @test occ_outside_allowed(state, lp, allowed) < 1e-3
        end
        # all 4 nev energies should be far below the 1e6 barrier; the true
        # sub-barrier spectrum has binomial(4,2)=6 states so nev=4 is safe
        @test all(nrgs .< 1e4)
    end
end

# --- tuple and linear-index dispatch must produce identical eigenstates ---
# position_state accepts either Vector{Int} (linear indices) or
# Vector{Tuple{Int,Int}} (coordinates). Both should build the same Hamiltonian
# and therefore return the same eigenvalues and eigenvectors (up to global phase).
# Linear index convention on a 4x4 lattice: index = (x-1)*Ly + y, so
# (1,1)→1, (2,1)→2, (1,2)→5.
if false
    @testset "position_state: tuple vs linear-index dispatch agree on 4x4 N=2" begin
        pdict = base_pdict(4, 4, 2, 3)
        allowed_lin = [1, 2, 5]
        allowed_tup = [(1, 1), (2, 1), (1, 2)]

        states_lin, nrgs_lin, _, _ = position_state(allowed_lin, pdict; output_level=0)
        states_tup, nrgs_tup, _, _ = position_state(allowed_tup, pdict; output_level=0)

        # eigenvalues must match exactly (same H constructed both ways)
        @test isapprox(nrgs_lin, nrgs_tup, atol=1e-8)
        # |⟨ψ_lin|ψ_tup⟩| = 1 confirms same ground state up to global phase
        @test isapprox(abs(dot(states_lin[1], states_tup[1])), 1.0, atol=1e-6)
    end
end

# --- two particles confined to bottom row of a 6x2 lattice ---
# All 6 sites in y=1 are allowed; the y=2 row is barriered.
# binomial(6,2)=15 restricted states, so nev=6 is comfortably in the sub-barrier sector.
# We directly sum the occupancy over the top row rather than using occ_outside_allowed
# to demonstrate the row-indexed check explicitly.
# Lattice choice: Lx=6, Ly=2 gives alpha=1/3 and 2N=4 divisible by Ly=2 ✓
if false
    @testset "position_state: N=2 confined to bottom strip on 6x2" begin
        pdict = base_pdict(6, 2, 2, 6)
        allowed = [(x, 1) for x in 1:6]
        states, nrgs, lp, _ = position_state(allowed, pdict; output_level=0)

        for state in states
            occ = get_occupancy(state, lp; if_plot=false)
            # occ[2,:] is the top row (y=2); all weight should be in y=1
            @test sum(occ[2, :]) < 1e-3
        end
        @test all(nrgs .< 1e4)
    end
end

# --- three particles restricted to a 3x2 block on 4x6 lattice ---
# 6 allowed sites → binomial(6,3)=20 restricted states; nev=6 requested.
# Lattice choice: Lx=4, Ly=6 gives alpha=0.25 and 2N=6 divisible by Ly=6 ✓
if false
    @testset "position_state: N=3 confined to 3x2 block on 4x6" begin
        pdict = base_pdict(4, 6, 3, 6)
        allowed = [(x, y) for x in 1:3, y in 1:2] |> vec
        states, nrgs, lp, _ = position_state(allowed, pdict; output_level=0)

        for state in states
            @test occ_outside_allowed(state, lp, allowed) < 1e-3
        end
        @test all(nrgs .< 1e4)
    end
end

# --- four particles restricted to left half of 6x4 lattice ---
# x in 1:3, y in 1:4 gives 12 allowed sites → binomial(12,4)=495 restricted states.
# Lattice choice: Lx=6, Ly=4 gives alpha=1/3 and 2N=8 divisible by Ly=4 ✓
if false
    @testset "position_state: N=4 confined to left half on 6x4" begin
        pdict = base_pdict(6, 4, 4, 8)
        allowed = [(x, y) for x in 1:3, y in 1:4] |> vec
        states, nrgs, lp, _ = position_state(allowed, pdict; output_level=0)

        for state in states
            @test occ_outside_allowed(state, lp, allowed) < 1e-3
        end
        @test all(nrgs .< 1e4)
    end
end

# --- three particles in bottom strip of 6x3 lattice ---
# Only y=1 is allowed; rows y=2 and y=3 are barriered.
# binomial(6,3)=20 restricted states, nev=8 requested.
# Lattice choice: Lx=6, Ly=3 gives alpha=1/3 and 2N=6 divisible by both Lx=6 and Ly=3 ✓
if false
    @testset "position_state: N=3 confined to bottom strip on 6x3" begin
        pdict = base_pdict(6, 3, 3, 8)
        allowed = [(x, 1) for x in 1:6]
        states, nrgs, lp, _ = position_state(allowed, pdict; output_level=0)

        for state in states
            occ = get_occupancy(state, lp; if_plot=false)
            # occ[2:3,:] covers both non-allowed rows; sum should be near zero
            @test sum(occ[2:3, :]) < 1e-3
        end
        @test all(nrgs .< 1e4)
    end
end






























"fin"
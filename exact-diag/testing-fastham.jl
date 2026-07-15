#####################################################
#=

This file is for testing a way of building more efficient ED Hamiltonians

Depends on:
    exact-diag/execute-ed.jl
    exact-diag/reading-hamiltonian.jl

=#
######################################################

include("execute-ed.jl")
include("reading-hamiltonian.jl")
using Test

do_all::Bool = true

if do_all || false
@testset "Only Hopping with Anis/Fluxes" begin;

    configs = [(4,4,2),(4,3,2),(4,4,3),(3,2,2),(4,3,3)]
    tw1 = 0.0
    tw2 = 0.0
    #tws = range(0.0,1.0,length=5)
    hanises = [1.0,0.5,2.0]
    for if_periodic_x in [true,false]
        for if_periodic_y in [true,false]
            for config in configs
                lx,ly,n = config
                for hanis in hanises
                    #println("Working on Lx = ",lx," Ly = ",ly," N = ",n," tw1 = ",tw1," tw2 = ",tw2," hanis = ",hanis)
                    pdict = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("tw1",tw1),("tw2",tw2),("if_check_fluxes",false),("if_pinning",false),("if_periodic_x",if_periodic_x),("if_periodic_y",if_periodic_y),("hopping_anisotropy",hanis),("interaction_strength",0.0),("lr",0),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
                    lattice_params,hamilt_params,running_args = get_normal_model_params_ed(pdict)    
                    basis_dataloc = running_args.basis_dataloc
                    full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
                    lattice_params["full_basis"] = full_basis
                    ham_correct = buildHam(lattice_params,hamilt_params; output_level=0)

                    ham_new = buildHopping(lattice_params; output_level=0,if_save=false)
                    ham_new = dressHopping(hamilt_params,lattice_params,ham_new; output_level=0)

                    @test isapprox(ham_correct,ham_new,atol=1e-10)
                end
            end
        end
    end

end;
end

if do_all || false
@testset "Interaction Part" begin;
    configs = [(4,4,2),(4,3,2),(4,4,3),(4,3,3)]
    tw1 = 0.0
    tw2 = 0.0
    #tws = range(0.0,1.0,length=5)
    hanises = [1.0,0.5,2.0]
    intstrens = [1.0,2.0,3.0]

    for if_periodic_x in [true,false]
        for if_periodic_y in [true,false]
            for config in configs
            lx,ly,n = config
                for hanis in hanises
                    for intstren in intstrens
                        #println("Working on Lx = ",lx," Ly = ",ly," N = ",n," tw1 = ",tw1," tw2 = ",tw2," hanis = ",hanis)
                        pdict = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("tw1",tw1),("tw2",tw2),("if_check_fluxes",false),("if_pinning",false),("if_periodic_x",if_periodic_x),("if_periodic_y",if_periodic_y),("hopping_anisotropy",hanis),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
                        lattice_params,hamilt_params,running_args = get_normal_model_params_ed(pdict)    
                        basis_dataloc = running_args.basis_dataloc
                        full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
                        lattice_params["full_basis"] = full_basis
                        ham_correct = buildHam(lattice_params,hamilt_params; output_level=0)

                        ham_new = buildHopping(lattice_params; output_level=0,if_save=false)
                        ham_new = dressHopping(hamilt_params,lattice_params,ham_new; output_level=0)
                        intham = buildInteraction(lattice_params; output_level=0,if_save=false)
                        intham = dressInteraction(hamilt_params,lattice_params,intham; output_level=0)
                        ham_new .+= intham

                        @test isapprox(ham_correct,ham_new,atol=1e-10)
                    end
                end
            end
        end
    end
end;
end

if do_all || true
@testset "High level getHamiltonian function (with pinning)" begin;
    lx,ly,n = 4,4,2
    tw1 = 0.0
    tw2 = 0.0
    hanis = 1.0
    intstren = 1.0
    if_periodic_x = true
    if_periodic_y = true
    pdict = Dict([("output_level",0),("Lx",lx),("Ly",ly),("N",n),("tw1",tw1),("tw2",tw2),("if_pinning",true),("if_periodic_x",if_periodic_x),("if_periodic_y",if_periodic_y),("hopping_anisotropy",hanis),("interaction_strength",intstren),("lr","all"),("filling",0.5),("nev",10),("if_find_data",false),("if_save_data",false)])
    lattice_params,hamilt_params,running_args = get_normal_model_params_ed(pdict)
    basis_dataloc = running_args.basis_dataloc
    full_basis = n_particle_basis(lattice_params; output_level=running_args.output_level,dataloc=basis_dataloc)
    lattice_params["full_basis"] = full_basis
    ham_correct = buildHam(lattice_params,hamilt_params; output_level=0)
    ham_new = getHamiltonian(lattice_params,hamilt_params; output_level=0)

    #println("Checking Hamiltonians: ",isapprox(ham_correct,ham_new,atol=1e-10))
    @test isapprox(ham_correct,ham_new,atol=1e-10)
end;
end

#=if do_all || false
@testset "Twistings" begin;
    
end;
end=#




#####################################################
#=

Tests checking that the dressed-matrix Hamiltonians from reading-hamiltonian.jl
(getHamiltonian, used by find_eigenstates when if_reading=true) match the
usual ED construction (buildHam/applyHam)

run test_reading_hamiltonian_all() for everything

=#
######################################################

function make_reading_test_params(lx::Int64,ly::Int64,n::Int64; kwargs...)
    pdict = Dict{String,Any}([("output_level",0),
                    ("Lx",lx),("Ly",ly),("N",n),
                    ("tw1",0.0),("tw2",0.0), # getHamiltonian errors on finite twists
                    ("if_check_fluxes",false),
                    ("if_periodic_x",get(kwargs,:if_periodic_x,true)),
                    ("if_periodic_y",get(kwargs,:if_periodic_y,true)),
                    ("hopping_anisotropy",get(kwargs,:hopping_anisotropy,1.0)),
                    ("interaction_strength",get(kwargs,:interaction_strength,1.0)),
                    ("scaling_type",get(kwargs,:scaling_type,"flat")),
                    ("if_pinning",get(kwargs,:if_pinning,false)),
                    ("lr","all"),
                    ("filling",get(kwargs,:filling,0.5)),
                    ("nev",1),
                    ("if_find_data",false),("if_save_data",false)])

    lattice_params,hamilt_params,running_args = get_normal_model_params_ed(pdict)
    lattice_params["full_basis"] = n_particle_basis(lattice_params; output_level=0,dataloc=running_args.basis_dataloc)
    return lattice_params,hamilt_params
end

function sparse_max_deviation(m1::SparseMatrixCSC{ComplexF64},m2::SparseMatrixCSC{ComplexF64})
    difference = m1 - m2
    return nnz(difference) == 0 ? 0.0 : maximum(abs.(nonzeros(difference)))
end

function test_reading_ham_matrix(; kwargs...)
    output_level::Int64 = get(kwargs,:output_level,1)
    atol::Float64 = get(kwargs,:atol,1e-10)

    # empty dataloc so getHamiltonian always builds fresh instead of pulling saved matrices
    empty_dataloc::String = mktempdir()

    # label => (lx,ly,n,pbc_x,pbc_y,hanis,intstren,scaling,if_pinning,filling)
    cases = [("hopping only, no flux",(4,4,2,true,true,1.0,0.0,"flat",false,0.0)),
             ("hopping + flux",(4,4,2,true,true,1.0,0.0,"flat",false,0.5)),
             ("anisotropy + flat U",(4,4,2,true,true,2.0,1.0,"flat",false,0.5)),
             ("odd Ly + flat U",(4,3,2,true,true,0.5,2.5,"flat",false,0.5)),
             ("flat U + pinning",(4,4,3,true,true,1.0,1.0,"flat",true,0.5)),
             ("open boundaries",(4,4,2,false,false,1.0,1.0,"flat",false,0.5)),
             ("cylinder",(4,4,2,true,false,1.0,1.0,"flat",false,0.5)),
             ("periodic Ly=2 + flux",(3,2,2,true,true,1.0,0.0,"flat",false,0.5)),
             ("exp-scaled U",(4,4,2,true,true,1.0,1.0,"exp",false,0.5)),
             ("gaussian-scaled U",(4,4,2,true,true,1.0,1.0,"gaussian",false,0.5))]

    if_all_passed::Bool = true
    for (label,(lx,ly,n,px,py,hanis,intstren,scaling,pinning,filling)) in cases
        lattice_params,hamilt_params = make_reading_test_params(lx,ly,n; if_periodic_x=px,if_periodic_y=py,hopping_anisotropy=hanis,interaction_strength=intstren,scaling_type=scaling,if_pinning=pinning,filling=filling)

        ham_usual = buildHam(lattice_params,hamilt_params; output_level=0)
        ham_reading = getHamiltonian(lattice_params,hamilt_params; output_level=0,if_save=false,dataloc=empty_dataloc)

        deviation = sparse_max_deviation(ham_usual,ham_reading)
        if_passed = deviation < atol
        if_all_passed &= if_passed
        output_level > 0 && println(rpad(label,22)," ",lx,"x",ly," N=",n," : ",if_passed ? "PASS" : "FAIL"," (max deviation ",round(deviation,sigdigits=3),")")
    end

    return if_all_passed
end

function test_reading_ham_cache(; kwargs...)
    output_level::Int64 = get(kwargs,:output_level,1)
    atol::Float64 = get(kwargs,:atol,1e-10)

    # save the undressed matrices to a scratch folder so the real cache is untouched
    tmp_dataloc::String = mktempdir()
    lattice_params,hamilt_params = make_reading_test_params(4,3,2; interaction_strength=1.5)

    hopping = buildHopping(lattice_params; output_level=0,if_save=false)
    interaction = buildInteraction(lattice_params; output_level=0,if_save=false)
    saveHopping(copy(hopping),tmp_dataloc,lattice_params; output_level=0)
    saveInteraction(copy(interaction),tmp_dataloc,lattice_params; output_level=0)

    # make sure getHamiltonian will actually pull from disk below
    if_found_hopping::Bool = findHopping(lattice_params; dataloc=tmp_dataloc,output_level=0)[1]
    if_found_interaction::Bool = findInteraction(lattice_params; dataloc=tmp_dataloc,output_level=0)[1]
    if !(if_found_hopping && if_found_interaction)
        output_level > 0 && println("cache roundtrip 4x3 N=2 : FAIL (saved matrices not found on disk)")
        return false
    end

    # dressHopping/dressInteraction mutate their input, so dress copies
    ham_fresh = dressHopping(hamilt_params,lattice_params,copy(hopping); output_level=0) + dressInteraction(hamilt_params,lattice_params,copy(interaction); output_level=0)
    ham_from_disk = getHamiltonian(lattice_params,hamilt_params; output_level=0,if_save=false,dataloc=tmp_dataloc)

    deviation = sparse_max_deviation(ham_fresh,ham_from_disk)
    if_passed::Bool = deviation < atol
    output_level > 0 && println("cache roundtrip 4x3 N=2 : ",if_passed ? "PASS" : "FAIL"," (max deviation ",round(deviation,sigdigits=3),")")
    return if_passed
end

function test_reading_ham_spectrum(; kwargs...)
    output_level::Int64 = get(kwargs,:output_level,1)
    nev::Int64 = get(kwargs,:nev,6)
    nrg_atol::Float64 = get(kwargs,:nrg_atol,1e-8)
    residual_atol::Float64 = get(kwargs,:residual_atol,1e-6)

    lattice_params,hamilt_params = make_reading_test_params(4,4,2; interaction_strength=1.0)
    shared_args = (if_densmat=false,if_save_data=false,if_exact=false,output_level=0)

    # NOTE: the if_reading path pulls from (and on first run saves to) the real
    # cluster-data/exact-diag cache since find_eigenstates hardcodes that location
    states_usual,nrgs_usual,_,ham_usual = find_eigenstates(nev,lattice_params,hamilt_params; shared_args...,if_function=false,if_reading=false)
    states_reading,nrgs_reading,_,_ = find_eigenstates(nev,lattice_params,hamilt_params; shared_args...,if_function=false,if_reading=true)
    states_func,nrgs_func,_ = find_eigenstates(nev,lattice_params,hamilt_params; shared_args...,if_function=true,if_reading=false)

    if_all_passed::Bool = true
    for (label,nrgs_other,states_other) in [("reading vs buildHam",nrgs_reading,states_reading),("function vs buildHam",nrgs_func,states_func)]
        nrg_deviation = maximum(abs.(nrgs_other .- nrgs_usual))

        # exact degeneracies make state-by-state overlaps ill-defined, so instead check every
        # state is an eigenstate of the usual Hamiltonian (catches any gauge/matrix mismatch);
        # also allow conjugated states since buildHam stores the transpose of the applyHam operator
        residual = maximum([norm(ham_usual*states_other[i] .- nrgs_other[i].*states_other[i]) for i in 1:nev])
        conj_residual = maximum([norm(ham_usual*conj.(states_other[i]) .- nrgs_other[i].*conj.(states_other[i])) for i in 1:nev])

        if_passed = nrg_deviation < nrg_atol && min(residual,conj_residual) < residual_atol
        if_all_passed &= if_passed
        output_level > 0 && println(rpad(label,22)," 4x4 N=2 nev=",nev," : ",if_passed ? "PASS" : "FAIL"," (nrg deviation ",round(nrg_deviation,sigdigits=3),", eigenstate residual ",round(residual,sigdigits=3),")")
        output_level > 0 && residual > residual_atol && conj_residual < residual_atol && println("    NOTE: states only match after complex conjugation (buildHam stores the transposed operator)")
    end

    return if_all_passed
end

function test_reading_hamiltonian_all(; kwargs...)
    output_level::Int64 = get(kwargs,:output_level,1)

    output_level > 0 && println("--- getHamiltonian vs buildHam matrix elements ---")
    if_matrix::Bool = test_reading_ham_matrix(; kwargs...)
    output_level > 0 && println("--- save/load roundtrip of undressed matrices ---")
    if_cache::Bool = test_reading_ham_cache(; kwargs...)
    output_level > 0 && println("--- find_eigenstates if_reading vs usual paths ---")
    if_spectrum::Bool = test_reading_ham_spectrum(; kwargs...)

    if_all_passed::Bool = if_matrix && if_cache && if_spectrum
    output_level > 0 && println("All reading-hamiltonian tests ",if_all_passed ? "PASSED" : "FAILED")
    return if_all_passed
end













































"fin"
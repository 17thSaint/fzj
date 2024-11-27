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


















































"fin"
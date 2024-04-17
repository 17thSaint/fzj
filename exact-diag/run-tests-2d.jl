using Test
include("two-dimensions.jl")

if_all = true

#=if false || if_all
    @testset "linear index" begin
        Lx,Ly = 4,4
        @test linear_index((1,1),Lx,Ly) == 1
        @test linear_index((2,1),Lx,Ly) == 2
        @test linear_index((1,2),Lx,Ly) == 5
        @test linear_index((2,2),Lx,Ly) == 6
        @test linear_index((4,4),Lx,Ly) == 16
        @test linear_index((1,4),Lx,Ly) == 13
    end
    
    @testset "coordinate index" begin
        Lx,Ly = 4,4
        @test coordinate(1,Lx,Ly) == (1,1)
        @test coordinate(2,Lx,Ly) == (2,1)
        @test coordinate(5,Lx,Ly) == (1,2)
        @test coordinate(6,Lx,Ly) == (2,2)
        @test coordinate(16,Lx,Ly) == (4,4)
        @test coordinate(13,Lx,Ly) == (1,4)
    end
end

if false || if_all
    @testset "hopping matrices exact" begin
        Lx,Ly = 2,2
        site = (1,1)

        m_right_x = right_x_hopping_matrix(site,Lx,Ly)
        @test m_right_x == [0 0 0 0; 1 -1 0 0; 0 0 1 0; 0 0 0 1]

        m_left_x = left_x_hopping_matrix(site,Lx,Ly)
        @test m_left_x == [0 0 0 0; 1 -1 0 0; 0 0 1 0; 0 0 0 1]

        m_up_y = up_y_hopping_matrix(site,Lx,Ly)
        @test m_up_y == [0 0 0 0; 0 1 0 0; 1 0 -1 0; 0 0 0 1]

        m_down_y = down_y_hopping_matrix(site,Lx,Ly)
        @test m_down_y == [0 0 0 0; 0 1 0 0; 1 0 -1 0; 0 0 0 1]
    end

    @testset "hopping matrices operation" begin
        Lx,Ly = 6,6
        site = (6,3)
        config_vector = zeros(Int64,Lx*Ly)
        config_vector[linear_index(site,Lx,Ly)] = 1

        right_x_config_vector = zeros(Int64,Lx*Ly)
        right_x_config_vector[linear_index((1,3),Lx,Ly)] = 1
        m_right_x = right_x_hopping_matrix(site,Lx,Ly)
        @test m_right_x * config_vector == right_x_config_vector

        left_x_config_vector = zeros(Int64,Lx*Ly)
        left_x_config_vector[linear_index((5,3),Lx,Ly)] = 1
        m_left_x = left_x_hopping_matrix(site,Lx,Ly)
        @test m_left_x * config_vector == left_x_config_vector

        up_y_config_vector = zeros(Int64,Lx*Ly)
        up_y_config_vector[linear_index((6,4),Lx,Ly)] = 1
        m_up_y = up_y_hopping_matrix(site,Lx,Ly)
        @test m_up_y * config_vector == up_y_config_vector

        down_y_config_vector = zeros(Int64,Lx*Ly)
        down_y_config_vector[linear_index((6,2),Lx,Ly)] = 1
        m_down_y = down_y_hopping_matrix(site,Lx,Ly)
        @test m_down_y * config_vector == down_y_config_vector
    end

    @testset "hop onto another particle" begin
        Lx,Ly = 6,6
        site = (3,3)
        config_vector = zeros(Int64,Lx*Ly)
        config_vector[linear_index(site,Lx,Ly)] = 1
        config_vector[linear_index((4,3),Lx,Ly)] = 1

        m_right_x = x_hopping_matrix(site,1,Lx,Ly)
        @test m_right_x * config_vector == zeros(Int64,Lx*Ly)
    end
end

if false || if_all
    @testset "number operator" begin
        Lx,Ly = 5,5
        site = (2,2)
        m = number_operator_matrix(site,Lx,Ly)
        
        single_particle_config_vector = zeros(Int64,Lx*Ly)
        single_particle_config_vector[linear_index(site,Lx,Ly)] = 1
        @test m * single_particle_config_vector == single_particle_config_vector

        two_particle_config_vector = zeros(Float64,Lx*Ly)
        two_particle_config_vector[linear_index(site,Lx,Ly)] = 1/sqrt(2)
        two_particle_config_vector[linear_index((3,3),Lx,Ly)] = 1/sqrt(2)
        @test isapprox(transpose(two_particle_config_vector) * m * two_particle_config_vector,0.5,atol=1e-3)
    end
end=#

if false || if_all
    @testset "operator change of basis old with naive" begin
        Lx,Ly = 4,4
        N = 2
        if_periodic_x,if_periodic_y = false,false
        full_basis,basis_dict = generate_basis_naive(Lx,Ly,N; output_level=0)
        lattice_params::Dict{String,Any} = Dict("Lx"=>Lx,
                              "Ly"=>Ly,
                              "N"=>N,
                              "if_periodic_x"=>if_periodic_x,
                              "if_periodic_y"=>if_periodic_y,
                              "full_basis"=>full_basis,
                              "basis_dict"=>basis_dict)

        change_of_basis,change_of_basis_inv = change_of_basis_matrix(lattice_params)
        lattice_params["change_of_basis"] = change_of_basis
        lattice_params["change_of_basis_inv"] = change_of_basis_inv
        alpha = 0.0
        hamilt_params = Dict("alpha"=>alpha,
                             "tx"=>1.0,
                             "ty"=>1.0,
                             "U"=>zeros(Ly),
                             "interaction_cutoff"=>1e-5)
        
        H = buildHam_naive(lattice_params,hamilt_params; output_level=0)
        nev = 3
        rez = eigsolve(H,nev)
        gs = rez[2][findfirst(x->x==minimum(rez[1]),rez[1])]
        
        # test if the hopping amplitude is the same between old and new basis
        s1 = (1,1)
        other_site = (3,2)
        old_basis_gs = change_of_basis * gs
        for s2 in [other_site,s1]
            new_basis_hopping = hopping_matrix_naive(s1,s2,lattice_params)
            new_basis_expval = abs(conj(transpose(gs)) * new_basis_hopping * gs)
            
            old_basis_hopping = get_old_basis_version([0 0; 1 0],[0 1; 0 0],s1,s2,Lx,Ly)
            old_basis_expval = abs(conj(transpose(old_basis_gs)) * old_basis_hopping * old_basis_gs)

            @test isapprox(new_basis_expval,old_basis_expval,atol=1e-5)
        end

        # test if the occupancy is the same between old and new basis
        old_basis_occupancy = get_old_basis_version([0 0;0 1],other_site,Lx,Ly)
        old_basis_occ_expval = abs(conj(transpose(old_basis_gs)) * old_basis_occupancy * old_basis_gs)
        new_basis_occupancy = change_of_basis_inv * old_basis_occupancy * change_of_basis
        new_basis_occ_expval = abs(conj(transpose(gs)) * new_basis_occupancy * gs)
        @test isapprox(new_basis_occ_expval,old_basis_occ_expval,atol=1e-5)


        # test if efficient occupancy matches the full basis method
        rho = density_matrix_naive(gs,lattice_params)
        occs_densmat = get_occupancy(rho,lattice_params; if_plot=false, plot_title="Density Matrix Occupancy")
        efficient_occs = get_occupancy_naive(gs,lattice_params; if_plot=false, plot_title="Efficient Occupancy")
        @test all((efficient_occs .- occs_densmat) ./ occs_densmat .< 1e-1)
    
    end
end

if false || if_all
    @testset "hopping probability function old way" begin
        Lx,Ly = 4,4
        N = 2
        if_periodic_x,if_periodic_y = false,false
        full_basis,basis_dict = generate_basis_naive(Lx,Ly,N; output_level=0)
        lattice_params::Dict{String,Any} = Dict("Lx"=>Lx,
                              "Ly"=>Ly,
                              "N"=>N,
                              "if_periodic_x"=>if_periodic_x,
                              "if_periodic_y"=>if_periodic_y,
                              "full_basis"=>full_basis,
                              "basis_dict"=>basis_dict)
        change_of_basis,change_of_basis_inv = change_of_basis_matrix(lattice_params)
        lattice_params["change_of_basis"] = change_of_basis
        lattice_params["change_of_basis_inv"] = change_of_basis_inv

        alpha = 0.0
        hamilt_params = Dict("alpha"=>alpha,
                             "tx"=>1.0,
                             "ty"=>1.0,
                             "U"=>zeros(Ly),
                             "interaction_cutoff"=>1e-5)
        
        H = buildHam_naive(lattice_params,hamilt_params; output_level=0)
        nev = 3
        rez = eigsolve(H,nev)
        gs = rez[2][findfirst(x->x==minimum(rez[1]),rez[1])]

        # test that occupancy matrix matches the hopping probability on itself for every site
        direct_occupancy = get_occupancy_naive(gs,lattice_params; if_plot=false, plot_title="Direct Occupancy")
        handmade_occs = zeros(Float64,Ly,Lx)
        for j in 1:Lx
            for s in 1:Ly
                site = (j,s)
                handmade_occs[s,j] = hopping_probability_old(gs,site,site,lattice_params)
            end
        end
        percent_diff = abs.((direct_occupancy .- handmade_occs) ./ direct_occupancy)
        # less than 3% difference
        @test all(percent_diff .< 0.03)

        # test hopping probability makes same density matrix
        rho = density_matrix_naive(gs,lattice_params)
        rho_efficient = density_matrix_old(gs,lattice_params; output_level=0)
        @test isapprox(rho,rho_efficient,atol=1e-5)
    end
end


if false || if_all
    @testset "New faster Index numbering" begin
        @test find_basis_index([4,3,2]) == 4
        @test find_basis_index([3,2,1]) == 1
        @test find_basis_index([4,2,1]) == 2
        @test find_basis_index([4,3,1]) == 3

        @test find_basis_index([2,1]) == 1
        @test find_basis_index([3,1]) == 2
        @test find_basis_index([3,2]) == 3
        @test find_basis_index([4,1]) == 4

        @test find_basis_index([4,3,2,1]) == 1
    end
end

if false || if_all
    @testset "Generating new basis" begin
        n2_basis = n_particle_basis(2,2,2; output_level=0)
        @test n2_basis == [2 3 3 4 4 4; 1 1 2 1 2 3]
        for i in 1:size(n2_basis,2)
            @test find_basis_index(n2_basis[:,i]) == i
        end

        n3_basis = n_particle_basis(3,2,2; output_level=0)
        @test n3_basis == [3 4 4 4; 2 2 3 3; 1 1 1 2]
        for i in 1:size(n3_basis,2)
            @test find_basis_index(n3_basis[:,i]) == i
        end

        n4_basis = n_particle_basis(4,2,2; output_level=0)
        @test n4_basis == [4; 3; 2; 1;;]

    end
end

if false || if_all
    @testset "hopping probability function" begin
        Lx,Ly = 4,4
        N = 2
        if_periodic_x,if_periodic_y = false,false
        full_basis = n_particle_basis(N,Lx,Ly; output_level=0)
        lattice_params::Dict{String,Any} = Dict("Lx"=>Lx,
                              "Ly"=>Ly,
                              "N"=>N,
                              "if_periodic_x"=>if_periodic_x,
                              "if_periodic_y"=>if_periodic_y,
                              "full_basis"=>full_basis)
    

        alpha = 0.0
        hamilt_params = Dict("alpha"=>alpha,
                             "tx"=>1.0,
                             "ty"=>1.0,
                             "U"=>zeros(Ly),
                             "interaction_cutoff"=>1e-5)
        
        H = buildHam(lattice_params,hamilt_params; output_level=0)
        nev = 3
        rez = eigsolve(H,nev)
        gs = rez[2][findfirst(x->x==minimum(rez[1]),rez[1])]

        # test that occupancy matrix matches the hopping probability on itself for every site
        direct_occupancy = get_occupancy(gs,lattice_params; if_plot=false, plot_title="Direct Occupancy")
        handmade_occs = zeros(Float64,Ly,Lx)
        for j in 1:Lx
            for s in 1:Ly
                site = (j,s)
                handmade_occs[s,j] = hopping_probability(gs,site,site,lattice_params; output_level=0)
            end
        end
        percent_diff = abs.((direct_occupancy .- handmade_occs) ./ direct_occupancy)
        # less than 0.1% difference
        @test all(percent_diff .< 0.001)
        
    end
end

if false || if_all
    @testset "faster density matrix calculation OBC" begin
        Lx,Ly = 4,4
        N = 2
        if_periodic_x,if_periodic_y = false,false
        full_basis = n_particle_basis(N,Lx,Ly; output_level=0)
        lattice_params::Dict{String,Any} = Dict("Lx"=>Lx,
                              "Ly"=>Ly,
                              "N"=>N,
                              "if_periodic_x"=>if_periodic_x,
                              "if_periodic_y"=>if_periodic_y,
                              "full_basis"=>full_basis)
    

        alpha = 0.0
        hamilt_params = Dict("alpha"=>alpha,
                             "tx"=>1.0,
                             "ty"=>1.0,
                             "U"=>zeros(Ly),
                             "interaction_cutoff"=>1e-5)
        
        H = buildHam(lattice_params,hamilt_params; output_level=0)
        nev = 3
        rez = eigsolve(H,nev)
        gs = rez[2][findfirst(x->x==minimum(rez[1]),rez[1])]

        rho_slow = density_matrix_slow(gs,lattice_params; output_level=0)
        rho_fast = density_matrix(gs,lattice_params; output_level=0)

        @test isapprox(rho_slow,rho_fast,atol=1e-5)
    end
end

if false || if_all
    @testset "faster density matrix calculation Periodic BC" begin
        Lx,Ly = 4,4
        N = 2
        if_periodic_x,if_periodic_y = true,true
        full_basis = n_particle_basis(N,Lx,Ly; output_level=0)
        lattice_params::Dict{String,Any} = Dict("Lx"=>Lx,
                              "Ly"=>Ly,
                              "N"=>N,
                              "if_periodic_x"=>if_periodic_x,
                              "if_periodic_y"=>if_periodic_y,
                              "full_basis"=>full_basis)
    

        alpha = 0.0
        hamilt_params = Dict("alpha"=>alpha,
                             "tx"=>1.0,
                             "ty"=>1.0,
                             "U"=>zeros(Ly),
                             "interaction_cutoff"=>1e-5)
        
        H = buildHam(lattice_params,hamilt_params; output_level=0)
        nev = 3
        rez = eigsolve(H,nev)
        gs = rez[2][findfirst(x->x==minimum(rez[1]),rez[1])]

        rho_slow = density_matrix_slow(gs,lattice_params; output_level=0)
        rho_fast = density_matrix(gs,lattice_params; output_level=0)

        @test isapprox(rho_slow,rho_fast,atol=1e-5)
    end
end
































"fin"
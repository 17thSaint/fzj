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

if true || if_all
    @testset "operator change of basis" begin
        Lx,Ly = 2,2
        N = 2
        if_periodic_x,if_periodic_y = true,true
        full_basis,basis_dict = generate_basis(Lx,Ly,N; output_level=0)
        
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
        
        H = buildHam(lattice_params,hamilt_params; output_level=0)
        
        nev = 3
        rez = eigsolve(H,nev)
        gs = rez[2][findfirst(x->x==minimum(rez[1]),rez[1])]
        
        new_basis_hopping = hopping_matrix((1,1),(2,1),lattice_params)
        new_basis_expval = abs(conj(transpose(gs)) * new_basis_hopping * gs)
        
        old_basis_hopping = get_old_basis_version([0 1; 0 0],[0 0; 1 0],(1,1),(2,1),Lx,Ly)
        old_basis_gs = change_of_basis * gs
        old_basis_expval = abs(conj(transpose(old_basis_gs)) * old_basis_hopping * old_basis_gs)

        @test isapprox(new_basis_expval,old_basis_expval,atol=1e-3)

        
        
        efficient_occupancy = get_occupancy(gs,lattice_params; if_plot=false)[1,1]

        old_basis_occupancy = get_old_basis_version([0 0;0 1],(1,1),Lx,Ly)
        old_basis_occ_expval = abs(conj(transpose(old_basis_gs)) * old_basis_occupancy * old_basis_gs)
        #println("Old Basis Occupancy: ",old_basis_occ_expval)
        @test isapprox(efficient_occupancy,old_basis_occ_expval,atol=1e-5)

        new_basis_occupancy = change_of_basis_inv * old_basis_occupancy * change_of_basis
        new_basis_occ_expval = abs(conj(transpose(gs)) * new_basis_occupancy * gs)
        #println("New Basis Occupancy: ",new_basis_occ_expval)
        @test isapprox(efficient_occupancy,new_basis_occ_expval,atol=1e-5)
    
        
    end
end

































"fin"
using Test
include("laughlin-wavefunc.jl")

# testing without excitation added
if true
@testset "Jain-Kamila close to Girvin-Jach" begin;

particles = 5
con = start_rand_config(particles,3)

jkver = laughlin_wavefunction(con)[1]
gjver = laughlin_wavefunction_girvinjach(con)[1]

@test abs(real(jkver) - real(gjver)) / abs(real(gjver)) < 0.001

end;
end

if true
@testset "1QH exact for 3 particles" begin;

particles = 3
con = start_rand_config(particles,3)

jkver = laughlin_wavefunction(con; qe_loc = 0.0, qe_cutoff = particles)[1]
gjver = laughlin_wavefunction_girvinjach(con; qe_loc = 0.0)[1]

exactmat = zeros(particles,particles) .* im
for i in 1:particles
	for j in 1:particles
		exactmat[i,j] = jastrow(con,con[i]) * (con[i]^(j))
	end
end
exver = log(Complex(det(exactmat)))
exver += -0.25*sum(abs2.(con))


@test abs(real(jkver) - real(exver)) / abs(real(exver)) < 0.001
@test abs(real(gjver) - real(exver)) / abs(real(exver)) < 0.001

end;
end

# this test is a quasiparticle which is different from the quasihole
# it will not work with the current version of laughlin_wavefunction function 
if false
@testset "1QP exact for 3 particles" begin;

particles = 3
con = start_rand_config(particles,3)

jkver = laughlin_wavefunction(con; qe_loc = 0.0, qe_cutoff = particles)[1]

exactmat = zeros(particles,particles) .* im
dj = zeros(3) .* im
dj[1] = 2*(con[1] - con[3]) + 2*(con[1] - con[2])
dj[2] = 2*(con[2] - con[1]) + 2*(con[2] - con[3])
dj[3] = 2*(con[3] - con[2]) + 2*(con[3] - con[1])
exactmat[1,:] = dj
exactmat[2,:] = [2*jastrow(con,con[i]) + con[i]*dj[i] for i in 1:particles]
exactmat[3,:] = [2*2*con[i]*jastrow(con,con[i]) + (con[i]^2)*dj[i] for i in 1:particles]
exver = log(Complex(det(exactmat)))
exver += -0.25*sum(abs2.(con))


@test abs(real(jkver) - real(exver)) / abs(real(exver)) < 0.001

end;
end































"fin"

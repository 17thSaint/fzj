using Test
include("fqh-thesis/cf-wavefunc.jl")
include("reverse-flux.jl")


if true
@testset "derivatives" begin

#particles = 3
for particles in [3,5,10,20,50]
con = start_rand_config(particles,3)

autodiff_result = nth_deriv_j2(con,con[1],1)
exact_result = 2*sum([(con[1] - con[k])*prod([(con[1] - con[j])^2 for j in deleteat!([l for l in 2:length(con)],k-1)]) for k in 2:length(con)])
@test abs(autodiff_result-exact_result)/abs(exact_result) <= 0.001

end

end;
end


if true
@testset "analytical Jain-Kamila" begin

particles = 3
con = start_rand_config(particles,3)

jkexact = test_3parts_jainkamila(con)[1]
myver = reverse_flux_wavefunction(con)[1]
# check if the analytical Jain-Kamila for 3 particles matches my AutoDiff calculation
@test abs(real(jkexact) - real(myver))/abs(real(jkexact)) <= 0.001

# the MSc function has some factors of 2 missing so this checks the edited
# analytical Jain-Kamila function against the MSc result
# this shows that the new AutoDiff matches the MSc up to a few factors of 2
exact_with_correction = test_3parts_jainkamila_msc(con)

all_deriv_orders = get_deriv_orders_matrix(particles)
all_pascal = [get_pascals_triangle(i)[2] for i in 1:particles]
acc_sets_matrix = get_full_acc_matrix(particles)
msc = get_rf_wavefunc(con,acc_sets_matrix,all_pascal,all_deriv_orders,[0,[0]],true)

@test abs(real(exact_with_correction) - real(msc))/abs(real(exact_with_correction)) <= 0.001

end;
end

if true
@testset "quasielectron" begin

particles = 3
con = start_rand_config(particles,3)

# testing that the first derivative with the quasielectron matches the exact
firstderiv_zj2 = jastrow_squared(con,con[1]) + con[1]*2*((con[1] - con[2])*((con[1] - con[3])^2) + (con[1] - con[3])*((con[1] - con[2])^2))
myver = nth_deriv_j2(con, con[1], 1; if_qe = true)
@test abs(real(firstderiv_zj2) - real(myver))/abs(real(firstderiv_zj2)) <= 0.0001

# check that my version with quasielectron matches the exact result
for c in 1:3
	jkexact = test_3parts_jainkamila(con,3; qe_cutoff = c)[1]
	myver = reverse_flux_wavefunction(con,3; qe_cutoff = c)[1]
	@test abs(real(jkexact) - real(myver))/abs(real(jkexact)) <= 0.001
end


end;
end





































"fin"

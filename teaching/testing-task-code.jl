using Test

include("task-code.jl")

do_all = false

count = 2
org = [0 for i in 1:count]
n = length(org)
js = 1.0
hz = 1.0
hx = 1.0
dt = 0.5
sites = [1,2]

if do_all | false
@testset "single qubit many-body matrix" begin

	dub_args = (xx,sites,js,dt)
	dmat = build_matrix_from_elements(get_exp_two_qubit_elem,dub_args,2)

	sing_args = (x,sites[1],js,dt)
	smat = build_matrix_from_elements(get_exp_single_qubit_elem,sing_args,2)

end;
end

if do_all | false
@testset "elementwise multiply matrices" begin

	args1 = (x,sites[1],hx,dt)
	args2 = (x,sites[2],hx,dt)
	fullmat1 = build_matrix_from_elements(get_exp_single_qubit_elem,args1,count)
	fullmat2 = build_matrix_from_elements(get_exp_single_qubit_elem,args2,count)
	correct_result = fullmat1 * fullmat2

	calced_result = build_matrix_from_elements(elem_mult_matrices,(get_exp_single_qubit_elem,args1,get_exp_single_qubit_elem,args2),count)
	@test calced_result == correct_result
end;
end

if do_all | false
@testset "local hamiltonian parts" begin

	og_fullz, ogxx, ogzz = get_exp_zpart(sites,n,js,hz,dt)
	calced_xx = build_matrix_from_elements(get_xx_elem,(sites,js,dt),count)
	calced_zz = build_matrix_from_elements(get_zz_elem,(sites,hz,dt),count)
	calced_fullz = build_matrix_from_elements(get_fullz_elem,(sites,js,hz,dt),count)
	
	@test isapprox(calced_xx,ogxx,atol=10^-5)
	@test isapprox(calced_zz,ogzz,atol=10^-5)
	@test isapprox(calced_fullz,og_fullz,atol=10^-5)
	
	ogx = get_exp_xpart(sites[1],count,dt,hx)
	calced_x = build_matrix_from_elements(get_x_elem,(sites[1],hx,dt),count)
	@test isapprox(calced_x,ogx,atol=10^-5)
	
	localham_args = (sites,js,hx,hz,dt)
	calced_localham = build_matrix_from_elements(get_localham_elem,localham_args,count)
	oglocalham = get_mbham_local(n,sites,js,hz,hx,dt)
	@test isapprox(calced_localham,oglocalham,atol=10^-5)

end;
end

if do_all | false
@testset "full hamiltonian" begin

	for counts in 2:3
		ogham = get_full_ham(counts,js,hz,hx,dt)
	
		calced_ham = build_matrix_from_elements(get_ham_elem_finished,(counts,js,hx,hz,dt),counts)
		@test isapprox(ogham,calced_ham,atol=10^-5)
		@test isapprox(conj(transpose(ogham)) * ogham,I,atol=10^-5)
	end

end;
end

if do_all | false
@testset "expectation value" begin
	
	downdown_wavefunc = get_matwavefunc_givenorg(org)
	op_args = (id,1)
	calced_result = get_expectation(get_single_qubit_elem,op_args,downdown_wavefunc)
	@test isapprox(calced_result,1.0,atol=10^-5)
	
end;
end

if true

num_sites = 2
org = [0 for i in 1:num_sites]
org[Int(ceil(num_sites/2))] = 1
js = 5.0
hx = 0.0
hz = 0.25
input_args = (num_sites,js,hx,hz)
time_end = 0.5
timesteps = 100
starting_wavefunc = get_matwavefunc_givenorg(org)
rez = do_time_evolution(time_end,timesteps,starting_wavefunc; arguments=input_args,corr=true)
#plot_correlations(rez["corrs"])
plot([rez["corrs"][i][1] for i in 1:timesteps+1])
#plot_site_mag_time_ev([i for i in 1:num_sites],rez["local_mag"],10 .* real.(rez["times"][1]),num_sites,timesteps)



end



























"fin"

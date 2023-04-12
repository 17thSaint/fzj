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
#ogham = get_mbham_local(n,sites,js,hz,hx,dt)

if do_all | false
@testset "single qubit many-body matrix" begin

	dub_args = (xx,sites,js,dt)
	dmat = build_matrix_from_elements(get_exp_two_qubit_elem,dub_args,2)

	sing_args = (x,sites[1],js,dt)
	smat = build_matrix_from_elements(get_exp_single_qubit_elem,sing_args,2)

end;
end

if do_all | true
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


#oginter = get_exp_zpart(sites,n,js,hz,dt)
#newinter = im.*zeros(2^n,2^n)
#zz = im.*zeros(2^n,2^n)
































"fin"

using Test

include("task-code.jl")

do_all = false

count = 3
org = [0 for i in 1:count]
n = length(org)
js = 1.0
hz = 1.0
hx = 1.0
dt = 0.5
sites = [1,2]
#ogham = get_mbham_local(n,sites,js,hz,hx,dt)

if do_all | false
@testset "multiplying matrices get element" begin

	corr_multip = get_exp_xpart(sites[1],n,dt,hx) * get_exp_xpart(sites[1],n,dt,hx)
	multip = im.*zeros(2^n,2^n)
	for i in 1:2^n
		for j in 1:2^n
			bs = int_to_binary(i-1,n)
			bps = int_to_binary(j-1,n)
			multip[i,j] = elem_mult_matrices(get_exp_single_qubit_elem,get_exp_single_qubit_elem,bs,bps,sites[1],sites[1],hx,hx,dt,x,x)
		end
	end
	@test multip == corr_multip

end;
end


#oginter = get_exp_zpart(sites,n,js,hz,dt)
#newinter = im.*zeros(2^n,2^n)
#zz = im.*zeros(2^n,2^n)
































"fin"

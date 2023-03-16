using Test

include("proper-mps.jl")

do_all = true

function make_upupup()
	d1 = Index(1)
	d2 = Index(2)
	d3 = Index(2)
	d4 = Index(2)
	n_L = Index(2)
	n_R = Index(2)
	n_C = Index(2)

	left_tensor = ITensor(d1,d2,n_L)
	left_tensor[1,1,1] = 0
	right_tensor = ITensor(d3,d1,n_R)
	right_tensor[1,1,1] = 0
	center_tensor = ITensor(d2,d3,n_C)
	center_tensor[1] = 0
	left_tensor[:,:,1] = [0; 0]
	left_tensor[:,:,2] = [0; 1]
	center_tensor[:,:,1] = [0 0; 0 0]
	center_tensor[:,:,2] = [0 0; 0 1]
	right_tensor[:,:,2] = [0; 1]
	right_tensor[:,:,1] = [0; 0]

	rez_amp_tensor = left_tensor * center_tensor * right_tensor
	
	upupup_wavefunc = get_full_wavefunc(rez_amp_tensor,3)[1]
	
	return upupup_wavefunc
end

if do_all | false
@testset "upupup state" begin
	
	calced_upupup_wavefunc = make_upupup()
	
	mat_part1 = [1 0; 0 0]
	mat_part2 = [0 0; 0 0]
	for i in 1:2
		for j in 1:2
			@test calced_upupup_wavefunc[i,j,1] == mat_part1[i,j]
			@test calced_upupup_wavefunc[i,j,2] == mat_part2[i,j]
		end
	end

end;
end

if do_all | false
@testset "2 site onsite ham" begin

calced_onsite_ham = make_onsite_ham(2)

creat,annih = get_creat_annih_tens()
iden = make_identity_tensor()
local_annih1 = replaceind(annih,inds(annih)[1],inds(creat)[2])
comp1 = (creat * local_annih1) * iden
local_annih2 = replaceind(annih,inds(annih)[1],inds(creat)[2])
comp2 = iden * (creat * local_annih1)
replaceinds!(comp2,inds(comp2),inds(comp1))
theory_ham = comp1 + comp2

for i in 1:prod(size(theory_ham))
	@test theory_ham[i] == calced_onsite_ham[i]
end

ct_h1 = conj(transpose(calced_onsite_ham))
replaceinds!(ct_h1,inds(ct_h1),inds(calced_onsite_ham))
@test ct_h1 == calced_onsite_ham

end;
end

if do_all | false
@testset "2 site neighbor ham" begin

calced_nn_ham = make_nearest_neighbor_ham(2)
calced_nn_ham_periodic = make_nearest_neighbor_ham(2,true)

creat,annih = get_creat_annih_tens()
comp1 = annih * creat
periodic_comp = creat * annih

theory_nn_ham = comp1
replaceinds!(comp1,inds(comp1),inds(periodic_comp))
theory_nn_ham_periodic = comp1 + periodic_comp

for i in 1:prod(size(theory_nn_ham))
	@test theory_nn_ham[i] == calced_nn_ham[i]
	@test theory_nn_ham_periodic[i] == calced_nn_ham_periodic[i]
end

ct_h = conj(transpose(calced_nn_ham))
replaceinds!(ct_h,inds(ct_h),inds(calced_nn_ham))
@test ct_h == calced_nn_ham

ct_hp = conj(transpose(calced_nn_ham_periodic))
replaceinds!(ct_hp,inds(ct_hp),inds(calced_nn_ham_periodic))
@test ct_hp == calced_nn_ham_periodic

end;
end

if do_all | false
@testset "3 site onsite energy" begin

onsite_ham_n3 = make_onsite_ham(3)
wavefunc_n3 = make_upupup()
calced_nrg = get_expect_ham_val(onsite_ham_n3,wavefunc_n3,3)
@test calced_nrg == 3.0

end;
end

if do_all | false
@testset "svd returns same wavefunc" begin

num_sites = 5
num_states = 3
rand_wavefunc = make_random_wavefunc(num_sites,num_states)
as,cs = make_As(rand_wavefunc,num_sites,num_states,[true])
remade_wavefunc = prod(as) * cs[end]

for i in 1:prod(size(rand_wavefunc))
	@test isapprox(imag(remade_wavefunc[i]),imag(rand_wavefunc[i]),atol=10^-10)
	@test isapprox(real(remade_wavefunc[i]),real(rand_wavefunc[i]),atol=10^-10)
end

end;
end










"fin"

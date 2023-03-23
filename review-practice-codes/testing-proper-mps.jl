using Test

include("proper-mps.jl")

do_all = false

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

if do_all | true
@testset "svd returns same wavefunc" begin

num_sites = 3
num_states = 3
keep_all = true
rand_wavefunc = make_random_wavefunc(num_sites,num_states)
as,cs,dims = make_As(rand_wavefunc,num_sites,num_states,keep_all,false)

for j in 1:num_sites-1
	remade_wavefunc = prod(as[1:j]) * cs[j+1]
	for i in 1:prod(size(cs[1]))
		@test isapprox(remade_wavefunc[i],cs[1][i],atol=10^-5)
	end
end
end;
end

if do_all | false
@testset "mats to wavefunc to mats rebuild" begin

num_states = 3
num_sites = 2
keep_all = true
all_amps,each_tensor = get_random_mps_coeffs(num_sites,num_states)
wavefunc_from_rand_matrices, local_coeff = get_full_wavefunc(all_amps,num_sites)
all_mats,all_cs,dims = make_As(wavefunc_from_rand_matrices,num_sites,num_states,keep_all,false)
rebuilt_wavefunc = prod(all_mats) * all_cs[end]

for i in 1:prod(size(rebuilt_wavefunc))
	@test isapprox(imag(rebuilt_wavefunc[i]),imag(wavefunc_from_rand_matrices[i]),atol=10^-5)
	@test isapprox(real(rebuilt_wavefunc[i]),real(wavefunc_from_rand_matrices[i]),atol=10^-5)
end

end;
end

if do_all | false

num_sites = 4
num_states = 2
keeping = "all"
keep_type = get_keeping_type(keeping)
rand_wavefunc = make_random_wavefunc(num_sites,num_states)
as,cs,dims = make_As(rand_wavefunc,num_sites,num_states,keep_type,true)

@testset "orthogonality of A1" begin

for a1 in 1:size(as[1])[2]
	for a1p in 1:size(as[1])[2]
		local_val = sum(as[1][i,a1] * conj(as[1])[i,a1p] for i in 1:size(as[1])[1])
		if a1 == a1p
			@test isapprox(1.0,local_val,atol=10^-5)
		else
			@test isapprox(0.0,local_val,atol=10^-5)
		end
	end
end

end;

@testset "orthogonality of As" begin

for a2 in 1:size(as[2])[2]
	for a2p in 1:size(as[2])[2]
		local_val = sum(as[2][i,a2] * conj(as[2])[i,a2p] for i in 1:size(as[2])[1])
		if a2 == a2p
			@test isapprox(1.0,local_val,atol=10^-5)
		else
			@test isapprox(0.0,local_val,atol=10^-5)
		end
	end
end


end;

@testset "check norm" begin

inner_prod = (cs[end] * conj(cs[end]))[1]
@test isapprox(1.0,inner_prod,atol=10^-5)

end;

end

if do_all | false

num_sites = 3
num_states = 4
keeping = "count"
keep_count = 3
keep_type = get_keeping_type(keeping,keep_count)
rand_wavefunc = make_random_wavefunc(num_sites,num_states)
as,cs,dims = make_As(rand_wavefunc,num_sites,num_states,keep_type,true)

@testset "orthogonality of A1 w/SVD" begin

for a1 in 1:size(as[1])[2]
	for a1p in 1:size(as[1])[2]
		local_val = sum(as[1][i,a1] * conj(as[1])[i,a1p] for i in 1:size(as[1])[1])
		if a1 == a1p
			@test isapprox(1.0,local_val,atol=10^-5)
		else
			@test isapprox(0.0,local_val,atol=10^-5)
		end
	end
end

end;

@testset "orthogonality of A2 w/SVD" begin

for a2 in 1:size(as[2])[2]
	for a2p in 1:size(as[2])[2]
		local_val = sum(as[2][i,a2] * conj(as[2])[i,a2p] for i in 1:size(as[2])[1])
		if a2 == a2p
			@test isapprox(1.0,local_val,atol=10^-5)
		else
			@test isapprox(0.0,local_val,atol=10^-5)
		end
	end
end

end;

@testset "check norm" begin

inner_prod = (cs[end] * conj(cs[end]))[1]
@test isapprox(1.0,inner_prod,atol=10^-5)

end;

end







"fin"

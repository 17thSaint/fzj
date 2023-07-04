include("fqh_effective.jl")
include("long-range-ttn.jl")


params_dict = Dict([("L",4),("nflavors",4),("nbosons",4),("alpha",1/8),("if_2ord_pert",true)])
L = get(params_dict, "L", 4)
nflavors = get(params_dict, "nflavors", Int(L/2))
nbosons = get(params_dict, "nbosons", nflavors)
t1 = get(params_dict, "t1", 1.0)
t2 = get(params_dict, "t2", 1.0)
U = get(params_dict, "U", 100)
U1 = 4*t1^2/U
U2 = U1/2
conserve_qns = true
if_nn_int = true
if_2ord_pert = true
nsweeps = 5
mdim = get(params_dict, "mdim", 100)
noise = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0]
if_save_data = false

#
change = 0.0001
count = 1
alpha_start  = 1/8
alpha_end = 1/8
derivs = zeros(3,count)
errs = zeros(3,count)
alphas = count == 1 ? [alpha_start] : [alpha_start + (i-1)*(alpha_end-alpha_start)/(count-1) for i in 1:count]
all_wavefuncs = []
println(alphas)

phi  = 2π*alpha_start
model_paras = (t1 = t1, t2 = t2, phi = phi, U1 = U1, U2 = U2, L = L,
					   nflavors = nflavors, nbosons = nbosons, if_nn_int = if_nn_int, if_2ord_pert = if_2ord_pert, mdim = mdim, noise = noise, if_save_data = if_save_data)
other_params_dict = Dict([("nbosons",nbosons),("U",U),("conserve_qns",conserve_qns),("nsweeps",nsweeps),("mdim",mdim),("noise",noise)])
metadata_dict = merge(named_tuple_to_dict(model_paras),other_params_dict)
datafile_name = make_parameters_filename(params_dict)

for i in 1:count
	alpha = alphas[i]
	if true
	#phi  = 2π*alpha
	phi_right  = 2π*(alpha+change)
	phi_left  = 2π*(alpha-change)

	psi = execute_mps(; model_paras...)
	append!(all_wavefuncs,[psi])
	println("Done Center")
	
	#model_paras_left = copy(model_paras)
	#model_paras_left[:phi] = phi_left
	psi_left = execute_mps(; model_paras..., phi = phi_left, psi_guess = psi)
	
	#E_left,psi_left = dmrg(H_left, psi; maxdim = 100, nsweeps = 3, noise = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0])
	
	#E_right,psi_right = dmrg(H_right, psi; maxdim = 100, nsweeps = 3, noise = [1E-2, 1E-2, 1E-2, 1E-2, 1E-2,0])
	
	end
	occ_mat_c = get_occupancy(psi)
	#occ_mat_left = get_occupancy(psi_left)
	#occ_mat_right = get_occupancy(psi_right)
	
	#bd_left = bulk_density(nothing; occ_mat=occ_mat_left)
	#bd_c = bulk_density(nothing; occ_mat=occ_mat_c)
	#bd_right = bulk_density(nothing; occ_mat=occ_mat_right)
	#=println("Bulk Densities = $bd_left, $bd_c, $bd_right")
	for w in 1:2
		deriv12 = deriv_bulk_dens(nothing,nothing,-change,w; occ_mat1=occ_mat_left,occ_mat2=occ_mat_c)
		deriv23 = deriv_bulk_dens(nothing,nothing,-change,w; occ_mat1=occ_mat_c,occ_mat2=occ_mat_right)
		deriv13 = deriv_bulk_dens(nothing,nothing,-2*change,w; occ_mat1=occ_mat_left,occ_mat2=occ_mat_right)
		deriv = mean([deriv12,deriv23,deriv13])
		derr = std([deriv12,deriv23,deriv13])
		derivs[w,i] = deriv
		errs[w,i] = derr
		println("Streda Value = $deriv +/- $derr")
	end
	=#
end
#
rez = get_densdens_corrs(all_wavefuncs[1]; avgs=false,if_plot=true)
#println(derivs,", ",errs)
#for i in 1:2
#errorbar(alphas,derivs[i,:],yerr=[errs[i,:],errs[i,:]],label="$i")
#end
#legend()




























"fin"

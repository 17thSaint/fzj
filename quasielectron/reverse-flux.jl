using SpecialMatrices,LinearAlgebra,LsqFit,TaylorDiff

include("general-funcs.jl")

function nth_deriv_jastrow(z, which, n; kwargs...)
    if_qe = get(kwargs, :if_qe, false)
    qe_loc = get(kwargs, :qe_loc, 0.0+im*0.0)
    if n == 0
    	if if_qe
    		return (which - qe_loc)*jastrow(z, which; kwargs...)
    	else
        	return jastrow(z, which; kwargs...)  # 0th derivative is the function itself
        end
    else
    	#println("Working on $(n-1)-th order")
    	if if_qe
    		#println("using QE at $qe_loc")
    		return derivative(c -> (c - qe_loc)*jastrow(z,c; kwargs...),which,n)
    	else
		#println("No QE")
	        return derivative(c -> jastrow(z,c; kwargs...),which,n)
	end
        #return conj(ForwardDiff.derivative(c -> real(nth_deriv_j2(z, c, n - 1)), which)[1])
    end
end

function test_jastrow_firstderiv(z)
	autodiff_result = nth_deriv_jastrow(z,z[1],1;power=2)
	exact_result = 2*sum([(z[1] - z[k])*prod([(z[1] - z[j])^2 for j in deleteat!([l for l in 2:length(z)],k-1)]) for k in 2:length(z)])
	if abs(autodiff_result-exact_result)/abs(exact_result) <= 0.001
		return true
	else
		return false,autodiff_result,exact_result
	end
end

function reverse_flux_wavefunction(z,m=3; kwargs...)
	num_parts = length(z)
	qe_cutoff = get(kwargs, :qe_cutoff, 0)
	jast_pow = Int((m+1)/2)
	full_mat = im.*zeros(num_parts,num_parts)
	if_qe = qe_cutoff > 0
	for i in 1:num_parts
		#println("Making Row $i")
		if i > qe_cutoff
			if_qe = false
		end
		full_mat[i,:] = [log(Complex(nth_deriv_jastrow(z,z[j],i-1; if_qe = if_qe,power=jast_pow,kwargs...))) + log(2)*(i-1) for j in 1:num_parts]
	end
	
	result = get_log_det(full_mat)[1]
	
	result += -0.25*sum(abs2.(z))
	
	return result,[0,0]
end

function test_3parts_girvinjach(z,m=3)
	num_parts = length(z)
	p = Int((m+1)/2)
	
	dz = im .* zeros(3)
	dz[1] = 2*p*dist_btw(z,2,3,2*p) * (dist_btw(z,1,2,2*p-1)*dist_btw(z,1,3,2*p) + dist_btw(z,1,2,2*p)*dist_btw(z,1,3,2*p-1))
	dz[2] = 2*p*dist_btw(z,1,3,2*p) * (-dist_btw(z,1,2,2*p-1)*dist_btw(z,2,3,2*p) + dist_btw(z,1,2,2*p)*dist_btw(z,2,3,2*p-1))
	dz[3] = 2*p*dist_btw(z,1,2,2*p) * (-dist_btw(z,1,3,2*p-1)*dist_btw(z,2,3,2*p) - dist_btw(z,1,3,2*p)*dist_btw(z,2,3,2*p-1))
	
	result = 0.0
	for j in 1:num_parts
		for i in 1:j-1
			result += log(Complex(-2*(dz[i] - dz[j])))
		end
	end
	#=
	for i in 1:num_parts
		result += -0.25*abs2(z[i])
	end
	=#
	
	return result
end

function analytical_3parts_girvinjach(z)
	lin_result = 0.0
	
	lin_result += -48 * dist_btw(z,1,2,1) * dist_btw(z,1,3,4) * dist_btw(z,2,3,4)
	lin_result += 48 * dist_btw(z,1,2,4) * dist_btw(z,1,3,2) * dist_btw(z,2,3,3)
	lin_result += -48 * dist_btw(z,1,2,2) * dist_btw(z,1,3,3) * dist_btw(z,2,3,4)
	lin_result += -48 * dist_btw(z,1,2,4) * dist_btw(z,1,3,4) * dist_btw(z,2,3,1)
	lin_result += 48 * dist_btw(z,1,2,4) * dist_btw(z,1,3,1) * dist_btw(z,2,3,4)
	lin_result += 48 * dist_btw(z,1,2,2) * dist_btw(z,1,3,4) * dist_btw(z,2,3,3)
	lin_result += 48 * dist_btw(z,1,2,4) * dist_btw(z,1,3,3) * dist_btw(z,2,3,2)
	lin_result += 48 * dist_btw(z,1,2,3) * dist_btw(z,1,3,2) * dist_btw(z,2,3,4)
	lin_result += 48 * dist_btw(z,1,2,3) * dist_btw(z,1,3,4) * dist_btw(z,2,3,2)
	
	result = log(Complex(lin_result)) + log(8)
	#=
	for i in 1:3
		result += -0.25*abs2(z[i])
	end
	=#
	return result
end

function test_3parts_jainkamila(z,m=3; kwargs...)
	qe_cutoff = get(kwargs, :qe_cutoff, 0)
	num_parts = length(z)
	p = Int((m+1)/2)
	
	fulljs = im .* zeros(num_parts,num_parts)
	fulljs[1,1] = dist_btw(z,1,2,1)*dist_btw(z,1,3,1)
	fulljs[1,2] = dist_btw(z,2,1,1)*dist_btw(z,2,3,1)
	fulljs[1,3] = dist_btw(z,3,1,1)*dist_btw(z,3,2,1)
	
	fulljs[2,1] = 1*(dist_btw(z,1,2,1-1)*dist_btw(z,1,3,1) + dist_btw(z,1,2,1)*dist_btw(z,1,3,1-1))
	fulljs[2,2] = 1*(dist_btw(z,2,1,1-1)*dist_btw(z,2,3,1) + dist_btw(z,2,1,1)*dist_btw(z,2,3,1-1))
	fulljs[2,3] = 1*(dist_btw(z,3,1,1-1)*dist_btw(z,3,2,1) + dist_btw(z,3,1,1)*dist_btw(z,3,2,1-1))
	

	fulljs[3,:] = [1,1,1] .* 2
	
	matver = im .* zeros(num_parts,num_parts)
	
	matver[1,:] = (fulljs[1,:]) .^ p
	if qe_cutoff > 0
		matver[1,:] .*= z
	end
	
	matver[2,:] = (2*p) .* (fulljs[1,:] .* fulljs[2,:])
	if qe_cutoff > 1
		matver[2,:] .*= z
		matver[2,:] += 2*(fulljs[1,:]) .^ p
	end
	
	matver[3,:] = (4*p) .* (((fulljs[2,:] .^ 2) ./ 1) .+ (fulljs[1,:] .* fulljs[3,:]))
	if qe_cutoff > 2
		matver[3,:] .*= z
		matver[3,:] += (8*p) .* (fulljs[1,:] .* fulljs[2,:])
	end

	return log(Complex(det(matver))) - 0.25*sum(abs2.(z)),log.(matver)
end

function test_3parts_jainkamila_msc(z,m=3)
	num_parts = length(z)
	p = Int((m+1)/2)
	
	fulljs = im .* zeros(num_parts,num_parts)
	fulljs[1,1] = dist_btw(z,1,2,1)*dist_btw(z,1,3,1)
	fulljs[1,2] = dist_btw(z,2,1,1)*dist_btw(z,2,3,1)
	fulljs[1,3] = dist_btw(z,3,1,1)*dist_btw(z,3,2,1)
	
	fulljs[2,1] = 1*(dist_btw(z,1,2,1-1)*dist_btw(z,1,3,1) + dist_btw(z,1,2,1)*dist_btw(z,1,3,1-1))
	fulljs[2,2] = 1*(dist_btw(z,2,1,1-1)*dist_btw(z,2,3,1) + dist_btw(z,2,1,1)*dist_btw(z,2,3,1-1))
	fulljs[2,3] = 1*(dist_btw(z,3,1,1-1)*dist_btw(z,3,2,1) + dist_btw(z,3,1,1)*dist_btw(z,3,2,1-1))
	

	fulljs[3,:] = [1,1,1] .* 1
	
	matver = im .* zeros(num_parts,num_parts)

	matver[1,:] = (fulljs[1,:]) .^ p
	matver[2,:] = (2*p) .* (fulljs[1,:] .* fulljs[2,:])
	matver[3,:] = (4*p) .* (((fulljs[2,:] .^ 2) ./ 2) .+ (fulljs[1,:] .* fulljs[3,:]))

	return log(Complex(det(matver))) - 0.25*sum(abs2.(z))
end


function test_4parts_girvinjach(z,m=3)
	num_parts = length(z)
	p = Int((m+1)/2)
	
	dz = im .* zeros(4)
	dz[1] = log(Complex(2*p*dist_btw(z,2,3,2*p)*dist_btw(z,2,4,2*p)*dist_btw(z,3,4,2*p))) + log(Complex(dist_btw(z,1,2,2*p-1)*dist_btw(z,1,3,2*p)*dist_btw(z,1,4,2*p) + dist_btw(z,1,2,2*p)*dist_btw(z,1,3,2*p-1)*dist_btw(z,1,4,2*p) + dist_btw(z,1,2,2*p)*dist_btw(z,1,3,2*p)*dist_btw(z,1,4,2*p-1)))
	dz[2] = log(Complex(2*p*dist_btw(z,1,3,2*p)*dist_btw(z,1,4,2*p)*dist_btw(z,3,4,2*p))) + log(Complex(-dist_btw(z,1,2,2*p-1)*dist_btw(z,2,3,2*p)*dist_btw(z,2,4,2*p) + dist_btw(z,1,2,2*p)*dist_btw(z,2,3,2*p-1)*dist_btw(z,2,4,2*p) + dist_btw(z,1,2,2*p)*dist_btw(z,2,3,2*p)*dist_btw(z,2,4,2*p-1)))
	dz[3] = log(Complex(2*p*dist_btw(z,1,2,2*p)*dist_btw(z,1,4,2*p)*dist_btw(z,2,4,2*p))) + log(Complex(-dist_btw(z,1,3,2*p-1)*dist_btw(z,2,3,2*p)*dist_btw(z,3,4,2*p) - dist_btw(z,1,3,2*p)*dist_btw(z,2,3,2*p-1)*dist_btw(z,3,4,2*p) - dist_btw(z,1,3,2*p)*dist_btw(z,2,3,2*p)*dist_btw(z,3,4,2*p-1)))
	dz[4] = log(Complex(2*p*dist_btw(z,1,2,2*p)*dist_btw(z,1,3,2*p)*dist_btw(z,2,3,2*p))) + log(Complex(-dist_btw(z,1,4,2*p-1)*dist_btw(z,2,4,2*p)*dist_btw(z,3,4,2*p) - dist_btw(z,1,4,2*p)*dist_btw(z,2,4,2*p-1)*dist_btw(z,3,4,2*p) - dist_btw(z,1,4,2*p)*dist_btw(z,2,4,2*p)*dist_btw(z,3,4,2*p-1)))
	
	result = 0.0
	for j in 1:num_parts
		for i in 1:j-1
			result += log(Complex(-2)) + get_log_subtract(dz[i],dz[j])
		end
	end
	#
	for i in 1:num_parts
		result += -0.25*abs2(z[i])
	end
	#
	
	return result
end

function linfit_matrix(mat1,mat2; kwargs...)
	direction = get(kwargs, :direction, "col")
	if_plot = get(kwargs, :if_plot, false)
	title_string = get(kwargs, :title_string, "")
	
	linfit_func(x,p) = p[1] .* x .+ p[2]
	direction == "col" ? xs = [mat1[i,:] for i in 1:size(mat1)[1]] : xs = [mat1[:,i] for i in 1:size(mat1)[1]]
	direction == "col" ? ys = [mat2[i,:] for i in 1:size(mat2)[1]] : ys = [mat2[:,i] for i in 1:size(mat2)[1]]
	fitparams = [[] for i in 1:size(mat1)[1]]
	fitparams_error = [[] for i in 1:size(mat1)[1]]
	for i in 1:size(mat1)[1]
		localfit = curve_fit(linfit_func,xs[i],ys[i],[0.5,0.5])
		fitparams[i] = localfit.param
		fitparams_error[i] = stderror(localfit)
	end
	
	if if_plot
		multip = [fitparams[i][1] for i in 1:size(mat1)[1]]
		multip_error = [fitparams_error[i][1] for i in 1:size(mat1)[1]]
		fig = figure()
		errorbar([i for i in 1:size(mat1)[1]],multip,yerr=[multip_error,multip_error])
		title("Multiplicative Part, $title_string")
		
		constant = [fitparams[i][2] for i in 1:size(mat1)[1]]
		constant_error = [fitparams_error[i][2] for i in 1:size(mat1)[1]]
		fig = figure()
		errorbar([i for i in 1:size(mat1)[1]],constant,yerr=[constant_error,constant_error])
		title("Constant Part, $title_string")
	end
	
	return fitparams,fitparams_error
end

#=
include("fqh-thesis/cf-wavefunc.jl")
particles = 3
con = start_rand_config(particles,3)
rm = sqrt(2*particles*3)

newver = reverse_flux_wavefunction(con)
oldver = 
#correct = test_3parts_jainkamila(con)
println("My ver = $newver and Exact = $correct")
=#

#=
using PyPlot
data_count = 50
xs = [-3*rm + i*(2*3*rm)/data_count for i in 1:data_count]
jkrfs = zeros(data_count,data_count)
exacts = zeros(data_count,data_count)
for i in 1:data_count
	new_x = xs[i]
	for j in 1:data_count
		println(i,", ",j)
		new_y = xs[j]
		start_con[1] = new_x - im*new_y
		jkrf = test_3parts_jainkamila(start_con)
		exact = test_3parts_girvinjach(start_con)
		jkrfs[j,i] = real(jkrf)
		exacts[j,i] = real(exact)
	end
end
#
fig = figure()
imshow(jkrfs) #.- oldrfs)
colorbar()
title("JK Exact")
#
fig = figure()
imshow(exacts)
colorbar()
title("GJ Exact")
#

linfit_matrix(exacts,jkrfs; if_plot=true,title_string="GJ and JK")
=#
#=m = 3
for particles in [3,4]
	con = start_rand_config(particles,m)
	
	rm = sqrt(2*particles*m)
	data_count = 50
	xs = [-10*rm + i*(2*10*rm)/data_count for i in 1:data_count]
	myver = fill(0.0,(data_count,data_count))
	for i in 1:length(xs)
		local_x = xs[i]
		for j in 1:length(xs)
			local_y = xs[j]
			#append!(xs_plot,[local_x])
			#append!(ys_plot,[local_y])
			con[1] = local_x - im*local_y
			
			myver[i,j] = 2*real(reverse_flux_wavefunction(con,m)[1])
		end
	end
	
	fig = figure()
	imshow(myver)
	colorbar()
	title("My RF")
	
	#mscver = get_rf_wavefunc(con,allowed_sets_matrix,full_pasc_tri,full_derivs,[0,[0]],true)
	#=if particles == 3
		exactver = test_3parts(con,m)
	elseif particles == 4
		exactver = test_4parts(con,m)
	end
	=#
	#println(isapprox(real(myver),real(exactver),atol=10^-3))
	#println("Mine = $(round(myver,digits=5)), Masters = $(round(mscver,digits=5)), Exact = $(round(exactver,digits=5))")
#end
=#






































"fin"

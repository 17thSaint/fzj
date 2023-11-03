using SpecialMatrices,LinearAlgebra,LsqFit,Zygote

function start_rand_config(num_parts::Int, m::Int)
    # Calculate the filling fraction
    filling = 1 / m

    # Calculate the characteristic length scale
    rm = sqrt(2 * num_parts / filling)

    # Generate random real and imaginary parts in one step
    real_parts = rand(Float64, num_parts) .* rand(-1:2:1, num_parts) .* rm
    imag_parts = rand(Float64, num_parts) .* rand(-1:2:1, num_parts) .* rm

    # Combine real and imaginary parts to create complex numbers
    config = real_parts .- im .* imag_parts

    return config
end

function jastrow_squared(z,which)
	result = 1.0
	for i in 1:length(z)
		if z[i] != which
			result *= (which - z[i])^2
		end
	end
	return result
end

function firstderiv_j2(z,which)
	return conj(ForwardDiff.gradient(c -> real(jastrow_squared(z,c)),which)[1])
end

function nth_deriv_j2(z, which, n)
    if n == 0
        return jastrow_squared(z, which)  # 0th derivative is the function itself
    else
    	#println("Working on $(n-1)-th order")
        return conj(ForwardDiff.derivative(c -> real(nth_deriv_j2(z, c, n - 1)), which)[1])
    end
end

function test_jastrow_firstderiv(z)
	autodiff_result = firstderiv_j2(z,z[1])
	exact_result = 2*sum([(z[1] - z[k])*prod([(z[1] - z[j])^2 for j in deleteat!([l for l in 2:length(z)],k-1)]) for k in 2:length(z)])
	if isapprox(abs(autodiff_result-exact_result)/abs(exact_result),0.001;atol=10^-3)
		return true
	else
		return false,autodiff_result,exact_result
	end
end

function reverse_flux_wavefunction(z,m=3)
	num_parts = length(z)
	full_mat = im.*zeros(num_parts,num_parts)
	for i in 1:num_parts
		#println("Making Row $i")
		full_mat[i,:] = [log(Complex(nth_deriv_j2(z,z[j],i-1))) + log(2)*(i-1) for j in 1:num_parts]
	end
	
	result = get_log_det(full_mat)[1]
	
	result += -0.25*sum(abs2.(z))
	
	return result,[0,0]
end

function get_log_add(a,b)
	if real(a) > real(b)
		ordered::Vector{typeof(a)} = [b,a]
	else
		ordered = [a,b]
	end
	result::ComplexF64 = ordered[2] + log(Complex(1 + exp(ordered[1] - ordered[2])))
	return Complex(result)
end

function get_log_subtract(a,b)
	ordered = [a,b]
	result::ComplexF64 = ordered[1] + log(Complex(1 - exp(ordered[2] - ordered[1])))
	return Complex(result)
end

function get_log_det(matrix::Matrix{ComplexF64},reg_input=false)
	num_parts::Int64 = size(matrix)[1]
	maxes::Vector{ComplexF64} = [0.0+0.0*im for i in 1:num_parts]
	rejected_indices::Vector{Int64} = []
	
	changed::Matrix{ComplexF64} = matrix + fill(0.0,(num_parts,num_parts))
	for i in 1:num_parts
		row::Vector{Float64} = real(changed[i,:])
		validation::Bool = true
		j::Int64 = 0
		index::Int64 = 0
		while validation
			val::ComplexF64 = sort(row)[end - j]
			guess_index::Int64 = findfirst(real.(changed[i,:]) .== val )
			if any(rejected_indices .== index)
				j += 1
			else
				index = guess_index
				validation = false
			end
		end
		maxes[i] = changed[i,index]
		changed[i,:] .-= maxes[i]
		append!(rejected_indices,index)
	end
	reduced_logdet::ComplexF64 = sum(maxes) + log(Complex(det(exp.(changed))))

	return reduced_logdet,maxes,changed
end


function dist_btw(z,i,j,pow)
	return (z[i] - z[j])^pow
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

function test_3parts_jainkamila(z,m=3)
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
	matver[2,:] = (2*p) .* (fulljs[1,:] .* fulljs[2,:])
	matver[3,:] = (4*p) .* (((fulljs[2,:] .^ 2) ./ 1) .+ (fulljs[1,:] .* fulljs[3,:]))

	return log(Complex(det(matver))) - 0.25*sum(abs2.(z))
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

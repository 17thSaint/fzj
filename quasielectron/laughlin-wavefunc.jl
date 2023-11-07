using TaylorDiff,LinearAlgebra

include("general-funcs.jl")

function particles_within_radius(particle_positions, center_point, radius)
    """
    Find the positions of particles within a given radius from a specified center particle (efficient version).

    Parameters:
        particle_positions: Array{Complex{Float64}} - An array of complex positions of particles.
        center_index: Int - The index of the center particle.
        radius: Float64 - The radius within which to find particles.

    Returns:
        Array{Complex{Float64}} - An array of complex positions of particles within the specified radius.
    """
    #center = particle_positions[center_index]
    return [particle for particle in particle_positions if abs(particle - center_point) <= radius]
end

function laughlin_wavefunction_girvinjach(z, m=3; kwargs...)
    """
    Compute the logarithm of the Laughlin wavefunction for a given set of particle positions.
    
    Parameters:
        z: Array{Complex{Float64}} - An array of complex positions of particles.
        m: Int - An integer parameter that determines the Laughlin state (e.g., 1 for Laughlin ν=1/3 state).

    Returns:
        Float64 - The natural logarithm of the Laughlin wavefunction at the given particle positions.
    """  
    qe_loc = get(kwargs, :qe_loc, nothing)
    if_qe = !isnothing(qe_loc)
      
    N = length(z)  # Number of particles
    log_norm_factor = 0.0
    log_qh_part = 0.0
    for i in 1:N
        for j in 1:i-1
            log_norm_factor += m*log(Complex(z[i] - z[j]))
        end
        if if_qe
	        log_qh_part += log(Complex(z[i] - qe_loc))
	end
    end
    
    log_exponent = -0.25 * sum(abs2.(z))
    
    return log_norm_factor + log_exponent + log_qh_part,[]
end

function firstderiv_znJ1(z,which,slater_order; kwargs...)
	if_qe = get(kwargs, :if_qe, false)
	qe_loc = get(kwargs, :qe_loc, 0.0+im*0.0)
	
	if if_qe
		result = 2*derivative(c -> (c^slater_order) * jastrow(z,c),which,1)
		if qe_loc != 0.0
			result -= qe_loc * (which^slater_order) * jastrow(z,which)
		end
		return result
	else
		return (which^slater_order) * jastrow(z,which)
	end
end

# the qe here is of the form of a quasihole
function laughlin_wavefunction(z, m=3; kwargs...)
	num_parts = length(z)
	qe_cutoff = get(kwargs, :qe_cutoff, 0)
	qe_loc = get(kwargs, :qe_loc, 0.0+im*0.0)
	full_mat = im.*zeros(num_parts,num_parts)
	if_qe = qe_cutoff > 0
	for i in 1:num_parts
		#println("Making Row $i")
		if i > qe_cutoff
			if_qe = false
		end
		full_mat[i,:] = [jastrow(z,z[j]; if_log=true) + (i-1)*log(Complex(z[j])) for j in 1:num_parts]
		if if_qe
			full_mat[i,:] .+= [log(Complex(z[j] - qe_loc)) for j in 1:num_parts]
		end
	end
	
	result = get_log_det(full_mat)[1]
	
	result += -0.25*sum(abs2.(z))
	
	return result,[0,0]
end












"fin"

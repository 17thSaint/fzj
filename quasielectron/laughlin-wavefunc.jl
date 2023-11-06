
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

function log_laughlin_wavefunction(z, m; kwargs...)
    """
    Compute the logarithm of the Laughlin wavefunction for a given set of particle positions.
    
    Parameters:
        z: Array{Complex{Float64}} - An array of complex positions of particles.
        m: Int - An integer parameter that determines the Laughlin state (e.g., 1 for Laughlin ν=1/3 state).

    Returns:
        Float64 - The natural logarithm of the Laughlin wavefunction at the given particle positions.
    """
    N = length(z)  # Number of particles
    log_norm_factor = 0.0
    for i in 1:N
        for j in 1:i-1
            log_norm_factor += m*log(Complex(z[i] - z[j]))
        end
    end
    
    log_exponent = -0.25 * sum(abs2.(z))
    
    return log_norm_factor + log_exponent,[]
end












"fin"

#####################################################
#=

This file is for calculating observables for effective 1D MPS

Depends on:
    

=#
######################################################

function make_density_correlations(wavefunc::MPS; kwargs...)
    Lphys::Int64,Lsynth::Int64 = get_mps_dims(wavefunc)

    density_correlations::Array{ComplexF64,4} = zeros(Lphys,Lphys,Lsynth,Lsynth)

    for s in 1:Lsynth
        for ss in 1:Lsynth
            corr_mat::Matrix{ComplexF64} = correlation_matrix(wavefunc,"Ns$(s)","Ns$(ss)")
            density_correlations[:,:,s,ss] = corr_mat
        end
    end

    return density_correlations
    
end



denscorrs = make_density_correlations(psis[1])
rez = range_cdwsf_radii(100,denscorrs)










































"fin"
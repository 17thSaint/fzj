using Pkg
Pkg.activate(".")
include("fqh_effective.jl")
include("time_evolution.jl")
include("../other-funcs/data-storage-funcs.jl")
using Statistics,Observers,LsqFit



function hamiltonian_remapped(L,tp=1.0; kwargs...)
        if_periodic_phys = get(kwargs, :if_periodic_phys, true)
        tilt_strength = get(kwargs, :tilt_strength, 0.0)
        centralflux_strength = get(kwargs, :centralflux_strength, 0.0)
        if_remapping = get(kwargs, :if_remapping, true)
        
        remap = if_remapping ? remapping_nnn(L) : nothing

        ampo = OpSum()
        for j in 1:L
            # physical dimension hopping
            next_site = j+1
            if j == L
                if if_periodic_phys
                    next_site = 1
                else
                    continue
                end
            end
            if if_remapping
                next_site = remap[next_site]
                starting_site = remap[j]
            else
                starting_site = j
            end

            
            coeff = -tp * exp(im*2*pi*centralflux_strength/L)
            ampo += (coeff, "Adag", starting_site, "A", next_site)
            ampo += (conj(coeff), "A", starting_site, "Adag", next_site)
        end
        
        if_tilt = tilt_strength != 0.0
        if if_tilt
            for j in 1:L
                phys_site = if_remapping ? remap[j] : j
                ampo += (-tilt_strength*j, "N", phys_site)
            end
        end
        
        return ampo
end

nrgvar_tol = 1E-6
Ls = range(6,step=4,length=10)
linear_bonddims = zeros(length(Ls))
remapped_bonddims = zeros(length(Ls))

for (idx,L) in enumerate(Ls)
    remap = remapping_nnn(L)
    for if_remapping in [true,false]
        ham = hamiltonian_remapped(L; if_remapping=if_remapping)
        obs = NRGVarObserver(nrgvar_tol,ham)

        part_count = Int(floor(L/2))
        mdim = 250

        dmrg_params = (if_densmat=false,nsweeps=50,if_parton=false,particle_type="Boson",conserve_qns=true,ham=ham,mdim=mdim,observer=obs)

        psi_gs = execute_mps(nothing,nothing,nothing,L,1,part_count; dmrg_params...)

        if_remapping ? remapped_bonddims[idx] = maxlinkdim(psi_gs) : linear_bonddims[idx] = maxlinkdim(psi_gs)
    end
end

plot(Ls,linear_bonddims,"-p",label="Linear")
plot(Ls,remapped_bonddims,"-p",label="Remapped")
legend()
xlabel("System size")
ylabel("Max bond dimension")

#=
corrmat = correlation_matrix(psi_gs, "Adag", "A")
dists = collect(0:L-1)
corrs = [mean(diag(corrmat,i)) for i in 0:L-1] ./ corrmat[1,1]

remapped_corrmat = zeros(L,L)
for i in 1:L
    for j in 1:L
        #println("Site $i,$j goes to $(remap[i]),$(remap[j])")
        remapped_corrmat[i,j] = corrmat[remap[i],remap[j]]
    end
end
remapped_corrs = [mean(diag(remapped_corrmat,i)) for i in 0:L-1] ./ remapped_corrmat[1,1]

plot(dists,corrs,"-p",label="Original")
plot(dists,remapped_corrs,"-p",label="Remapped")
legend()
xlabel("Distance")
ylabel("Correlation")
=#





































"fin"
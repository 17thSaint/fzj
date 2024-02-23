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

        for j in 1:L
            if if_remapping
                phys_site = remap[j]
            else
                phys_site = j
            end
            ampo += (1000000.0, "Adag * A * Adag * A", phys_site)
            ampo -= (1000000.0, "Adag * A", phys_site)
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

open_cores = "all"#get(params_dict, "open_cores", "all")
if typeof(open_cores) != String
	BLAS.set_num_threads(open_cores)	
end

nrg_tol = 1E-10
Ls = range(10,step=2,length=10)
linear_bonddims = zeros(length(Ls))
remapped_bonddims = zeros(length(Ls))
mdim = [100,250,500]
L = 16
dists = collect(0:L-1)

println("Working on L = $L")
remap = remapping_nnn(L)
part_count = Int(floor(L/2))

states = make_states(L,part_count,1)
sidx = siteinds("Boson", L; conserve_qns = true, dim=2)
psi0 = randomMPS(sidx, states)

rerun = true

if rerun
ham_remapped = hamiltonian_remapped(L; if_remapping=true)
obs_remapped = NRGObserver(nrg_tol,ham_remapped)
dmrg_params_remapped = (outputlevel=1,psi_guess=psi0,if_densmat=false,nsweeps=50,if_parton=false,particle_type="Boson",conserve_qns=true,ham=ham_remapped,mdim=mdim,observer=obs_remapped)
psi_remapped = execute_mps(nothing,nothing,nothing,L,1,part_count; dmrg_params_remapped...)
end
corrmat_remap = correlation_matrix(psi_remapped, "Adag", "A")
remapped_corrmat = zeros(L,L)
for i in 1:L
    for j in 1:L
        remapped_corrmat[i,j] = corrmat_remap[remap[i],remap[j]]
    end
end
remapped_corrs = [mean(diag(remapped_corrmat,i)) for i in 0:L-1] #./ remapped_corrmat[1,1]
plot(dists,remapped_corrs,"-p",label="Remap")

if rerun
ham_linear = hamiltonian_remapped(L; if_remapping=false)
obs_linear = NRGObserver(nrg_tol,ham_linear)
dmrg_params_linear = (outputlevel=1,psi_guess=psi0,if_densmat=false,nsweeps=50,if_parton=false,particle_type="Boson",conserve_qns=true,ham=ham_linear,mdim=mdim,observer=obs_linear)
psi_linear = execute_mps(nothing,nothing,nothing,L,1,part_count; dmrg_params_linear...)
end
corrmat = correlation_matrix(psi_linear, "Adag", "A")
corrs = [mean(diag(corrmat,i)) for i in 0:L-1] #./ corrmat[1,1]
plot(dists,corrs,"-p",label="Linear")

xlabel("Distance")
ylabel("Correlations")
yscale("log")
legend()

#end
#=
linear_bonddims = [8.0,32.0,111.0,259.0,416.0,592.0,772.0,958.0]
remapped_bonddims = [8.0,32.0,114.0,259.0,415.0,595.0,782.0,967.0]
=#







































"fin"
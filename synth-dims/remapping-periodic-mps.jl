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

function get_1D_gs(L,if_remapping,nrg_tol,mdim; kwargs...)
    part_count = Int(floor(L/2))
    
    states = make_states(L,part_count,1)
    sidx = siteinds("Boson", L; conserve_qns = true, dim=2)
    psi0 = randomMPS(sidx, states)

    ham = hamiltonian_remapped(L; if_remapping=if_remapping)
    obs = NRGObserver(nrg_tol,ham)
    dmrg_params = (outputlevel=1,psi_guess=psi0,if_densmat=false,nsweeps=50,if_parton=false,particle_type="Boson",conserve_qns=true,ham=ham,mdim=mdim,observer=obs)
    psi = execute_mps(nothing,nothing,nothing,L,1,part_count; dmrg_params...)
    
    return psi
end

function get_corrs(psi::MPS,if_remapping::Bool; kwargs...)
    if_plot = get(kwargs, :if_plot, true)

    L = length(psi)
    remap = remapping_nnn(L)
    dists = collect(0:L-1)
    corrmat = correlation_matrix(psi, "Adag", "A")

    if if_remapping
        remapped_corrmat = zeros(L,L)
        for i in 1:L
            for j in 1:L
                remapped_corrmat[i,j] = corrmat[remap[i],remap[j]]
            end
        end
        remapped_corrs = [mean(diag(remapped_corrmat,i)) for i in 0:L-1] ./ remapped_corrmat[1,1]
        if if_plot
            plot(dists,remapped_corrs,"-p",label="Remap")
            xlabel("Distance")
            ylabel("Correlations")
            yscale("log")
            legend()
        end
    else
        corrs = [mean(diag(corrmat,i)) for i in 0:L-1] ./ corrmat[1,1]
        if if_plot
            plot(dists,corrs,"-p",label="Linear")
            xlabel("Distance")
            ylabel("Correlations")
            yscale("log")
            legend()
        end
    end

    if if_remapping
        return dists,remapped_corrs
    else
        return dists,corrs
    end
end


    

open_cores = "all"#get(params_dict, "open_cores", "all")
if typeof(open_cores) != String
	BLAS.set_num_threads(open_cores)	
end

nrg_tol = 1E-10
Ls = range(10,step=2,length=20)
linear_bonddims = zeros(length(Ls))
remapped_bonddims = zeros(length(Ls))
mdim = [100,250,500]

corrdiffs_linear = zeros(length(Ls))
corrdiffs_remap = zeros(length(Ls))

for (idx,L) in enumerate(Ls)

    println("Working on L = $L")
    psi_linear = get_1D_gs(L,false,nrg_tol,mdim)
    mld = maxlinkdim(psi_linear)
    while mld < 2
        psi_linear = get_1D_gs(L,false,nrg_tol,mdim)
        mld = maxlinkdim(psi_linear)
    end
    #dists_linear,corrs_linear = get_corrs(psi_linear,false;if_plot=false)

    psi_remap = get_1D_gs(L,true,nrg_tol,mdim)
    mldr = maxlinkdim(psi_remap)
    while mldr < 2
        psi_remap = get_1D_gs(L,true,nrg_tol,mdim)
        mldr = maxlinkdim(psi_remap)
    end
    #dists_remap,corrs_remap = get_corrs(psi_remap,true;if_plot=false)
    
    #middle = Int(ceil(3*L/4))
    corrdiffs_linear[idx] = abs(correlation_matrix(psi_linear,"Adag","A";sites=(1,L))[1,2] - correlation_matrix(psi_linear,"Adag","A";sites=(1,2))[1,2])
    corrdiffs_remap[idx] = abs(correlation_matrix(psi_remap,"Adag","A";sites=(1,L))[1,2] - correlation_matrix(psi_remap,"Adag","A";sites=(1,2))[1,2])

    println("Linear is ",corrdiffs_linear[idx]," and Remap is ",corrdiffs_remap[idx])

    if idx > 1
        plot(Ls[idx-1:idx],corrdiffs_linear[idx-1:idx],"-p",c="b")
        plot(Ls[idx-1:idx],corrdiffs_remap[idx-1:idx],"-p",c="r")
        yscale("log")
    end
end

fig2 = figure()
plot(Ls,corrdiffs_linear,"-p",label="Linear")
plot(Ls,corrdiffs_remap,"-p",label="Remap")
yscale("log")
xlabel("System Size")
ylabel("AVG(Adag_i A_i+L-1) - AVG(Adag_i A_i+)")
title("Long Range Correlation Difference Linear vs NNN Remap")








































"fin"
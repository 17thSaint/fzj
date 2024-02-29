using Pkg
Pkg.activate(".")
include("fqh_effective.jl")
include("time_evolution.jl")
include("../other-funcs/data-storage-funcs.jl")
using Statistics,Observers,LsqFit,SpecialFunctions

function entropy_mps(offcenter_wavefunc::MPS,center_site::Int; kwargs...)
	if_spectrum = get(kwargs, :if_spectrum, false)
	wavefunc = orthogonalize(offcenter_wavefunc,center_site)
	u,s,v = svd(wavefunc[center_site], (linkind(wavefunc,center_site)))#, siteind(wavefunc,center_site)))
	entropy = 0.0
	for i in 1:dim(s,1)
		p = s[i,i]^2
		entropy -= p * log(p)
	end
	if if_spectrum
		return entropy,diag(s)
	else
		return entropy
	end
end

function hamiltonian_remapped(L,tp=1.0; kwargs...)
        if_periodic_phys = get(kwargs, :if_periodic_phys, true)
        tilt_strength = get(kwargs, :tilt_strength, 0.0)
        centralflux_strength = get(kwargs, :centralflux_strength, 0.0)
        if_remapping = get(kwargs, :if_remapping, true)
        chemical_strength = get(kwargs, :chemical_strength, 0.0)
        onsite_strength = get(kwargs, :onsite_strength, 1.0)
        
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

            
            coeff = -tp
            ampo += (coeff, "Adag", starting_site, "A", next_site)
            ampo += (conj(coeff), "A", starting_site, "Adag", next_site)
        end

        # onsite interaction
        for j in 1:L
            if if_remapping
                phys_site = remap[j]
            else
                phys_site = j
            end
            ampo += (onsite_strength, "N * N", phys_site)
            ampo -= (onsite_strength, "N", phys_site)
        end
    
        if chemical_strength != 0.0
            for j in 1:L
                phys_site = if_remapping ? remap[j] : j
                ampo += (chemical_strength, "N", phys_site)
            end
        end
        
        return ampo
end

function get_1D_gs(L,if_remapping,nrg_tol,mdim; kwargs...)
    part_count = Int(floor(L/2))
    maxocc = get(kwargs, :maxocc, 1)
    psi_ortho = get(kwargs, :psi_ortho, nothing)
    cutoff = get(kwargs, :cutoff, 1E-10)
    hopping_strength = get(kwargs, :hopping_strength, 1.0)
    
    states = make_states(L,part_count,1)
    sidx = isnothing(psi_ortho) ? siteinds("Boson", L; conserve_qns = true, dim=maxocc+1) : siteinds(psi_ortho)
    psi0 = randomMPS(sidx, states)

    ham = hamiltonian_remapped(L,hopping_strength; if_remapping=if_remapping,kwargs...)
    obs = NRGObserver(nrg_tol,ham)
    dmrg_params = (cutoff=cutoff,psi_ortho=psi_ortho,outputlevel=1,psi_guess=psi0,if_densmat=false,nsweeps=100,if_parton=false,particle_type="Boson",conserve_qns=true,ham=ham,mdim=mdim,observer=obs)
    psi, nrg = execute_mps(nothing,nothing,nothing,L,1,part_count; dmrg_params...,if_nrg=true)
    
    return psi, nrg
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

cutoff = 1E-8
mu = 0.0
ts = 1.0
us = 0.1
#locL = 10
#ustrens = range(0.01,stop=10.0,length=10)
Ls = collect(10:4:50)
#intparts = zeros(length(Ls))
#zeromomoccs = zeros(length(Ls))
#for (idx,us) in enumerate(ustrens)
for (idx,locL) in enumerate(Ls)
    dists = collect(0:locL-1)
    psi_loc, gs_nrg_loc = get_1D_gs(locL,false,1E-4,100; onsite_strength = us,hopping_strength=ts,chemical_strength=mu,if_periodic_phys=false,cutoff=cutoff,maxocc=5)
    corrmat = correlation_matrix(psi_loc, "Adag", "A")
    siteocc = expect(psi_loc,"N")
    normalization_mat = sqrt.(siteocc * transpose(siteocc))
    corrmat ./= normalization_mat
    corrs = [mean(diag(corrmat,i)) for i in 0:locL-1]
    #fig = figure()
    #plot(corrs)
    #yscale("log")
    decayfit(x,p) = (p[1].* exp.(-x ./ p[2]) ./ (x .^ 0.5)) .+ p[3]
    starting_point = 5
    decayfithere = curve_fit(decayfit,dists[starting_point:end],corrs[starting_point:end],[1.0,1.0,0.0])
    corrlength = decayfithere.param[2]

    integralpart = sqrt(corrlength) * gamma(1/2) #/ (locL^1.5)
    sumpart = sum(corrmat) / locL^2
    scatter(locL,integralpart,c="b")
    #scatter(us,sumpart,c="r")
end



#=
fig = figure()
nrg_diffs = []
Ls = collect(60:4:100)
for L in Ls
    psi, gs_nrg = get_1D_gs(L,false,1E-8,200; if_periodic_phys=false,cutoff=cutoff,hopping_strength=ts,chemical_strength=mu,maxocc=1)
    #=excited_psi, first_nrg = get_1D_gs(L,false,1E-8,200; if_periodic_phys=false,psi_ortho=psi,cutoff=cutoff)
    nrg_diff = first_nrg - gs_nrg
    append!(nrg_diffs,[nrg_diff])
    scatter(L,nrg_diff,c="b")
    xlabel("1 / System Size")
    ylabel("Energy Difference")=#
    append!(nrg_diffs,[entropy_mps(psi,Int(ceil(L/2)))])
    scatter(L,nrg_diffs[end],c="b")
end

xs = Ls
svnfunc(x,p) = p[1] .* log.((2/pi) .* x .* sin(pi/2)) .+ p[2]
svnfit = curve_fit(svnfunc,xs,nrg_diffs,[1.0,1.0])
plot(xs,svnfunc(xs,svnfit.param),label="Fit")
title("Central Charge = $(svnfit.param[1]*6)")
=#

#=
xs = collect(10:4:50)
plot(xs,nrg_diffs,"-p")
ratfit(x,p) = p[1] ./ x .+ p[2]
fit = curve_fit(ratfit,xs,nrg_diffs,[1.0,1.0])
plot(xs,ratfit(xs,fit.param),label="Fit")
title("Limit is $(fit.param[2])")
=#
#=

open_cores = "all"#get(params_dict, "open_cores", "all")
if typeof(open_cores) != String
	BLAS.set_num_threads(open_cores)	
end

nrg_tol = 1E-10
Ls = range(10,step=2,length=40)
linear_bonddims = zeros(length(Ls))
remapped_bonddims = zeros(length(Ls))
mdim = [100,250,500]
#
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
    remap = remapping_nnn(L)
    corrdiffs_linear[idx] = abs(correlation_matrix(psi_linear,"Adag","A";sites=(1,L))[1,2] - correlation_matrix(psi_linear,"Adag","A";sites=(1,2))[1,2]) / correlation_matrix(psi_linear,"Adag","A";sites=(1,2))[1,2]
    corrdiffs_remap[idx] = abs(correlation_matrix(psi_remap,"Adag","A";sites=(1,remap[L]))[1,2] - correlation_matrix(psi_remap,"Adag","A";sites=(1,remap[2]))[1,2]) / correlation_matrix(psi_remap,"Adag","A";sites=(1,remap[2]))[1,2]

    println("Linear is ",corrdiffs_linear[idx]," and Remap is ",corrdiffs_remap[idx])

    if idx > 1
        plot(Ls[idx-1:idx],corrdiffs_linear[idx-1:idx],"-p",c="b")
        plot(Ls[idx-1:idx],corrdiffs_remap[idx-1:idx],"-p",c="r")
        yscale("log")
    end
end
#

exp_fit(x,p) = p[1].* exp.(x ./ p[2]) .+ p[3]

fit_linear = curve_fit(exp_fit,Ls[4:end],corrdiffs_linear[4:end],[1.0,1.0,0.0])
fit_remap = curve_fit(exp_fit,Ls[3:end],corrdiffs_remap[3:end],[1.0,1.0,0.0])

fig2 = figure()
scatter(Ls,corrdiffs_linear,label="Linear")
plot(Ls[4:end],exp_fit(Ls[4:end],fit_linear.param))
scatter(Ls,corrdiffs_remap,label="Remap")
plot(Ls[3:end],exp_fit(Ls[3:end],fit_remap.param))
yscale("log")
xlabel("System Size")
ylabel("(Adag_1 A_2 - Adag_1 A_L) / Adag_1 A_2")
title("Percent Long Range Correlation Diff Linear vs NNN Remap")
legend()

=#






































"fin"
include("fqh_effective.jl")
include("time_evolution.jl")
include("../other-funcs/data-storage-funcs.jl")
using Statistics,PyPlot,Observers,ITensorTDVP,LsqFit

lin_model(x,p) = p[1].* x .+ p[2]

function firstderiv_ham(s,j,L,nflavors,chi,tp=1.0,ts=1.0; kwargs...)
    if_periodic_phys = get(kwargs, :if_periodic_phys, false)
    centralflux_strength = get(kwargs, :centralflux_strength, 0.0)
    if_s0 = get(kwargs, :if_s0, true)
    
    s0 = 0.0
    if if_s0
        s0 = nflavors/2
    end

    ampo = OpSum()
    next_site = j+1
    if j == L
        if if_periodic_phys
            next_site = 1
        end
    end
    ampo += (-tp * (im*2*pi/L) * exp(im*chi*(s-s0)/(nflavors)) * exp(im*2*pi*centralflux_strength/L), "Cr$s", j, "Anh$s", next_site)
    ampo += (-tp * (-im*2*pi/L) * exp(-im*chi*(s-s0)/(nflavors)) * exp(-im*2*pi*centralflux_strength/L), "Anh$s", j, "Cr$s", next_site)
    
    return ampo
end

function secderiv_ham(s,j,L,nflavors,chi,tp=1.0,ts=1.0; kwargs...)
    if_periodic_phys = get(kwargs, :if_periodic_phys, false)
    centralflux_strength = get(kwargs, :centralflux_strength, 0.0)
    if_s0 = get(kwargs, :if_s0, true)
    
    s0 = 0.0
    if if_s0
        s0 = nflavors/2
    end

    ampo = OpSum()
    next_site = j+1
    if j == L
        if if_periodic_phys
            next_site = 1
        end
    end
    ampo += (tp * (2*pi/L)^2 * exp(im*chi*(s-s0)/(nflavors)) * exp(im*2*pi*centralflux_strength/L), "Cr$s", j, "Anh$s", next_site)
    ampo += (tp * (2*pi/L)^2 * exp(-im*chi*(s-s0)/(nflavors)) * exp(-im*2*pi*centralflux_strength/L), "Anh$s", j, "Cr$s", next_site)
    
    return ampo
end

function calc_deriv(order,psi,s,j,nflavors,chi,ham_params)
    L = length(psi)
    if order == 1
        mpo_ver = MPO(firstderiv_ham(s,j,L,nflavors,chi; ham_params...),siteinds(psi))
    elseif order == 2
        mpo_ver = MPO(secderiv_ham(s,j,L,nflavors,chi; ham_params...),siteinds(psi))
    end
    return inner(psi',mpo_ver,psi) / inner(psi,psi)
end

function allsite_deriv(order,psi,nflavors,chi,ham_params)
    L = length(psi)
    currents = zeros(nflavors,L) .* im
    for j in 1:L
        for s in 1:nflavors
            if order == 1
                mpo_ver = MPO(firstderiv_ham(s,j,L,nflavors,chi; ham_params...),siteinds(psi))
            elseif order == 2
                mpo_ver = MPO(secderiv_ham(s,j,L,nflavors,chi; ham_params...),siteinds(psi))
            end
            currents[s,j] = inner(psi',mpo_ver,psi) / inner(psi,psi)
        end
    end
    return currents
end

function hamiltonian_universal(L,nflavors,chi,tp=1.0,ts=1.0; kwargs...)
        if_periodic_phys = get(kwargs, :if_periodic_phys, false)
        if_periodic_synth = get(kwargs, :if_periodic_synth, false)
        tilt_strength = get(kwargs, :tilt_strength, 0.0)
        centralflux_strength = get(kwargs, :centralflux_strength, 0.0)
        if_s0 = get(kwargs, :if_s0, true)
        
        s0 = 0.0
        if if_s0
            s0 = nflavors/2
        end

        ampo = OpSum()
        for j in 1:L
            for s in 1:nflavors
                # physical dimension hopping
                next_site = j+1
                if j == L
                    if if_periodic_phys
                        next_site = 1
                    else
                        continue
                    end
                end
                ampo += (-tp * exp(im*chi*(s-s0)/(nflavors)) * exp(im*2*pi*centralflux_strength/L), "Cr$s", j, "Anh$s", next_site)
                ampo += (-tp * exp(-im*chi*(s-s0)/(nflavors)) * exp(-im*2*pi*centralflux_strength/L), "Anh$s", j, "Cr$s", next_site)
            end
        end
        
        for j in 1:L
            for s in 1:nflavors
                # synthetic dimension hopping
                next_site = s+1
                if s == nflavors
                    if if_periodic_synth
                        next_site = 1
                    else
                        continue
                    end
                end
                ampo += (-ts * 1.0, "Cr$(next_site) * Anh$(s)", j)
                ampo += (-ts * 1.0, "Cr$(s) * Anh$(next_site)", j)
                #ampo += (-t2 * exp(-im*phi*j), "Anh$(next_site) * Cr$(s)", j)
            end
        end
        
        if_tilt = tilt_strength != 0.0
        if if_tilt
            for j in 1:L
                for s in 1:nflavors
                    ampo += (-tilt_strength*j, "Ns$(s)", j)
                end
            end
        end
        
        return ampo
end

#
if_save_data = true
dataloc = "/home/patrick/fzj/main-git/cluster-data/synth-dims/"

nrgvar_tol = 1E-10
mdim = 100

#geo_params = [()] # (L,nf,nb)
Ls = [4,5,6,7,8,9,10,11,12,13,14,15]

nu = 1.0
#
states = []
for L in Ls
#L = 8
nflavors = Int(ceil(0.75*L))
part_count = nflavors
#chi = 0.0#part_count / (nu*L*nflavors)
tilt = 0.0

if_per_phys = true
if_per_virt = false

centralflux_strength = 0.0

naming_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",part_count),("mdim",mdim),("centralflux_strength",centralflux_strength)])
metadata = merge(naming_dict,Dict([("if_periodic_phys",if_per_phys),("if_periodic_virt",if_per_virt),("tilt_strength",tilt),("location",dataloc),("if_save_data",if_save_data),("nrgvar_tol",nrgvar_tol),("mdim",mdim)]))



#current_strength = 0.00

#
if true
counting = 20
scaling = 64
strens = [part_count / (nu*L*nflavors)]#range(0.1,stop=0.15,length=counting)#0.5 .+ [sort([-i/scaling for i in 1:counting]); [0.0]; [i/scaling for i in 1:counting]]
nrgs = zeros(length(strens)) .* im
currents = zeros(nflavors,length(strens)) .* im
drudes = zeros(nflavors,length(strens)) .* im

println("Chi = ",part_count / (nu*L*nflavors))
#
for (i,chi) in enumerate(strens)
    metadata["chi"] = chi
    naming_dict["chi"] = round(chi,digits=5)
    #=
    metadata["centralflux_strength"] = centralflux_strength
    if centralflux_strength < 0.0
        naming_dict["centralflux_strength"] = "n" * string(-round(centralflux_strength,digits=5))
    else
        naming_dict["centralflux_strength"] = round(centralflux_strength,digits=5)
    end
    =#
    filename = make_parameters_filename(naming_dict)
    metadata["name"] = filename
    display(filename)

    ham_params = (if_periodic_phys=if_per_phys,if_periodic_synth=if_per_virt,centralflux_strength=centralflux_strength,tilt_strength=0.0)
    ham_start = hamiltonian_universal(L,nflavors,chi; ham_params...)
    obs = NRGVarObserver(nrgvar_tol,ham_start)

    if_exists,found_data = check_data_exists("mps-"*filename*".jld2","mps";location=dataloc)

    if if_exists
        psi_gs = found_data
    else
        dmrg_params = (ham=ham_start,mdim=mdim,if_save_data=if_save_data,metadata=metadata,name=filename,location=dataloc,observer=obs)
        psi_gs = execute_mps(nothing,nothing,chi,L,nflavors,part_count; dmrg_params...)
        println("Energy Variance = ",energy_variance(psi_gs,ham_start)," at Chi = ",chi)
    end

    #append!(states,[psi_gs])
    nrgs[i] = calculate_energy(psi_gs,ham_start)

    centersite = Int(ceil(L/2))
    u,s,v = svd(psi_gs[centersite],linkind(psi_gs,centersite))
    ents = sum([s[n,n]^2 * 2 * log(s[n,n]) for n in 1:size(s)[1]])
    scatter([L],[ents],c="b")

    #currents[:,i] = [calc_deriv(1,psi_gs,s,Int(L/2),nflavors,chi,ham_params) for s in 1:nflavors]
    #drudes[:,i] = [calc_deriv(2,psi_gs,s,Int(L/2),nflavors,chi,ham_params) for s in 1:nflavors]
    #=fig = figure()
    imshow(real.(drude))
    title("Drudes at Current Strength = $(centralflux_strength)")
    xlabel("Physical Site")
    ylabel("Synthetic Site")
    colorbar()
    =#
end
#
#=
fig = figure()
plot(strens,real.(nrgs),"-p",label="Energy")
xlabel("Phi")
legend()
fig2 = figure()
for i in 1:nflavors
    plot(strens,real.(currents[i,:]),"-p",label="Current $i")
end
xlabel("Phi")
legend()
fig3 = figure()
for i in 1:nflavors
    plot(strens,real.(drudes[i,:]),"-p",label="Drude $i")
end
xlabel("Phi")
legend()
=#
end
#
end
#

#
if false
    time_end = 10.0
    time_change = 0.1
    mdim_time = 50

    chi = strens[1]
    println("Chi = ",chi)
    #
    tilts = [0.00]
    time_states = [[] for i in 1:length(tilts)]
    time_occs = [[] for i in 1:length(tilts)]
    alltimes = [[] for i in 1:length(tilts)]
    time_currents = [zeros(nflavors,Int(time_end/time_change)+1) .* im for i in 1:length(tilts)]
    for (i,tilt) in enumerate(tilts)
        ham_params_evolve = (if_periodic_phys=if_per_phys,if_periodic_synth=if_per_virt,centralflux_strength=centralflux_strength,tilt_strength=tilt)
        ham_evolve = hamiltonian_universal(L,nflavors,chi; ham_params_evolve...)

        psi_gs = states[1]
        #time_currents[i][:,1] = [calc_deriv(1,psi_gs,s,Int(L/2),nflavors,chi,ham_params_evolve) for s in 1:nflavors]
        #rez0, otherham0 = evolve_in_time(psi_gs,time_end,time_change,ham_start; mdim=mdim_time,obs_measures=Dict("occs" => current_occ, "states" => return_state, "nrg_vars" => current_nrgvar, "nrgs" => current_nrg))
        rez, otherham = evolve_in_time(psi_gs,time_end,time_change,ham_evolve; mdim=mdim_time,obs_measures=Dict("occs" => current_occ, "states" => return_state))
        times = [[0.0]; rez["times"].results]
        time_states[i] = [[psi_gs]; rez["states"].results]
        time_occs[i] = [[get_occupancy(psi_gs;if_plot=false)]; rez["occs"].results]
        alltimes[i] = times
    end
    #

    times = alltimes[1]
    ham_params_evolve = (if_periodic_phys=if_per_phys,if_periodic_synth=if_per_virt,centralflux_strength=centralflux_strength,tilt_strength=0.0)
    for j in 1:length(times)-1
        #fig = figure()
        #imshow(real.(allsite_deriv(1,rez["states"].results[j],nflavors,chi,ham_params_evolve)))
        #colorbar()
        time_currents[1][:,j+1] = [calc_deriv(1,time_states[1][j],s,Int(L/2),nflavors,chi,ham_params_evolve) for s in 1:nflavors]
    end

    sumcurrents = [0.0 for i in 1:length(times)]
    for j in 1:length(times)
        sumcurrents[j] = real(sum(allsite_deriv(1,time_states[1][j],nflavors,chi,ham_params_evolve)))
    end
    fig = figure()
    plot(times,sumcurrents,"-p")
    xlabel("Time")
    ylabel("Total Current")

    #println(length(times),", ",length())
    fig = figure()
    for j in 1:nflavors
        plot(times,real.(time_currents[1][j,:]) .- real.(time_currents[1][j,1]),"-p",label="$j")
    end
    xlabel("Time")
    ylabel("Current")
    title("Tilt = $(round(tilt,digits=3)), Chi = $(round(chi,digits=3))")
    legend()

    denspols = [density_polarization(nothing,time_occs[1][j]) for j in 1:length(times)]
    fig2 = figure()
    plot(times,real.(denspols) .- real(denspols[1]),"-p")
    xlabel("Times")
    ylabel("Density Polarization")

    spacdenspols = [spacial_density_polarization(nothing,time_occs[1][j])[2] for j in 1:length(times)]
    fig2 = figure()
    for i in 1:nflavors
        yvals = real.([spacdenspols[j][i] for j in 1:length(times)])
        plot(times,yvals .- yvals[1],"-p",label="$i")
    end
    legend()
    xlabel("Times")
    ylabel("Spacial Density Polarization")

    #=fig2 = figure()
    plot(times,imag.(currents_null),"-p")
    title("Null")
    =#
end
#















"fin"
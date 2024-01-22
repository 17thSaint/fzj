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
                ampo += (-tp * exp(im*pi*chi*(s-s0)) * exp(im*2*pi*centralflux_strength/L), "Cr$s", j, "Anh$s", next_site)
                ampo += (-tp * exp(-im*pi*chi*(s-s0)) * exp(-im*2*pi*centralflux_strength/L), "Anh$s", j, "Cr$s", next_site)
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
if_save_data = false
dataloc = get_folder_location("cluster-data/synth-dims","geraghty")
if_densmat = true

nsweeps = 100
nrgvar_tol = 1E-7
mdim = 500
noise = [0.0]

#geo_params = [()] # (L,nf,nb)
#Ls = [4,5,6,7,8,9,10,11,12,13,14,15]

#nu = 1.0
#
states = []
#for L in Ls
L = 16
nflavors = 8#Int(ceil(0.75*L))
part_count = Int(L/2)
#chi = 0.0#part_count / (nu*L*nflavors)
tilt = 0.0

if_per_phys = true
if_per_virt = false

centralflux_strength = 0.0

naming_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",part_count),("centralflux_strength",centralflux_strength)])
metadata = merge(naming_dict,Dict([("if_periodic_phys",if_per_phys),("if_periodic_virt",if_per_virt),("tilt_strength",tilt),("location",dataloc),("if_save_data",if_save_data),("nrgvar_tol",nrgvar_tol),("mdim",mdim)]))

#
if true
counting = 20
#scaling = 64
strens = range(part_count/(0.3*L*nflavors),part_count/(0.37*L*nflavors),length=counting)#0.5 .+ [sort([-i/scaling for i in 1:counting]); [0.0]; [i/scaling for i in 1:counting]]
sf_orderparams = zeros(length(strens))
bonddims = zeros(length(strens))
distcorrs = zeros(length(strens))
ees = zeros(length(strens))
corrlengs = [zeros(length(strens)) for i in 1:nflavors]
#
#nrgs = zeros(length(strens)) .* im
#currents = zeros(nflavors,length(strens)) .* im
#drudes = zeros(nflavors,length(strens)) .* im
#
#=chi = 0.0
params_dict = Dict([("L",L),("nbosons",part_count),("nflavors",nflavors),("centralflux_strength",centralflux_strength)])
loc = "/home/patrick/fzj/main-git/cluster-data/synth-dims/"
all_files = find_data_file(params_dict,"mps",loc)
display(all_files)
#

#
strens = zeros(length(all_files))
#nrgs = zeros(length(all_files)) .* im
#currents = zeros(nflavors,length(all_files)) .* im
#drudes = zeros(nflavors,length(all_files)) .* im
sf_orderparams = zeros(length(all_files))
bonddims = zeros(length(all_files))
distcorrs = zeros(length(all_files))
ees = zeros(length(all_files))
corrlengs = [zeros(length(all_files)) for i in 1:nflavors]
=#

#println("Chi = ",part_count / (nu*L*nflavors))
for (idx,chi) in enumerate(strens)
#for (idx,f) in enumerate(all_files)
    if false
        found_data, found_metadata = read_data_jld2(f,loc)
        #centralflux_strength = found_metadata["centralflux_strength"]
        chi = found_metadata["chi"]
        strens[idx] = chi
        psi_gs = found_data["mps"]

        if part_count/(chi*nflavors*L) == 0.28 || isapprox(part_count/(chi*nflavors*L),0.38;atol=0.01)
            get_occupancy(psi_gs; plot_title="Chi = $(chi)")
        else
            continue
        end
        #=
        if "final_nrg_variance" in keys(found_metadata)
            fn_nrg_var = found_metadata["final_nrg_variance"]
            #println("Found in Metadata")
        else
            fn_nrg_var = energy_variance(psi_gs,found_metadata["ham"])
            #println("Didn't have it stored")
        end
        if fn_nrg_var > nrgvar_tol
            run_again(f; location=loc)
        end
        =#
        try
            densmat = found_data["densmat"]
        catch
            densmat = nothing
        end
    end
    #ham_params = (if_periodic_phys=if_per_phys,if_periodic_synth=if_per_virt,centralflux_strength=centralflux_strength,tilt_strength=0.0)
    #display(found_metadata)
    if true
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

        if_exists,found_data = check_data_exists(naming_dict,"mps";location=dataloc)
        #
        if if_exists
            psi_gs = found_data[1]["mps"]
            densmat = found_data[1]["densmat"]
        else
            dmrg_params = (ham=ham_start,mdim=mdim,if_save_data=if_save_data,metadata=metadata,name=filename,location=dataloc,observer=obs,if_densmat=if_densmat,nsweeps=nsweeps,noise=noise)
            psi_gs, densmat = execute_mps(nothing,nothing,chi,L,nflavors,part_count; dmrg_params...)
            println("Energy Variance = ",energy_variance(psi_gs,ham_start)," at Chi = ",chi)
        end
        
    end

    #append!(states,[psi_gs])

    bonddims[idx] = maxlinkdim(psi_gs)
    sf_orderparams[idx] = abs(momentum_occupation(psi_gs,1,0.0; densmat=densmat)[2][1])
    #distcorrs[idx] = minimum(abs.(distance_correlation(psi_gs; if_plot=false)[2]))
    ees[idx] = entanglement_entropy(psi_gs)
    virt_corr_lengths = physical_distance_correlation(psi_gs; if_plot=false)[3]
    for i in 1:nflavors
        corrlengs[i][idx] = virt_corr_lengths[i]
    end
    #
    
    
    
    
    #physical_distance_correlation(psi_gs)

    if true
    if idx > 1
        plot([part_count/(strens[idx-1]*nflavors*L),part_count/(chi*nflavors*L)],[corrlengs[1][idx-1],corrlengs[1][idx]],"-p",c="b")
        #plot([(part_count-1)/(strens[idx-1]*tot_sites),part_count/(alpha*tot_sites)],[centermoms[idx-1],centermoms[idx]],"-p",c="b")
    else
        #scatter([part_count/(strens[idx]*nflavors*L)],[sf_orderparams[idx]],c="b")
        #scatter([part_count/(strens[idx]*tot_sites)],[centermoms[idx]],c="b")
    end
    end
    
    if false
    nrgs[i] = calculate_energy(psi_gs,found_metadata["ham"])
    currents[:,i] = [calc_deriv(1,psi_gs,s,Int(L/2),nflavors,chi,ham_params) for s in 1:nflavors]
    drudes[:,i] = [calc_deriv(2,psi_gs,s,Int(L/2),nflavors,chi,ham_params) for s in 1:nflavors]
    end
    
    #=
    centersite = Int(ceil(L/2))
    u,s,v = svd(psi_gs[centersite],linkind(psi_gs,centersite))
    ents = sum([s[n,n]^2 * 2 * log(s[n,n]) for n in 1:size(s)[1]])
    scatter([L],[ents],c="b")
    =#
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
xlabel("Central Flux Strength")
legend()
#
fig2 = figure()
for i in 1:nflavors
    plot(strens,real.(currents[i,:]),"-p",label="Current $i")
end
xlabel("Central Flux Strength")
legend()

fig3 = figure()
for i in 1:nflavors
    plot(strens,real.(drudes[i,:]),"-p",label="Drude $i")
end
xlabel("Central Flux Strength")
legend()
=#
end


xvals = part_count ./ (strens .* (L*nflavors))
#
fig1 = figure()
scatter(xvals,bonddims)
xlabel("Filling Factor")
ylabel("Bond Dimension")

#
fig2 = figure()
plot(xvals,distcorrs,"-p")
xlabel("Filling Factor")
ylabel("Distance Correlation")

fig3 = figure()
plot(xvals,sf_orderparams,"-p")
xlabel("Filling Factor")
ylabel("SF Order Parameter")
#

fig4 = figure()
plot(xvals,ees,"-p")
plot(xvals,log.(bonddims),"-p")
xlabel("Filling Factor")
ylabel("Entanglement Entropy")

fig5 = figure()
for s in 1:nflavors
    plot(xvals,corrlengs[s],"-p",label="$s")
end
xlabel("Filling Factor")
ylabel("Correlation Length")
legend()
#
#

#
#end
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
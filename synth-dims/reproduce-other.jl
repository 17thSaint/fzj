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

if_save_data = false
dataloc = "/home/patrick/fzj/main-git/cluster-data/synth-dims/"

nrgvar_tol = 1E-12
mdim = 50

nu = 1.0
L = 8
nflavors = 5
part_count = 5
chi = part_count / (nu*L*nflavors)
tilt = 0.0

if_per_phys = true
if_per_virt = false

naming_dict = Dict([("L",L),("nflavors",nflavors),("nbosons",part_count),("chi",chi),("mdim",mdim),("centralflux_strength","n")])
metadata = merge(naming_dict,Dict([("if_periodic_phys",if_per_phys),("if_periodic_virt",if_per_virt),("tilt_strength",tilt),("location",dataloc),("if_save_data",if_save_data),("nrgvar_tol",nrgvar_tol),("mdim",mdim)]))



#current_strength = 0.00

#
if true
counting = 1
scaling = 16
strens = 0.0 .+ [sort([-i/scaling for i in 1:counting]); [0.0]; [i/scaling for i in 1:counting]]#range(0.00,stop=1.0,length=counting)
display(strens)
nrgs = zeros(length(strens)) .* im
currents = zeros(nflavors,length(strens)) .* im
drudes = zeros(nflavors,length(strens)) .* im
states = []

for (i,centralflux_strength) in enumerate(strens)
    metadata["centralflux_strength"] = centralflux_strength
    if centralflux_strength < 0.0
        naming_dict["centralflux_strength"] = "n" * string(-round(centralflux_strength,digits=5))
    else
        naming_dict["centralflux_strength"] = round(centralflux_strength,digits=5)
    end
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
        println("Energy Variance = ",energy_variance(psi_gs,ham_start)," at centralflux Strength = ",centralflux_strength)
    end

    append!(states,[psi_gs])
    nrgs[i] = calculate_energy(psi_gs,ham_start)

    currents[:,i] = [calc_deriv(1,psi_gs,s,Int(L/2),nflavors,chi,ham_params) for s in 1:nflavors]
    drudes[:,i] = [calc_deriv(2,psi_gs,s,Int(L/2),nflavors,chi,ham_params) for s in 1:nflavors]
    #=fig = figure()
    imshow(real.(drude))
    title("Drudes at Current Strength = $(centralflux_strength)")
    xlabel("Physical Site")
    ylabel("Synthetic Site")
    colorbar()
    =#
end
#
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
#
end
#


#
if false
    time_end = 10.0
    time_change = 0.2
    mdim_time = 50

    tilts = [0.1,0.075,0.05,0.025,0.01]
    for (i,tilt) in enumerate(tilts)
        ham_evolve = hamiltonian_universal(L,nflavors,chi; if_periodic_phys=if_per_phys,if_periodic_synth=if_per_virt,current_strength=current_strength,tilt_strength=tilt,if_s0=true)

        #rez0, otherham0 = evolve_in_time(psi_gs,time_end,time_change,ham_start; mdim=mdim_time,obs_measures=Dict("occs" => current_occ, "states" => return_state, "nrg_vars" => current_nrgvar, "nrgs" => current_nrg))
        rez, otherham = evolve_in_time(psi_gs,time_end,time_change,ham_evolve; mdim=mdim_time,obs_measures=Dict("occs" => current_occ, "states" => return_state))
        times = [[0.0]; rez["times"].results]

        #currents = [[jx0]; [get_current(rez["states"].results[i]; alpha=chi)[2][3] for i in 1:length(times)-1]] .-jx0
        #currents_null = [[jx0]; [get_current(rez0["states"].results[i]; alpha=chi)[2][3] for i in 1:length(times)-1]] .-jx0

        #fig = figure()
        #plot(times,-imag.(currents) ./ maximum(-imag.(currents)[1:10]),"-p",label="$(round(tilt,digits=3))")
        #legend()
    end

    #=fig2 = figure()
    plot(times,imag.(currents_null),"-p")
    title("Null")
    =#
end
#















"fin"
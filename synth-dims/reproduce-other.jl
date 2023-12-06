include("fqh_effective.jl")
include("time_evolution.jl")
using Statistics,PyPlot,Observers,ITensorTDVP,LsqFit

lin_model(x,p) = p[1].* x .+ p[2]

function hamiltonian_universal(L,nflavors,chi,tp=1.0,ts=1.0; kwargs...)
        if_periodic_phys = get(kwargs, :if_periodic_phys, false)
        if_periodic_synth = get(kwargs, :if_periodic_synth, false)
        tilt_strength = get(kwargs, :tilt_strength, 0.0)
        current_strength = get(kwargs, :current_strength, 0.0)
        if_s0 = get(kwargs, :if_s0, false)
        
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
                ampo += (-tp * exp(im*chi*(s-s0)/(nflavors)) * exp(im*current_strength/L), "Cr$s", j, "Anh$s", next_site)
                ampo += (-tp * exp(-im*chi*(s-s0)/(nflavors)) * exp(-im*current_strength/L), "Anh$s", j, "Cr$s", next_site)
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

nrgvar_tol = 1E-8
mdim = 50

nu = 1.0
L = 8
nflavors = 5
part_count = 5
chi = part_count / (nu*L*nflavors)
#tilt = 0.001

if_per_phys = true
if_per_virt = true

#current_strength = 0.00

#
if true
change = 0.001
counting = 5
strens = range(0.2,stop=0.3,length=counting)
strens = [strens; strens .+ change]
nrgs = zeros(length(strens)) .* im
currents = zeros(Int(length(strens)/2)) .* im
states = []

for (i,current_strength) in enumerate(strens)
ham_start = hamiltonian_universal(L,nflavors,chi; if_periodic_phys=if_per_phys,if_periodic_synth=if_per_virt,current_strength=current_strength,if_s0=true,tilt_strength=0.0)
obs = NRGVarObserver(nrgvar_tol,ham_start)
psi_gs = execute_mps(nothing,nothing,chi,L,nflavors,part_count; ham=ham_start,mdim=mdim,if_save_data=if_save_data,observer=obs)
println("Energy Variance = ",energy_variance(psi_gs,ham_start))

#jx0 = get_current(psi_gs; alpha=chi)[2][3]
append!(states,[psi_gs])
#occ0 = get_occupancy(psi_gs)
if (i+1) % 2 == 0
jx = get_current(psi_gs; alpha=chi, if_exp_part=true)
currents[Int((i+1)/2)] = jx[2][3]#jx[1]
end
nrgs[i] = calculate_energy(psi_gs,ham_start)
end
scatter(strens,real.(nrgs),"-p")
#
#jx_theory = real.([(nrgs[i+count] - nrgs[i])/change for i in 1:count]) #.* L
#new_strens = [strens[i] + change/2 for i in 1:count]
#
#plot(strens[1:count],imag.(currents),"-p")
#plot(new_strens,jx_theory,label="Theory")
#xlabel("Phi")
#ylabel("Current")
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

currents = [[jx0]; [get_current(rez["states"].results[i]; alpha=chi)[2][3] for i in 1:length(times)-1]] .-jx0
#currents_null = [[jx0]; [get_current(rez0["states"].results[i]; alpha=chi)[2][3] for i in 1:length(times)-1]] .-jx0

#fig = figure()
plot(times,-imag.(currents) ./ maximum(-imag.(currents)[1:10]),"-p",label="$(round(tilt,digits=3))")
legend()
end

#=fig2 = figure()
plot(times,imag.(currents_null),"-p")
title("Null")
=#
end
#















"fin"
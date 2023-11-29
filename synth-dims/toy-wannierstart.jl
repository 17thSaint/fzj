include("fqh_effective.jl")
include("time_evolution.jl")

function wannierstark_ham(L::Int,field_strength::Number; kwargs...)
    if_periodic = get(kwargs, :if_periodic, false)
    hopping_strength = get(kwargs, :hopping_strength, 1.0)

    ampo = OpSum()

    # hopping term
    for j in 1:L
        next_site = j+1
        if j == L
            if if_periodic
                next_site = 1
            else
                continue
            end
        end
        ampo += (hopping_strength, "adag", j, "a", next_site)
        ampo += (conj(hopping_strength), "adag", next_site, "a", j)
    end

    # applied field term
    for j in 1:L
        ampo += (-field_strength*j, "N", j)
    end

    return ampo
end

function avg_pos(psi)
    occs = expect(psi,"N")
    occs ./= sum(occs)
    result = 0.0
    for (i,o) in enumerate(occs)
        result += i*o
    end
    return result
end

function current_avgpos(; psi, bond, half_sweep, current_time)
	if bond == 1 && half_sweep == 2
	  pos = avg_pos(psi)
      scatter([-imag(current_time)],[pos],c="r")
      return pos
	end
	return nothing
end

L = 10
mdim = 100
mdim_time = 100
nbosons = 3
time_end = 200.0

ham_start = wannierstark_ham(L,0.0)

sites = siteinds("Boson", L)
states = make_states(L,nbosons,1)
psi0 = MPS(sites,states)
gs_psi = execute_mps(nothing,nothing,nothing,L,nothing,nbosons; psi_guess=psi0,ham=ham_start,mdim=mdim)
println("Initial Energy Variance = ",energy_variance(gs_psi,MPO(ham_start,siteinds(gs_psi))))

#=occs = expect(gs_psi,"N")
fig = figure()
plot(collect(1:L),occs,"-p")
xlabel("Position")
ylabel("Occupation")
title("Initial Occupation")
=#

time_change = 0.5

estrength = 0.01
ham_evolve = wannierstark_ham(L,estrength)
rez, last_ham = evolve_in_time(gs_psi,time_end,time_change,ham_evolve; mdim=mdim_time,obs_measures=Dict("current_pos" => current_avgpos))
times = rez["times"].results

#=
fig2 = figure()
plot(times,rez["current_pos"].results,"-p")
xlabel("Time")
ylabel("Position")
=#

nrg_var = [energy_variance(rez["states"].results[i],MPO(ham_evolve,siteinds(rez["states"].results[i]))) for i in 1:length(times)]
nrgvars[i] = mean(nrg_var)
println("Average Energy Variance = ",mean(nrg_var))
#end











"fin"
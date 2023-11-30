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
      #scatter([-imag(current_time)],[pos],c="r")
      return pos
	end
	return nothing
end

function get_spectrum(number_of_levels,ham,psi0,L,nbosons; kwargs...)
    mdim = get(kwargs, :mdim, 100)

    gs_psi = execute_mps(nothing,nothing,nothing,L,nothing,nbosons; psi_guess=psi0,ham=ham,mdim=mdim)
    states = [gs_psi]
    energies = [calculate_energy(gs_psi,ham)]
    psi_ortho = [gs_psi]
    for i in 2:number_of_levels
        println("Working on level $i")
        wavefunc = execute_mps(nothing,nothing,nothing,L,nothing,nbosons; psi_ortho=psi_ortho,psi_guess=psi0,ham=ham_start,mdim=mdim)
        append!(states,[wavefunc])
        append!(energies,[calculate_energy(wavefunc,ham)])
        append!(psi_ortho,[wavefunc])
    end
    return states,energies
end

function positional_nrg(L,pos,estregnth=0.01)
    sites = siteinds("Boson", L)
    states = ["0" for i in 1:L]
    states[pos] = "1"
    ham_evolve = wannierstark_ham(L,estrength)
    psi = MPS(sites,states)
    return calculate_energy(psi,ham_evolve)
end

function find_decay_rate(positions,times; kwargs...)
    if_plot = get(kwargs, :if_plot, true)
    model(x,p) = p[5] .* exp.(p[1] .* x) .* (sin.(p[2] .* x .+ p[4])).^1 .+ p[3]
    fit = curve_fit(model,times,positions,[-1.0,1.0,positions[1],0.0,positions[1]])
    decay_rate = -1/fit.param[1]
    freq = -fit.param[2]
    if if_plot
        fig = figure()
        scatter(times,positions,c="b")
        plot(times,model(times,fit.param),c="r")
        xlabel("Time")
        ylabel("Position")
        title("Decay Rate = $(round(decay_rate,digits=4)), with Freq = $(round(freq,digits=4))")
    end
    return decay_rate,freq
end

#
L = 10
mdim = 100
mdim_time = 100
nbosons = 1
time_end = 50.0

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
#=
e_start = 0.001
e_end = 0.2
e_count = 50
strens = [e_start + (i-1)*(e_end-e_start)/(e_count-1) for i in 1:e_count]
rates = [0.0 for i in 1:length(strens)]
freqs = [0.0 for i in 1:length(strens)]
time_change = 0.5
for (i,estrength) in enumerate(strens)
#estrength = 0.1
ham_evolve = wannierstark_ham(L,estrength)
time_gs_psi = execute_mps(nothing,nothing,nothing,L,nothing,nbosons; psi_guess=gs_psi,ham=ham_evolve,mdim=mdim)

rez, last_ham = evolve_in_time(gs_psi,time_end,time_change,ham_evolve; mdim=mdim_time,obs_measures=Dict("current_pos" => current_avgpos))
times = rez["times"].results

decay_params = find_decay_rate(rez["current_pos"].results,times; if_plot=false)
rates[i] = decay_params[1]
freqs[i] = decay_params[2]
end
=#
fig = figure()
plot(strens,abs.(rates),"-p")
yscale("log")
xlabel("Field Strength")
ylabel("Decay Time")

fig2 = figure()
plot(strens,abs.(freqs),"-p")
xlabel("Field Strength")   
ylabel("Frequency")

#=
fig2 = figure()
plot(times,rez["current_pos"].results,"-p")
xlabel("Time")
ylabel("Position")
=#

#nrg_var = [energy_variance(rez["states"].results[i],MPO(ham_evolve,siteinds(rez["states"].results[i]))) for i in 1:length(times)]
#nrgvars[i] = mean(nrg_var)
#println("Average Energy Variance = ",mean(nrg_var))
#end











"fin"
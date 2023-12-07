include("fqh_effective.jl")
include("time_evolution.jl")
using Statistics,PyPlot,Observers,ITensorTDVP,LsqFit

function wannierstark_ham(L::Int,field_strength::Number; kwargs...)
    if_periodic = get(kwargs, :if_periodic, true)
    hopping_strength = get(kwargs, :hopping_strength, 1.0)
    current_strength = get(kwargs, :current_strength, 0.0)
    if_current = current_strength != 0.0
    if_field = field_strength != 0.0

    println("Whether or not Periodic = ",if_periodic)

    if_current ? hopping_strength *= exp(im*2*pi*current_strength/L) : nothing

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

    if if_field
        for j in 1:L
            ampo += (-field_strength*j, "N", j)
        end
    end

    return ampo
end

function momentum_ham(L,nbosons,current_strength; kwargs...)
    ampo = OpSum()

    for i in 0:L-1
        coeff = cos((pi*nbosons - current_strength + 2*pi*i)/L) / L
        ampo += (coeff, "adag", i+1, "a", i+1)
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

function toy_current(psi,phi,L)
	fullmat = correlation_matrix(psi,"adag","a")
    #fig = figure()
    #imshow(abs.(fullmat))
	currents = zeros(L) .* im
    for i in 1:L
        minus = i-1
        if i == 1
            minus = L
        end
        plus = i+1
        if i == L
            plus = 1
        end
        inflow = exp(im*phi*2*pi/L) * fullmat[plus,i] + exp(-im*phi*2*pi/L) * fullmat[minus,i]
        outflow = exp(-im*phi*2*pi/L) * fullmat[i,minus] + exp(im*phi*2*pi/L) * fullmat[i,plus]
        #display(inflow)
        #display(outflow)
        currents[i] = inflow - outflow #exp(im*phi*2*pi/L) * (fullmat[plus,i] - fullmat[i,minus]) + exp(-im*phi*2*pi/L) * (fullmat[minus,i] - fullmat[i,plus])
    end
    
    return -imag.(currents)
end

function working_ham(L,phi)
    ampo = OpSum()
    coeff = exp(im*2*pi*phi/L)
    for i in 1:L
        ampo += (coeff, "S+", i, "S-", mod1(i+1,L))
        ampo += (conj(coeff), "S+", mod1(i+1,L), "S-", i)
    end
    return ampo
end

function current_calculate(psi::MPS,site::Int,L,phi)
    coeff = exp(im*phi*2*pi/L)
    sites = siteinds(psi)
    opform = im*coeff*op(sites[site],"S+") * op(sites[mod1(site+1,L)],"S-") - im*conj(coeff)*op(sites[mod1(site+1,L)],"S+") * op(sites[site],"S-")
    
    result = inner(psi,apply(opform,psi))
    return result
end

function current_calculate(psi::MPS,L,phi)
    return [current_calculate(psi,i,L,phi) for i in 1:L]
end

function differentiate_state(psi::MPS,psi_shifted::MPS,shift)
    return abs2.((inner(psi,psi_shifted) - 1) / shift)
end

function differentiate_state(psi::Vector,psi_shifted::Vector,shift)
    return [differentiate_state(psi[i],psi_shifted[i],shift) for i in 1:length(psi)]
end

linmod(x,p) = p[1] .* x .+ p[2]

function find_deriv_overlap(psi::MPS,observer,local_strength,psi0)
    shifted_states = []
    new_strengths = []
    derivoverlaps = []
    new_nrgs = []
    change = local_strength/100
    change_limit = local_strength/10
    next_strength = local_strength + change
    reshiftcount = 0
    fitcount = 5
    attempts = 1
    while length(shifted_states) < 2*fitcount
        if next_strength > local_strength + change_limit
            reshiftcount += 1
            # this calculates new strengths which are shifted by 1/2, 1/4, 1/8 if havent found enough good fits
            next_strength = local_strength + change*(0.5^reshiftcount) + change
        end
        next_ham = working_ham(length(psi),next_strength)
        next_psi = execute_mps(nothing,nothing,nothing,nothing,nothing,nothing; psi_guess=psi0,ham=next_ham,mdim=100,outputlevel=0)#,observer=observer)
        derivoverlap = differentiate_state(psi,next_psi,next_strength-local_strength)
        println(next_strength,", ",local_strength,", ",derivoverlap,", ",reshiftcount)
        if derivoverlap <= 2.0
            append!(shifted_states,[next_psi])
            append!(new_strengths,[next_strength])
            append!(derivoverlaps,[derivoverlap])
            append!(new_nrgs,[calculate_energy(next_psi,next_ham)])

            # need mirrored strength to get second derivative of energy
            otherdirection_strength = 2*local_strength - next_strength
            otherdirection_ham = working_ham(length(psi),otherdirection_strength)
            otherdirection_psi = execute_mps(nothing,nothing,nothing,nothing,nothing,nothing; psi_guess=psi0,ham=otherdirection_ham,mdim=100,outputlevel=0)#,observer=observer)
            append!(shifted_states,[otherdirection_psi])
            append!(new_strengths,[otherdirection_strength])
            append!(new_nrgs,[calculate_energy(otherdirection_psi,otherdirection_ham)])
        end
        attempts += 1
        next_strength += change
    end
    forward_strengths = [new_strengths[2*i+1] for i in 0:fitcount-1]
    cf = curve_fit(linmod,forward_strengths .- local_strength,derivoverlaps,[-1.0,0.5])
    fig = figure()
    scatter(forward_strengths .- local_strength,derivoverlaps)
    plot(forward_strengths .- local_strength,linmod(forward_strengths .- local_strength,cf.param))
    title("Overlap: Limit Value = $(cf.param[2])")
    println("Found Limit Derivative with $(attempts-fitcount) extra attempts and $(reshiftcount) reshifts")
    return shifted_states,new_strengths,new_nrgs,cf.param[2]
end

function nrg_limitderiv(nrgs,real_strengths,og_strengths,order)
    fitcount = 5
    nrg_sets = [nrgs[(2*fitcount+1)*i+1:(2*fitcount+1)*i+1+2*fitcount] for i in 0:length(og_strengths)-1]
    strength_sets = [real_strengths[(2*fitcount+1)*i+1:(2*fitcount+1)*i+1+2*fitcount] for i in 0:length(og_strengths)-1]
    if order == 1
        values = []
        for (i,nrg_set) in enumerate(nrg_sets)
            forward_nrgs = [[nrg_set[1]]; [nrg_set[2*j] for j in 1:fitcount]]
            forward_strengths = [[strength_sets[i][1]]; [strength_sets[i][2*j] for j in 1:fitcount]]
            append!(values,first_deriv_nrgs(forward_nrgs,forward_strengths))
        end
    elseif order == 2
        println("Not done yet")
        values = []
        for (i,nrg_set) in enumerate(nrg_sets)
            center_nrg = nrg_set[1]
            center_strength = strength_sets[i][1]
            forward_nrgs = [nrg_set[2*j] for j in 1:fitcount]
            forward_strengths = [strength_sets[i][2*j] for j in 1:fitcount]
            reverse_nrgs = [nrg_set[2*j+1] for j in 1:fitcount]
            append!(values,second_deriv_nrgs(forward_nrgs,forward_strengths,reverse_nrgs,center_nrg,center_strength))
        end
    end

    return values
end

function second_deriv_nrgs(forward_nrgs,forward_strens,reverse_nrgs,center_nrg,center_stren)
    full_derivs = (forward_nrgs .- (2 .* center_nrg) .+ reverse_nrgs) ./ (forward_strens .- center_stren).^2
    derivs = real.(full_derivs)
    println("Imag Part: ",mean(imag.(derivs)))
    cf = curve_fit(linmod,forward_strens .- center_stren,derivs,[-1.0,0.5])
    fig = figure()
    scatter(forward_strens .- center_stren,derivs)
    plot(forward_strens .- center_stren,linmod(forward_strens .- center_stren,cf.param))
    title("Second Derivative")
    return cf.param[2]
end

function first_deriv_nrgs(nrgs,strengths)
    full_derivs = (nrgs[2:end] .- nrgs[1]) ./ (strengths[2:end] .- strengths[1])
    derivs = real.(full_derivs)
    println("Imag Part: ",mean(imag.(derivs)))
    cf = curve_fit(linmod,strengths[2:end] .- strengths[1],derivs,[-1.0,0.5])
    fig = figure()
    scatter(strengths[2:end] .- strengths[1],derivs)
    plot(strengths[2:end] .- strengths[1],linmod(strengths[2:end] .- strengths[1],cf.param))
    title("First Derivative")
    return cf.param[2]
end


        

#
#L = 10
mdim = 100
mdim_time = 100
nbosons = 1
time_end = 50.0
nrgvar_tol = 1e-8

#change = 0.01
counting = 2
strens = range(0.1,stop=0.25,length=counting) #0.25 .+ [i/1000 for i in 0:100]
#strens = [strens; strens .+ change]
Ls = [10]
real_strens = [[] for i in 1:length(Ls)]
nrgs = [[] for i in 1:length(Ls)]
direct_jx = [[] for i in 1:length(Ls)]
secderiv_nrg = [[] for i in 1:length(Ls)]
deriv_jx = [[] for i in 1:length(Ls)]
states = [[] for i in 1:length(Ls)]
deriv_states = [[] for i in 1:length(Ls)]
theory = [[] for i in 1:length(Ls)]
fromham = [[] for i in 1:length(Ls)]

for (j,L) in enumerate(Ls)
sites = siteinds("Qubit", L; conserve_number = true)
making_states = make_states(L,nbosons,1)
psi0 = MPS(sites,making_states)
#nrgs[j] = zeros(length(strens)) .* im
#direct_jx[j] = zeros(length(strens)) .* im
for (i,stren) in enumerate(strens)
    append!(real_strens[j],[stren])
    hamhere = working_ham(L,stren)
    obs = NRGVarObserver(nrgvar_tol,hamhere)

    gs_psi = execute_mps(nothing,nothing,nothing,nothing,nothing,nothing; psi_guess=psi0,ham=hamhere,mdim=mdim)#,observer=obs)
    append!(states[j],[gs_psi])
    append!(nrgs[j],[calculate_energy(gs_psi,hamhere)])
    #direct_jx[j][i] = current_calculate(gs_psi,Int(L/2),L,stren)
    #
    new_states,new_strens,new_nrgs,deriv_overlap = find_deriv_overlap(gs_psi,obs,stren,psi0)
    append!(real_strens[j],new_strens)
    append!(states[j],new_states)
    append!(deriv_states[j],[deriv_overlap])
    append!(nrgs[j],new_nrgs)
    #
end
#
deriv_jx[j] = nrg_limitderiv(nrgs[j],real_strens[j],strens,1)
secderiv_nrg[j] = (L/(2*pi))^2 .* nrg_limitderiv(nrgs[j],real_strens[j],strens,2)
#deriv_states[j] = differentiate_state(states[j][1:counting],states[j][2*counting+1:end],change)
#
end
#

#rr = [differentiate_state(states[1][1],states[1][i],strens[i] - strens[1]) for i in 2:length(strens)]
#scatter(strens[2:end],rr)

#
if true
fig = figure()
for (j,L) in enumerate(Ls)
    #fig1 = figure()
scatter(strens,real.(deriv_jx[j]),label="Deriv $L")
#scatter(strens,real.(direct_jx[j]),label="Direct $L")
legend()
xlabel("Phi")
ylabel("Jx")
end
end

if true
fig = figure()
for (j,L) in enumerate(Ls)
#fig2 = figure()
scatter(real_strens,real.(nrgs[j]),label="$L")
xlabel("Phi")
ylabel("Energy")
legend()
end
end

if true
fig = figure()
for (j,L) in enumerate(Ls)
#secderiv_nrg[j] = L^2 .* ((nrgs[j][2*counting+1:end] .- (2 .* nrgs[j][counting+1:2*counting]) .+ nrgs[j][1:counting]) / (change/2)^2 / (2*pi)^2)
#fig3 = figure()
scatter(strens,real.(secderiv_nrg[j]),label="$L")
xlabel("Phi")
ylabel("Drude Weight")
legend()
end
end

if false
fig4 = figure()
ccs = [0.0 for i in 1:length(Ls)]
for (j,L) in enumerate(Ls)
#fig4 = figure()
theory[j] = (secderiv_nrg[j] .* (2*pi/L)^2) .+ ((2. * nrgs[j][counting+1:2*counting]) .* (deriv_states[j] .- 1))
fromham[j] = nrgs[j][1:counting]
what = filter(x -> x < 0, real.(theory[j] ./ fromham[j]))
display(what)
ccs[j] = mean(what)
scatter(strens[1:counting],real.(fromham[j]),label="From Ham $L")
scatter(strens[1:counting],real.(theory[j]),label="Theory $L")
xlabel("Phi")
ylabel("Exp of SecDeriv Hamiltonian")
legend()
end
end
#
#=
torus_current = 0.01
counting = 50
change = 0.0001
strens = range(0.01,stop=2.0,length=counting)
derivstrens = strens .+ change/2
#strens = [strens; strens .+ change]
nrgs = zeros(length(strens)) .* im
jxs = zeros(length(strens))
#direct_jxs = zeros(Int(length(strens)/2)) .* im
for (i,torus_current) in enumerate(strens)
ham_start = wannierstark_ham(L,0.0; current_strength=torus_current, if_periodic=true)
obs = NRGVarObserver(nrgvar_tol,ham_start)

gs_psi = execute_mps(nothing,nothing,nothing,L,nothing,nbosons; psi_guess=psi0,ham=ham_start,mdim=mdim,observer=obs,if_periodic=true)
#println("Initial Energy Variance = ",energy_variance(gs_psi,ham_start))
nrgs[i] = calculate_energy(gs_psi,ham_start)
loccurr = toy_current(gs_psi,torus_current,L)
#display(loccurr)
#println("Mean values = ",mean(loccurr))
#plot(loccurr)
jxs[i] = loccurr[Int(L/2)]
end

#=
for i in 1:counting
    direct_jxs[i] = (nrgs[i+counting] - nrgs[i])/change
end
plot(derivstrens,imag.(direct_jxs),"-p",label="Deriv")
=#
#
plot(strens,jxs,"-p",label="Operator")
xlabel("Current Strength")
ylabel("Jx")
legend()
=#

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

#=
fig = figure()
plot(strens,abs.(rates),"-p")
yscale("log")
xlabel("Field Strength")
ylabel("Decay Time")

fig2 = figure()
plot(strens,abs.(freqs),"-p")
xlabel("Field Strength")   
ylabel("Frequency")
=#

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
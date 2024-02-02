if true
#include("../other-funcs/data-storage-funcs.jl")
#include("long-range-ttn.jl")
#include("fqh_effective.jl")
using Observers, NumericalIntegration
using ITensorTDVP
end

function current_nrg(; psi, bond, half_sweep)
    if bond == 1 && half_sweep == 2
      return calculate_energy(psi,ham_evolve)
    end
    return nothing
end

function current_nrgvar(; psi, bond, half_sweep)
    if bond == 1 && half_sweep == 2
      nrgvar = energy_variance(psi,ham_evolve)
      println("Current Energy Variance = ",nrgvar)
      return nrgvar
    end
    return nothing
end

function get_time(; current_time, bond, half_sweep)
    if bond == 1 && half_sweep == 2
	  return -imag(current_time)
	end
	return nothing
end

function return_state(; psi, bond, half_sweep)
	if bond == 1 && half_sweep == 2
	  return psi
	end
	return nothing
end

function return_state_ttn(psi)
    return psi
end

function current_occ(; psi, bond, half_sweep, current_time)
	if bond == 1 && half_sweep == 2
      #if isapprox(-imag(current_time) % 2,0.0;atol=10^-3) || isapprox(-imag(current_time) % 2,2.0;atol=10^-3)
       #   return get_occupancy(psi; plot_title="$(-imag(current_time))")
      #else
      #    println("No Plot: ",-imag(current_time) % 2)
    	  return get_occupancy(psi; if_plot=false)
      #end
	end
	return nothing
end

function current_density_polarization(; psi, bond, half_sweep, current_time)
    if bond == 1 && half_sweep == 2
      result = density_polarization(psi)
      #println("Result: ",result)
      scatter([-imag(current_time)],[result],c="r")
      return result
    end
    return nothing
end

function current_spacial_density_polarization(; psi, bond, half_sweep, current_time)
    if bond == 1 && half_sweep == 2
      result = spacial_density_polarization(psi)
      #println("Result: ",result)
      scatter([-imag(current_time)],[result],c="b")
      return result
    end
    return nothing
end

function average_position(psi,occmat=nothing)
	occmat = isnothing(occmat) ? get_occupancy(psi; if_plot=false) : occmat
	xs = [-Int(size(occmat)[1]/2) + (i-1)*size(occmat)[1]/(size(occmat)[1]-1) for i in 1:size(occmat)[1]]
	ys = [-Int(size(occmat)[2]/2) + (i-1)*size(occmat)[2]/(size(occmat)[2]-1) for i in 1:size(occmat)[2]]
	n1 = integrate((xs,ys),occmat)
	occmat ./= n1
	position_matrix = zeros(size(occmat))
	avg_position = [0.0,0.0]
	for i in 1:size(occmat)[1]
		for j in 1:size(occmat)[2]
			avg_position += [i,j] .* occmat[i,j]
		end
	end
	return avg_position
end

function average_position(all_occmats::Vector{Matrix{Float64}})
    return [average_position(nothing,occmat) for occmat in all_occmats]
end

function density_polarization(psi,occmat=nothing)
    occmat = isnothing(occmat) ? get_occupancy(psi; if_plot=false) : occmat
    result = 0.0
    m0 = (size(occmat)[2]+1)/2
    for i in 1:size(occmat)[2]
        result += 2*sum(occmat[:,i] .* (i - m0))
    end
    result /= sum(occmat)*size(occmat)[2]
    return result
end

function density_polarization(all_occmats::Vector{Matrix{Any}})
    return [density_polarization(nothing,occmat) for occmat in all_occmats]
end

function spacial_density_polarization(psi,occmat=nothing)
    occmat = isnothing(occmat) ? get_occupancy(psi; if_plot=false) : occmat
    result = 0.0
    j0 = (size(occmat)[1]+1)/2
    virtual_results = [0.0 for i in 1:size(occmat)[2]]
    for s in 1:size(occmat)[2]
        virtual_results[s] = 2*sum(occmat[:,s] .* [i - j0 for i in 1:size(occmat)[1]])
    end
    for i in 1:size(occmat)[1]
        result += sum(occmat[i,:] .* (i - j0))
    end
    return 2 * result / (sum(occmat)*size(occmat)[1]),virtual_results
end

function spacial_density_polarization(all_occmats::Vector{Matrix{Float64}})
    return [spacial_density_polarization(nothing,occmat) for occmat in all_occmats]
end

function calculate_energy(psi::MPS,H)
    if typeof(H) != MPO
        H = MPO(H,siteinds(psi))
    end
	return inner(psi', H, psi) / inner(psi, psi)
end

function calculate_energy(psis::Vector,H)
    return [calculate_energy(psi,H) for psi in psis]
end

function energy_variance(psi,H)
    if typeof(H) != MPO
        H = MPO(H,siteinds(psi))
    end
    fo_nrg = calculate_energy(psi,H)
	return abs(fo_nrg^2 - inner(H,psi,H,psi) / inner(psi,psi))
end

mutable struct TTNObserver <: AbstractObserver
    energy::Vector{Float64}
    times::Vector{Float64}
    measurement_operators::Dict{String,Any}
    measurement_results::Dict{String,Any}

    function TTNObserver(what_to_measure::Dict{String,Any})
        measurements = Dict(key => Any[] for key in keys(what_to_measure))
        return new(Float64[], Float64[], what_to_measure, measurements)
    end
end

function ITensors.measure!(ob::TTNObserver; kwargs...)
    tdvp = kwargs[:sweep_handler]
    time_now = tdvp.current_time


    #((tdvp.dirloop !== :backward) || (pos !== tdvp.path[2])) && return

    topPos = (TTNKit.number_of_layers(TTNKit.network(tdvp.ttn)), 1)
    n_sites = TTNKit.number_of_sites(TTNKit.network(tdvp.ttn))

    action = TTNKit.∂A(tdvp.pTPO, topPos)
    T = tdvp.ttn[topPos]
    actionT = action(T)

    push!(ob.energy, real(ITensors.scalar(dag(T)*actionT)))
    push!(ob.times, time_now)

    for (key,val) in ob.measurement_operators
        params = ()
        func = val
        try
            if length(val) > 1
                func = val[1]
                params = val[2]
            else
                func = val[1]
            end
        catch
            func = val
        end
        result = func(tdvp.ttn; params...)
        if typeof(result) == Vector
            append!(ob.measurement_results[key],result)
        else
            append!(ob.measurement_results[key],[result])
        end
    end

end


function evolve_in_time(psi0,final_time,dt,ham; kwargs...)
    obs_measures::Dict{String,Any} = get(kwargs,:obs_measures,Dict{String,Any}())
    outputlevel = get(kwargs,:outputlevel,1)
    mdim = get(kwargs,:mdim,10)
    if_save_data = get(kwargs,:if_save_data,false)
    if_states = get(kwargs,:if_states,true)

    metadata = get(kwargs, :metadata, Dict())
    metadata["time_end"] = final_time
    metadata["dt"] = dt
    metadata["psi0"] = psi0
    metadata["time_ham"] = ham
    #display(metadata)

    time_steps = Int(ceil(final_time/dt))
    println("Time Steps: $time_steps")

    if typeof(psi0) == MPS
        H = MPO(ham,siteinds(psi0))
    else
        lat = TTNKit.physical_lattice(TTNKit.network(psi0))
        H = TTNKit.TPO(ham,lat)
    end
    println("Made Hamiltonian")

    #obs = get(kwargs, :observer, NoObserver())
    #
    if typeof(psi0) == MPS
        obs = Observer("times" => get_time)
        if_states ? obs["states"] = return_state : nothing
        if obs_measures != nothing
            for (key,val) in obs_measures
                obs[key] = val
            end
        end
    else
        if_states ? obs_measures["states"] = (return_state_ttn,()) : nothing
        obs = TTNObserver(obs_measures)
    end
    #obs = Observer("times" => current_time)#,"states" => return_state)#, "occs" => current_occ)
    println("Made Observer")
    #

    time_start = time()
    if typeof(psi0) == MPS
        psit = tdvp(H,psi0,-im*dt; nsweeps=time_steps,outputlevel=outputlevel,maxdim=mdim,(observer!)=obs)
    else
        time_ttnswp = TTNKit.tdvp(psi0,H; timestep=dt,finaltime=final_time-dt, observer=obs)
    end
    time_end = time()

    if if_save_data
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "tevo")
		filename = check_plot_label(filename,"tevo")
        metadata["runtime"] = time_end - time_start
		data_dict = Dict([("observer",obs)])
		write_data_jld2(filename,data_dict,location,metadata)
	end

    if typeof(psi0) == MPS
        return obs,H
    else
        return obs,time_ttnswp
    end
end

function check_if_GS(psi0,ham,metadata)
   new_gs = execute_mps(metadata["U1"],metadata["U2"],metadata["phi"],metadata["L"],metadata["nflavors"],metadata["nbosons"]; psi_guess=psi0,ham=ham,dict_to_symbols(metadata)...,if_save_data=false,mdim=metadata["maxlinkdim"],outputlevel=1)
   println("Initial Energy Variance = ",energy_variance(new_gs,MPO(ham,siteinds(new_gs))))
   println("Initial Occupation Variance = ",first(occupancy_variance(new_gs; if_plot=false)))
   return new_gs
end

function execute_tevo(psi0_filename,final_time,dt; kwargs...)
    location = get(kwargs, :location, pwd())
    if_current = get(kwargs, :if_current, false)
    if if_current
        current_strength = get(kwargs, :current_strength, 0.0)
    else
        current_strength = 0.0
    end
    mdim = get(kwargs, :mdim, 10)
    if_GScheck = get(kwargs, :if_GScheck, false)

    data0,metadata0 = read_data_jld2(psi0_filename,location)
    local_phi = 0.0#metadata0["phi"]
    psi0 = data0["mps"]
    if if_GScheck
        static_ham = hamiltonian(metadata0["t1"],metadata0["t2"],local_phi,metadata0["U1"],metadata0["U2"],metadata0["L"],metadata0["nflavors"]; dict_to_symbols(metadata0)...)
        psi0 = check_if_GS(psi0,static_ham,metadata0)
    end


    new_metadata = Dict([("if_applied_current",if_current),("current_strength",current_strength),("time_mdim",mdim)])
    metadata = merge(metadata0,new_metadata)
    metadata["phi"] = local_phi
    metadata["alpha"] = round(local_phi/(2*pi),digits=4)
    naming_dict = merge(get_params_dict_from_filename(psi0_filename),Dict([("if_current",if_current),("current_strength",current_strength)]))
    naming_dict["alpha"] = round(local_phi/(2*pi),digits=4)

    time_ham = hamiltonian(metadata0["t1"],metadata0["t2"],local_phi,metadata0["U1"],metadata0["U2"],metadata0["L"],metadata0["nflavors"]; dict_to_symbols(metadata0)...,if_applied_current=if_current,current_strength=current_strength)

    naming_dict["dt"] = dt
	naming_dict["time_end"] = final_time
	filename = make_parameters_filename(naming_dict)
	println(filename)

    if_exists,data = check_data_exists("tevo-"*filename*".jld2","observer";location=location)
    if if_exists
        return data
    end

	tevo_params = (metadata=metadata,location=loc,name=filename)

	obs,time_ham = evolve_in_time(psi0,final_time,dt,time_ham; tevo_params...,kwargs...)
    return obs,time_ham
end










































"fin"
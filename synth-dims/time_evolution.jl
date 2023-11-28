if false
include("../other-funcs/data-storage-funcs.jl")
include("long-range-ttn.jl")
include("fqh_effective.jl")
#using ITensorsTDVP, Observers
end

function current_time(; current_time, bond, half_sweep)
	if bond == 1 && half_sweep == 2
	  return current_time
	end
	return nothing
end

function return_state(; psi, bond, half_sweep)
	if bond == 1 && half_sweep == 2
	  return psi
	end
	return nothing
end

function current_occ(; psi, bond, half_sweep)
	if bond == 1 && half_sweep == 2
	  return get_occupancy(psi; if_plot=false)
	end
	return nothing
end

function current_density_polarization(; psi, bond, half_sweep, current_time)
    if bond == 1 && half_sweep == 2
      result = density_polarization(psi)
      #println("Result: ",result)
      scatter([current_time],[result],c="r")
      return result
    end
    return nothing
end

function current_spacial_density_polarization(; psi, bond, half_sweep, current_time)
    if bond == 1 && half_sweep == 2
      result = spacial_density_polarization(psi)
      println("Result: ",result)
      scatter([current_time],[result],c="b")
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
    result /= size(occmat)[2]
    return result
end

function density_polarization(all_occmats::Vector{Matrix{Float64}})
    return [density_polarization(nothing,occmat) for occmat in all_occmats]
end

function spacial_density_polarization(psi,occmat=nothing)
    occmat = isnothing(occmat) ? get_occupancy(psi; if_plot=false) : occmat
    result = 0.0
    j0 = (size(occmat)[1]+1)/2
    for i in 1:size(occmat)[1]
        result += sum(occmat[i,:] .* (i - j0))
    end
    return 2 * result / size(occmat)[1]
end

function spacial_density_polarization(all_occmats::Vector{Matrix{Float64}})
    return [spacial_density_polarization(nothing,occmat) for occmat in all_occmats]
end

function evolve_in_time(psi0,final_time,dt,ham; kwargs...)
    obs_measures = get(kwargs,:obs_measures,nothing)
    outputlevel = get(kwargs,:outputlevel,1)
    mdim = get(kwargs,:mdim,10)
    if_save_data = get(kwargs,:if_save_data,false)

    metadata = get(kwargs, :metadata, Dict())
    metadata["time_end"] = final_time
    metadata["dt"] = dt
    metadata["psi0"] = psi0
    metadata["time_ham"] = ham
    display(metadata)

    time_steps = Int(ceil(final_time/dt))
    println("Time Steps: $time_steps")

    obs = Observer("states" => return_state, "times" => current_time, "occs" => current_occ)
    if obs_measures != nothing
        for (key,val) in obs_measures
            obs[key] = val
        end
    end
    println("Made Observer")

    H = MPO(ham,siteinds(psi0))
    println("Made Hamiltonian")

    time_start = time()
    psit = tdvp(H,psi0,dt; nsweeps=time_steps,outputlevel=outputlevel,maxdim=mdim,(observer!)=obs)
    time_end = time()

    if if_save_data
		location = get(kwargs, :location, pwd())
		filename = get(kwargs, :name, "tevo")
		filename = check_plot_label(filename,"tevo")
        metadata["runtime"] = time_end - time_start
		data_dict = Dict([("observer",obs)])
		write_data_jld2(filename,data_dict,location,metadata)
	end

    return obs
end

function execute_tevo(psi0_filename,final_time,dt; kwargs...)
    location = get(kwargs, :location, pwd())
    if_current = get(kwargs, :if_current, false)
    current_strength = get(kwargs, :current_strength, 0.0)
    mdim = get(kwargs, :mdim, 10)

    data0,metadata0 = read_data_jld2(psi0_filename,location)
    psi0 = data0["mps"]
    new_metadata = Dict([("if_applied_current",if_current),("current_strength",current_strength),("time_mdim",mdim)])
    metadata = merge(metadata0,new_metadata)
    naming_dict = merge(get_params_dict_from_filename(psi0_filename),Dict([("if_current",if_current),("current_strength",current_strength)]))

    time_ham = hamiltonian(metadata0["t1"],metadata0["t2"],metadata0["phi"],metadata0["U1"],metadata0["U2"],metadata0["L"],metadata0["nflavors"]; dict_to_symbols(metadata0)...,if_applied_current=if_current,current_strength=current_strength)

    naming_dict["dt"] = dt
	naming_dict["time_end"] = final_time
	filename = make_parameters_filename(naming_dict)
	println(filename)

    if_exists,data = check_data_exists("tevo-"*filename*".jld2","observer";location=location)
    if if_exists
        return data
    end

	tevo_params = (metadata=metadata,location=loc,name=filename)

	obs = evolve_in_time(psi0,final_time,dt,time_ham; tevo_params...,kwargs...)
    return obs
end










































"fin"
#####################################################
#=

This file contains the functions for time evolution of ED

Depends on:
    

=#
######################################################

# build Hamiltonian for given parameters and given time
function timeham(timestep::Int,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    
    # initialize hamilt_params
    localtime_hamilt_params::Dict = hamilt_params

    # reset hamilt_params given the timestep from the t_evo_params
    for (k,v) in t_evo_params
        if k != "dt" && k != "nsteps" && k != "tmax" && k != "when_dt_ends" && k != "current_dt"
            hamilt_params[k] = v[timestep]
        end
    end

    # build the Hamiltonian
    ht = buildHam(lattice_params,localtime_hamilt_params; kwargs...,output_level=0)

    return ht
end

function k1(wavefunc::Vector{ComplexF64},ht::SparseMatrixCSC,timestep::Int,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)

    # build the Hamiltonian for the given timestep
    #ht = timeham(timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)

    # calculate the k1 vector
    k1 = -im * (ht * wavefunc)

    opl > 2 && println("Finished k1")

    return k1
end

function k2(wavefunc::Vector{ComplexF64},k1::Vector{ComplexF64},timestep::Int,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)

    # build the Hamiltonian for the given timestep
    ht = timeham(timestep+1,t_evo_params,lattice_params,hamilt_params; kwargs...)

    # calculate the k2 vector
    k2 = -im * (ht * (wavefunc + ((0.5 * t_evo_params["current_dt"]) .* k1)))

    opl > 2 && println("Finished k2")

    return k2,ht
end

function k3(wavefunc::Vector{ComplexF64},k2::Vector{ComplexF64},ht::SparseMatrixCSC,timestep::Int,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)

    # build the Hamiltonian for the given timestep
    #ht = timeham(timestep+1,t_evo_params,lattice_params,hamilt_params; kwargs...)

    # calculate the k3 vector
    k3 = -im * (ht * (wavefunc + ((0.5 * t_evo_params["current_dt"]) .* k2)))

    opl > 2 && println("Finished k3")

    return k3
end

function k4(wavefunc::Vector{ComplexF64},k3::Vector{ComplexF64},timestep::Int,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)

    # build the Hamiltonian for the given timestep
    ht = timeham(timestep+2,t_evo_params,lattice_params,hamilt_params; kwargs...)

    # calculate the k4 vector
    k4 = -im * (ht * (wavefunc + (t_evo_params["current_dt"] .* k3)))

    opl > 2 && println("Finished k4")

    return k4,ht
end

function runge_kutta_step(wavefunc::Vector{ComplexF64},ht_prev::SparseMatrixCSC,timestep::Int,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)

    t_evo_params["current_dt"] = t_evo_params["dt"][1]
    timestep > t_evo_params["when_dt_ends"][1] && (t_evo_params["current_dt"] = t_evo_params["dt"][2])

    # calculate the k1, k2, k3 and k4 vectors
    k1_val = k1(wavefunc,ht_prev,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)
    k2_val,ht_half = k2(wavefunc,k1_val,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)
    k3_val = k3(wavefunc,k2_val,ht_half,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)
    k4_val,ht_next = k4(wavefunc,k3_val,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)

    # update the wavefunction
    new_wavefunc = wavefunc + ((t_evo_params["current_dt"] / 6) .* (k1_val + (2 .* k2_val) + (2 .* k3_val) + k4_val))

    opl > 2 && println("Finished Runge-Kutta step")

    return new_wavefunc,ht_next
end

function get_tevo_filename(timeevo_dict::Dict,lattice_dict::Dict,hamilt_dict::Dict; kwargs...)

    filename_dict = Dict{String,Any}()

    filename_dict["Lx"] = lattice_dict["Lx"]
    filename_dict["Ly"] = lattice_dict["Ly"]
    filename_dict["N"] = lattice_dict["N"]
    filename_dict["alpha"] = hamilt_dict["alpha"][2]
    filename_dict["if_periodic_x"] = lattice_dict["if_periodic_x"]
    filename_dict["if_periodic_y"] = lattice_dict["if_periodic_y"]

    filename_dict["dt"] = timeevo_dict["dt"]

    if hamilt_dict["disorder_strength"] != 0.0
        filename_dict["disorder_strength"] = hamilt_dict["disorder_strength"]
    end
    if haskey(hamilt_dict,"if_pinning") && hamilt_dict["if_pinning"]
        filename_dict["if_pinning"] = hamilt_dict["if_pinning"]
        filename_dict["pinning_strength"] = hamilt_dict["pinning_strength"]
    end
    if haskey(hamilt_dict,"periodic_potential_strength") && hamilt_dict["periodic_potential_strength"] != 0.0
        filename_dict["periodic_potential_strength"] = hamilt_dict["periodic_potential_strength"]
    end

    dataloc = get_folder_location("cluster-data/exact-diag/time-evo")
    if_both = 0

    for (k,v) in timeevo_dict
        if k != "dt" && k != "nsteps" && k != "tmax" && k != "when_dt_ends" && k != "other_dt" && k != "when_change_dt"
            filename_dict["rampparam"] = k
            if string(v[1]) == "linear_ramp"
                filename_dict["ramptype"] = "linear"
            elseif string(v[1]) == "pulse_ramp"
                filename_dict["ramptype"] = "pulse"
            else
                error("Unknown ramp type for time evolution parameter $k")
            end
        end

        # set data storage location based on the ramp parameter
        if k == "tx"
            dataloc = get_folder_location("cluster-data/exact-diag/time-evo/tx-ramp")
            if_both += 1
        elseif k == "intstren"
            dataloc = get_folder_location("cluster-data/exact-diag/time-evo/intstren-ramp")
            if_both += 1
        end
    end
    dataloc = get(kwargs, :dataloc, dataloc)

    # if ramping multiple parameters, save in "mixed-ramp" folder
    if if_both > 1
        dataloc = get_folder_location("cluster-data/exact-diag/time-evo/mixed-ramp")
    end

    if !haskey(timeevo_dict,"tx")
        filename_dict["hopping_anisotropy"] = hamilt_dict["tx"] / hamilt_dict["ty"]
    end

    # still need to figure out naming of interaction strength ramp
    if hamilt_dict["U"][2] == 0.0
        filename_dict["interaction_strength"] = 0.0
    else
        filename_dict["interaction_strength"] = hamilt_dict["U"][1]
    end
    if hamilt_dict["scaling_type"] != "flat"
		filename_dict["scaling"] = hamilt_dict["scaling_type"]
		if hamilt_dict["scaling_type"] == "gaussian"
			filename_dict["sigma"] = hamilt_dict["sigma"]
		elseif hamilt_dict["scaling_type"] == "exp"
			filename_dict["corr_length"] = hamilt_dict["corr_length"]
		elseif hamilt_dict["scaling_type"] == "rydberg"
			filename_dict["blockade_radius"] = hamilt_dict["blockade_radius"]
		else
			error("ULR Scaling Type Not Recognized: $(hamilt_dict["scaling_type"])")
		end
	end
    
    return dataloc,"tevo-" * make_parameters_filename(filename_dict) * ".jld2"
end

function save_tevo_data(tevo_wavefunc::SparseMatrixCSC,metadata::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)
    dataloc::String = metadata["dataloc"]
    filename::String = metadata["filename"]

    println("Saving time evolution data to $dataloc")

    # save metadata line
    data = Dict()
    write_data(filename,data,dataloc,metadata; kwargs...,dataloc=dataloc)
    
    # save wavefunction data
    for i in 1:size(tevo_wavefunc,2)
        data[string("tevowavefunc_",i)] = tevo_wavefunc[:,i]
    end
    return write_data("wavefunc"*filename,data,dataloc,metadata; kwargs...,dataloc=dataloc,output_level=opl-1)
end

function save_tevo_data(tevo_wavefunc::Vector,metadata::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)
    dataloc::String = metadata["dataloc"]
    filename::String = metadata["filename"]

    println("Saving time evolution data to $dataloc")

    # save metadata line
    data = Dict()
    write_data(filename,data,dataloc,metadata; kwargs...,dataloc=dataloc)
    
    # save wavefunction data
    for j in 1:length(tevo_wavefunc)
        for i in 1:size(tevo_wavefunc[1],2)
            data[string("tevowavefunc_gs$(j)_",i)] = tevo_wavefunc[j][:,i]
        end
    end
    return write_data("wavefunc"*filename,data,dataloc,metadata; kwargs...,dataloc=dataloc,output_level=opl-1)
end

# possibly introduce saving interval 
function save_tevo_data_local(local_wavefunc::Vector,timestep::Int; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)
    dataloc::String = kwargs[:dataloc]
    filename::String = kwargs[:filename]
    filepath = joinpath(dataloc,filename)

    # save wavefunction data
    data = Dict(string("tevowavefunc_",Int((timestep+1)/2)) => local_wavefunc)
    modify_data(data,filepath; kwargs...,output_level=opl-1)
end

function save_tevo_data_local(local_wavefunc::Vector{Vector},timestep::Int; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)
    dataloc::String = kwargs[:dataloc]
    filename::String = kwargs[:filename]
    filepath = joinpath(dataloc,filename)

    # save wavefunction data
    data = Dict()
    for i in 1:length(local_wavefunc)
        data[string("tevowavefunc_gs$(i)_",Int((timestep+1)/2))] = local_wavefunc[i]
    end
    modify_data(data,filepath; kwargs...,output_level=opl-1)

end

function time_evolution(starting_wavefunc::Vector{ComplexF64},starting_ham::SparseMatrixCSC,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)
    if_instant_gs::Bool = get(kwargs, :if_instant_gs, true)
    if_save_data::Bool = kwargs[:if_save_data]
    if_continuous_saving::Bool = kwargs[:if_continuous_saving]

    # make full metadata
    metadata = merge(lattice_params,hamilt_params,t_evo_params,named_tuple_to_dict(kwargs))

    # initialize the wavefunction
    wavefunc = starting_wavefunc

    # initialize ht_prev
    ht_prev = starting_ham

    nsteps::Int = t_evo_params["nsteps"]
    tevo_wavefunc = spzeros(ComplexF64,length(wavefunc),Int(1+(nsteps+1)/2))

    # if continuous saving, save the initial state
    if_save_data && if_continuous_saving && (actual_filename = save_tevo_data(tevo_wavefunc,metadata))

    if if_instant_gs
        nev = get(kwargs, :nev, 10)
        instant_spec = Dict{String,SparseMatrixCSC}()
        for i in 1:nev
            instant_spec[string(i)] = spzeros(ComplexF64,length(wavefunc),Int(1+(nsteps-1)/2))
        end
        running_args = get_quick_running_args(nev)
    end

    opl > 0 && println("Starting time evolution")

    display(t_evo_params)

    # perform the time evolution
    for timestep in 1:2:nsteps

        wavefunc_gs1,ht_new = runge_kutta_step(wavefunc,ht_prev,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)
        
        normalize!(wavefunc)

        if if_instant_gs
            
            fulloverlap = 0.0

            states,nrgs,rhos,hh = find_eigenstates(running_args.nev,lattice_params,hamilt_params; running_args...)
            for i in 1:nev
                local_state = states[i]
                fulloverlap += abs2(dot(wavefunc,local_state))
                instant_spec[string(i)][:,Int((timestep+1)/2)] = local_state
            end
            opl > 1 && println("Found instantaneous eigenstates at step $timestep")

            if fulloverlap < 1e-6
                error("State is Lost! Overlap with instantaneous eigenstates: $fulloverlap")
            end
        end
        
        tevo_wavefunc[:,Int((timestep+1)/2)] = wavefunc
        ht_prev = ht_new

        # save data if continuous saving is enabled
        if_save_data && if_continuous_saving && save_tevo_data_local(wavefunc,timestep; kwargs...,filename=actual_filename)
        
        opl > 0 && println("Finished $(round(timestep / nsteps * 100, digits = 2))% of time evolution")
    end

    # save data if not continuous saving
    if_save_data && !if_continuous_saving && save_tevo_data(tevo_wavefunc,metadata; kwargs...)

    opl > 0 && println("Time evolution completed.")

    if if_instant_gs
        return tevo_wavefunc,instant_spec
    else
        return tevo_wavefunc, nothing
    end
end

# run time evolution for multiple initial states (e.g. ground state and first excited state)
function time_evolution(starting_wavefunc::Vector{Vector{ComplexF64}},starting_ham::SparseMatrixCSC,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)
    if_instant_gs::Bool = get(kwargs, :if_instant_gs, true)
    if_save_data::Bool = kwargs[:if_save_data]
    if_continuous_saving::Bool = kwargs[:if_continuous_saving]

    # make full metadata
    metadata = merge(lattice_params,hamilt_params,t_evo_params,named_tuple_to_dict(kwargs))

    # initialize the wavefunction
    wavefunc::Vector{Vector{ComplexF64}} = Vector{Vector{ComplexF64}}(undef,length(starting_wavefunc))
    for i in 1:length(starting_wavefunc)
        wavefunc[i] = starting_wavefunc[i]
    end

    # initialize ht_prev
    ht_prev = starting_ham

    nsteps::Int = t_evo_params["nsteps"]
    tevo_wavefunc = [spzeros(ComplexF64,length(wavefunc[1]),Int(1+(nsteps+1)/2)) for i in 1:length(wavefunc)]

    # if continuous saving, save the initial state
    if_save_data && if_continuous_saving && (actual_filename = save_tevo_data(tevo_wavefunc,metadata))

    if if_instant_gs
        nev = get(kwargs, :nev, 10)
        instant_spec = Dict{String,SparseMatrixCSC}()
        for i in 1:nev
            instant_spec[string(i)] = spzeros(ComplexF64,length(wavefunc),Int(1+(nsteps-1)/2))
        end
        running_args = get_quick_running_args(nev)
    end

    opl > 0 && println("Starting time evolution")

    display(t_evo_params)

    # perform the time evolution
    for timestep in 1:2:nsteps

        for i in 1:length(wavefunc)
            wavefunc[i],ht_new = runge_kutta_step(wavefunc[i],ht_prev,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)
            ht_prev = ht_new
        end
        
        for i in 1:length(wavefunc)
            normalize!(wavefunc[i])
        end

        if if_instant_gs
            
            fulloverlap = 0.0

            states,nrgs,rhos,hh = find_eigenstates(running_args.nev,lattice_params,hamilt_params; running_args...)
            for i in 1:nev
                local_state = states[i]
                fulloverlap += abs2(dot(wavefunc,local_state))
                instant_spec[string(i)][:,Int((timestep+1)/2)] = local_state
            end
            opl > 1 && println("Found instantaneous eigenstates at step $timestep")

            if fulloverlap < 1e-6
                error("State is Lost! Overlap with instantaneous eigenstates: $fulloverlap")
            end
        end
        
        for i in 1:length(wavefunc)
            tevo_wavefunc[i][:,Int((timestep+1)/2)] = wavefunc[i]
        end

        # save data if continuous saving is enabled
        if_save_data && if_continuous_saving && save_tevo_data_local(wavefunc,timestep; kwargs...,filename=actual_filename)
        
        opl > 0 && println("Finished $(round(timestep / nsteps * 100, digits = 2))% of time evolution")
    end

    # save data if not continuous saving
    if_save_data && !if_continuous_saving && save_tevo_data(tevo_wavefunc,metadata; kwargs...)

    opl > 0 && println("Time evolution completed.")

    if if_instant_gs
        return tevo_wavefunc,instant_spec
    else
        return tevo_wavefunc, nothing
    end
end

function make_tevo_params(given_parameters::Dict)
    # initialize the time evolution parameters
    t_evo_params::Dict = Dict{String,Any}()

    t_evo_params["dt"] = [given_parameters["dt"], given_parameters["other_dt"]]
    t_evo_params["when_dt_ends"] = [given_parameters["when_change_dt"],given_parameters["nsteps"]]
    
    t_evo_params["nsteps"] = 2*given_parameters["nsteps"] - 1
    t_evo_params["tmax"] = t_evo_params["nsteps"] * t_evo_params["dt"]

    for (k,v) in given_parameters
        if k != "dt" && k != "nsteps" && k != "tmax" && k != "when_change_dt" && k != "other_dt"
            t_evo_params[k] = v[1](t_evo_params["nsteps"],t_evo_params["dt"][1]; v[2]...)
        end
    end

    return t_evo_params
end

function linear_ramp(nsteps::Int,dt::Float64; kwargs...)
    
    starting_value::Float64 = kwargs[:starting_value]
    ending_value::Float64 = kwargs[:ending_value]

    starting_time::Float64 = 0.0#get(kwargs, :starting_time, 0.0)
    ending_time::Float64 = get(kwargs, :ending_time, nsteps * dt)

    #steps_until_start::Int = Int(ceil(starting_time / dt))
    steps_until_end::Int = Int(ceil(ending_time / dt))

    return vcat(range(starting_value, ending_value, length = steps_until_end + 1), ending_value .* ones(nsteps - steps_until_end + 1))
end

function find_when_change_dt(tmax::Float64,leastramptime::Float64; kwargs...)
    max_nsteps::Int = get(kwargs, :max_nsteps, 1e4)

    dt = leastramptime / 3
    current_nsteps = Int(ceil(tmax / dt))

    if current_nsteps > max_nsteps
        steps_to_10x_ramptime = Int(ceil(leastramptime * 10 / dt))
        midtime = steps_to_10x_ramptime * dt
        remaining_steps = max_nsteps - steps_to_10x_ramptime
        other_dt = (tmax - midtime) / remaining_steps
    else
        other_dt = dt
        steps_to_10x_ramptime = 1
    end

    return steps_to_10x_ramptime, other_dt
end

function make_times(dts::Vector{Float64},when_dt_ends::Vector{Int})
    alltimes = zeros(Float64,when_dt_ends[end])
    for i in 2:when_dt_ends[1]
        alltimes[i] = alltimes[i-1] + dts[1]
    end
    for i in when_dt_ends[1]+1:when_dt_ends[2]
        alltimes[i] = alltimes[i-1] + dts[2]
    end
    return alltimes
end

function get_dt(tmax::Float64,leastramptime::Float64; kwargs...)
    
    default_dt = get(kwargs, :default_dt, 0.0005)
    default_nsteps = Int(ceil(tmax / default_dt))
    when_change_dt::Int = 1

    dt = default_dt

    other_dt::Float64 = dt

    if dt >= leastramptime
        #println("Least ramp time $leastramptime is larger than dt $dt, using leastramptime / 3 instead")
        dt = leastramptime / 3
        when_change_dt,other_dt = find_when_change_dt(tmax,leastramptime; kwargs...)
    end

    #println("Using time step dt = $dt for tmax = $tmax and least ramp time = $leastramptime")

    return dt,when_change_dt,other_dt
end

function get_maxramptime(time_params::Dict)
    all_ramptimes = []
    for (k,v) in time_params
        if k != "dt" && k != "tmax"
            push!(all_ramptimes,v[end])
        end
    end
    return maximum(all_ramptimes)
end

function get_leastramptime(time_params::Dict)
    all_ramptimes = []
    for (k,v) in time_params
        if k != "dt" && k != "tmax"
            push!(all_ramptimes,v[end])
        end
    end
    return minimum(all_ramptimes)
end

function run_timeevo(starting_gs::Vector,time_params::Dict,lattice_dict::Dict,hamilt_dict::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)
    
    tmax_global = 25.0
    dt_global = 0.05

    tmax::Float64 = time_params["tmax"]
    dt::Float64 = time_params["dt"]
    max_nsteps::Int = Int(ceil(tmax / dt))
    when_change_dt::Int = max_nsteps + 1
    
    tevo_pdict::Dict{String,Any} = Dict([("dt",dt),("tmax",tmax),("nsteps",max_nsteps),("when_change_dt",when_change_dt),("other_dt",dt)])

    # structure the control parameter values
    for (k,v) in time_params
        if k != "dt" && k != "tmax"
            if length(v) == 4
                tevo_pdict[k] = (v[1],(starting_value=v[2],ending_value=v[3],starting_time=0.0,ending_time=v[4]))
            elseif length(v) == 5    
                tevo_pdict[k] = (v[1],(starting_value=v[2],ending_value=v[3],starting_time=v[4],ending_time=v[5]))
            elseif length(v) == 3
                tevo_pdict[k] = (v[1],(ending_time=v[2],pulse_ramp=v[3]))
            else
                error("Invalid length of time evolution parameter $k: $(length(v))")
            end
        end
    end

    if_save_data::Bool = get(kwargs, :if_save_data, false)
    if_continuous_saving::Bool = get(kwargs, :if_continuous_saving, if_save_data && size(hamilt_dict["H"],1) > 10000)
    dataloc::String, filename::String = get_tevo_filename(tevo_pdict,lattice_dict,hamilt_dict; kwargs...)
    saving_args = (if_save_data=if_save_data,if_continuous_saving=if_continuous_saving,dataloc=dataloc,filename=filename,)

    if opl > 0
        println("Starting time evolution for $(lattice_dict["Lx"])x$(lattice_dict["Ly"]) N=$(lattice_dict["N"])")
        println("Saving filepath is $(joinpath(saving_args[:dataloc],saving_args[:filename]))")
        display(lattice_dict)
        display(hamilt_dict)
        display(tevo_pdict)
    end

    tevo_dict = make_tevo_params(tevo_pdict)
    
    tevo_groundstate,instantaneous_spectrum = time_evolution(starting_gs,hamilt_dict["H"],tevo_dict,lattice_dict,hamilt_dict; saving_args...,kwargs...) #output_level=1, nev=speccount

    return tevo_groundstate,tevo_dict,instantaneous_spectrum,saving_args
end






































"fin"
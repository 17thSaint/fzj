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
        if k != "dt" && k != "nsteps" && k != "tmax"
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
    k2 = -im * (ht * (wavefunc + ((0.5 * t_evo_params["dt"]) .* k1)))

    opl > 2 && println("Finished k2")

    return k2,ht
end

function k3(wavefunc::Vector{ComplexF64},k2::Vector{ComplexF64},ht::SparseMatrixCSC,timestep::Int,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)

    # build the Hamiltonian for the given timestep
    #ht = timeham(timestep+1,t_evo_params,lattice_params,hamilt_params; kwargs...)

    # calculate the k3 vector
    k3 = -im * (ht * (wavefunc + ((0.5 * t_evo_params["dt"]) .* k2)))

    opl > 2 && println("Finished k3")

    return k3
end

function k4(wavefunc::Vector{ComplexF64},k3::Vector{ComplexF64},timestep::Int,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)

    # build the Hamiltonian for the given timestep
    ht = timeham(timestep+2,t_evo_params,lattice_params,hamilt_params; kwargs...)

    # calculate the k4 vector
    k4 = -im * (ht * (wavefunc + (t_evo_params["dt"] .* k3)))

    opl > 2 && println("Finished k4")

    return k4,ht
end

function runge_kutta_step(wavefunc::Vector{ComplexF64},ht_prev::SparseMatrixCSC,timestep::Int,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)

    # calculate the k1, k2, k3 and k4 vectors
    k1_val = k1(wavefunc,ht_prev,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)
    k2_val,ht_half = k2(wavefunc,k1_val,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)
    k3_val = k3(wavefunc,k2_val,ht_half,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)
    k4_val,ht_next = k4(wavefunc,k3_val,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)

    # update the wavefunction
    new_wavefunc = wavefunc + ((t_evo_params["dt"] / 6) .* (k1_val + (2 .* k2_val) + (2 .* k3_val) + k4_val))

    opl > 2 && println("Finished Runge-Kutta step")

    return new_wavefunc,ht_next
end

function time_evolution(starting_wavefunc::Vector{ComplexF64},starting_ham::SparseMatrixCSC,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)
    if_instant_gs::Bool = get(kwargs, :if_instant_gs, true)

    # initialize the wavefunction
    wavefunc = starting_wavefunc

    # initialize ht_prev
    ht_prev = starting_ham

    nsteps::Int = t_evo_params["nsteps"]
    tevo_wavefunc = spzeros(ComplexF64,length(wavefunc),Int(1+(nsteps-1)/2))

    if if_instant_gs
        nev = get(kwargs, :nev, 10)
        instant_spec = Dict{String,SparseMatrixCSC}()
        for i in 1:nev
            instant_spec[string(i)] = spzeros(ComplexF64,length(wavefunc),Int(1+(nsteps-1)/2))
        end
        running_args = get_quick_running_args(nev)
    end

    opl > 0 && println("Starting time evolution")

    # perform the time evolution
    for timestep in 1:2:nsteps

        wavefunc,ht_new = runge_kutta_step(wavefunc,ht_prev,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)
        
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
        
        if opl > 0
            println("Finished $(round(timestep / nsteps * 100, digits = 2))% of time evolution")
        end
    end

    opl > 0 && println("Time evolution completed.")

    if if_instant_gs
        return tevo_wavefunc,instant_spec
    else
        return tevo_wavefunc
    end
end

function make_tevo_params(given_parameters::Dict)
    # initialize the time evolution parameters
    t_evo_params::Dict = Dict{String,Any}()

    # set the time step and number of steps
    t_evo_params["dt"] = get(given_parameters, "dt", 0.01)
    if haskey(given_parameters, "nsteps")
        t_evo_params["nsteps"] = 2*given_parameters["nsteps"] - 1
        t_evo_params["tmax"] = t_evo_params["nsteps"] * t_evo_params["dt"]
    else
        t_evo_params["nsteps"] = 2*Int(ceil(given_parameters["tmax"] / t_evo_params["dt"])) - 1
    end

    for (k,v) in given_parameters
        if k != "dt" && k != "nsteps" && k != "tmax"
            t_evo_params[k] = v[1](t_evo_params["nsteps"],t_evo_params["dt"]; v[2]...)
        end
    end

    return t_evo_params
end

function linear_ramp(nsteps::Int,dt::Float64; kwargs...)
    
    starting_value::Float64 = kwargs[:starting_value]
    ending_value::Float64 = kwargs[:ending_value]

    starting_time::Float64 = get(kwargs, :starting_time, 0.0)
    ending_time::Float64 = get(kwargs, :ending_time, nsteps * dt)

    steps_until_start::Int = Int(ceil(starting_time / dt))
    steps_until_end::Int = Int(ceil(ending_time / dt))

    return vcat(starting_value .* ones(steps_until_start), range(starting_value, ending_value, length = steps_until_end - steps_until_start + 1), ending_value .* ones(nsteps - steps_until_end + 1))
end

function get_dt(tmax::Float64; kwargs...)
    default_dt = get(kwargs, :default_dt, 0.0005)
    default_nsteps = Int(ceil(tmax / default_dt))
    if default_nsteps < 100
        dt = tmax / 100
    else
        dt = default_dt
    end
    return dt
end

function run_timeevo(starting_gs::Vector{ComplexF64},time_params::Dict,lattice_dict::Dict,hamilt_dict::Dict; kwargs...)
    
    max_ramp_time::Float64 = time_evo_dict["ramptime"]
    tmax::Float64 = 5*max_ramp_time
    
    dt::Float64 = get_dt(tmax; kwargs...)

    tevo_pdict = Dict([("dt",dt),("tmax",tmax),("tx",(linear_ramp,(starting_value=params_dict["tx"],ending_value=end_tx,ending_time=ramptime)))])

    # structure the control parameter values
    for (k,v) in time_params
        if k != "dt" && k != "tmax"
            if length(v) == 4
                tevo_pdict[k] = (v[1],(starting_value=v[2],ending_value=v[3],starting_time=0.0,ending_time=v[4]))
            elseif length(v) == 5    
                tevo_pdict[k] = (v[1],(starting_value=v[2],ending_value=v[3],starting_time=v[4],ending_time=v[5]))
            else
                error("Invalid length of time evolution parameter $k: $(length(v))")
            end
        end
    end

    tevo_dict = make_tevo_params(tevo_pdict)
    
    tevo_groundstate,instantaneous_spectrum = time_evolution(starting_gs,hamilt_dict["H"],tevo_dict,lattice_dict,hamilt_dict; kwargs...) #output_level=1, nev=speccount

    return tevo_groundstate,tevo_dict,instantaneous_spectrum
end






































"fin"
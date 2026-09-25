#####################################################
#=

This file contains the functions for time evolution of ED

Depends on:
    

=#
######################################################

# accumulated response int_0^t B(t') sin(w(t-t'))/w dt' of the magnetic field gradient up to
# the given timestep: the Green's function of the trapped mode xddot + w^2 x = B(t) started
# from rest, so a constant gradient rings around the static displacement B/w^2 instead of
# growing without bound. The trap frequency w is t_evo_params["trap_frequency"], which is
# always present: make_tevo_params sets it for every run, falling back to 1.0 (the kernel this
# was originally written with) only when the caller names no frequency.
# The ramp arrays are sampled on the RK4 half-step grid (spacing dt/2, see linear_ramp)
# and timeham indexes them directly, so the trapezoid sum runs over that spacing; the
# two-dt structure of t_evo_params["dt"] is respected in case the step size changes
# partway through the run. Zero at timestep 1, where the U profile is the undisplaced one.
function get_magnetic_gradient_integral(timestep::Int,t_evo_params::Dict)
    timestep < 2 && return 0.0

    gradient_values = t_evo_params["magnetic_gradient_time"]
    when_dt_ends::Int = t_evo_params["when_dt_ends"][1]
    dt_early::Float64 = t_evo_params["dt"][1] / 2
    dt_late::Float64 = t_evo_params["dt"][2] / 2
    trap_frequency::Float64 = t_evo_params["trap_frequency"]

    # spacings of the sampled points, then their times by prefix sum (t_1 = 0)
    spacings::Vector{Float64} = ifelse.((1:timestep-1) .> when_dt_ends, dt_late, dt_early)
    times::Vector{Float64} = pushfirst!(cumsum(spacings),0.0)
    t::Float64 = times[timestep]

    # retarded integrand B(t') * sin(w(t-t')) at every sample point, then the trapezoid over
    # it; the 1/w of the oscillator response divides the finished integral
    integrand::Vector{Float64} = @views gradient_values[1:timestep] .* sin.(trap_frequency .* (t .- times))

    return 0.5 * sum(@views spacings .* (integrand[1:timestep-1] .+ integrand[2:timestep])) / trap_frequency
end

# The same response at every sample at once, plus the mode's velocity:
#     x(t_k) = (1/w) int_0^t_k B(t') sin(w(t_k - t')) dt',   xdot(t_k) = int_0^t_k B(t') cos(w(t_k - t')) dt',
# for gradient samples separated by the interval lengths `steps` (t_1 = 0), with the trapezoid
# weights of get_magnetic_gradient_integral. Splitting sin(w(t-t')) = sin(wt)cos(wt') - cos(wt)sin(wt')
# pulls the upper limit out of the integrand and leaves two cumulative trapezoids, so all n samples
# cost O(n) together instead of O(n) each. Pass a Float64 h for a uniform grid.
function magnetic_gradient_response(gradient_values::AbstractVector{<:Real},steps::AbstractVector{<:Real},trap_frequency::Float64)
    n = length(gradient_values)
    @assert length(steps) == n - 1 "need one step per interval: $(n) samples, $(length(steps)) steps"
    wt = trap_frequency .* pushfirst!(cumsum(Float64.(steps)),0.0)
    c, s = cos.(wt), sin.(wt)
    cumtrap(f) = pushfirst!(cumsum(@views 0.5 .* steps .* (f[1:n-1] .+ f[2:n])),0.0)
    ic, is = cumtrap(gradient_values .* c), cumtrap(gradient_values .* s)
    return (s .* ic .- c .* is) ./ trap_frequency, c .* ic .+ s .* is
end

magnetic_gradient_response(gradient_values::AbstractVector{<:Real},h::Float64,trap_frequency::Float64) = magnetic_gradient_response(gradient_values,fill(h,length(gradient_values)-1),trap_frequency)

# Inverse of get_magnetic_gradient_integral: the gradient B(t) that drives a prescribed spacing
# a(t), read off the same trapped-mode equation, B = addot + w^2 (a - a0), with a0 the spacing
# the mode is displaced from (the dd spacing of the starting state, a[1] by default).
# spacing_values is sampled on a uniform grid of spacing h (the RK4 half-step dt/2 when it is a
# ramp array) and B is returned on the same samples: central second differences inside, and the
# second-order one-sided stencils at the two ends, so a smooth a(t) gives O(h^2) everywhere.
# Also returns the initial kick adot(0+). The mode starts from rest and a finite B only sets the
# acceleration, so a(t) leaving a0 with a slope cannot come from B alone -- it needs an impulse of
# that area at t = 0 (response kick*sin(wt)/w on top of the integral of B). Zero, to O(h^2), for
# any a(t) a finite gradient can produce; e.g. 1.5 for a linear ramp 0.5 -> 2.0 over T = 1.
# Where a(t) has a kink inside the grid, the returned B shows it as a single-sample spike of area
# ~ the velocity jump (the discrete delta); where B itself jumps, the sample straddling the
# switch gets a stencil-averaged value.
function get_magnetic_gradient_from_spacing(spacing_values::AbstractVector{<:Real},h::Float64,trap_frequency::Float64; a0::Float64=Float64(spacing_values[1]))
    n = length(spacing_values)
    @assert n >= 4 "need at least 4 samples for the one-sided end stencils, got $(n)"
    a = Float64.(spacing_values)

    addot = similar(a)
    @views addot[2:n-1] .= (a[3:n] .- 2 .* a[2:n-1] .+ a[1:n-2]) ./ h^2
    addot[1] = (2a[1] - 5a[2] + 4a[3] - a[4]) / h^2
    addot[n] = (2a[n] - 5a[n-1] + 4a[n-2] - a[n-3]) / h^2

    gradient = addot .+ trap_frequency^2 .* (a .- a0)
    kick = (-3a[1] + 4a[2] - a[3]) / (2h)

    return gradient, kick
end

# build Hamiltonian for given parameters and given time
function timeham(timestep::Int,t_evo_params::Dict,lattice_params::Dict,hamilt_params::Dict; kwargs...)
    
    # initialize hamilt_params
    localtime_hamilt_params::Dict = hamilt_params

    # reset hamilt_params given the timestep from the t_evo_params
    if_rebuild_ulr = false
    for (k,v) in t_evo_params
        if !(k in ("dt","nsteps","tmax","when_dt_ends","current_dt","trap_frequency","magnetic_gradient_integral_time"))
            # these parameters only enter the Hamiltonian through the coupling vector U,
            # which must be rebuilt whenever one of them takes a new value (skipping the
            # rebuild while the value sits constant, e.g. the hold after a ramp ends)
            if k in ("interaction_strength","corr_length","sigma","blockade_radius","magnetic_spacing","magnetic_gradient_time") && get(hamilt_params,k,nothing) !== v[timestep]
                if_rebuild_ulr = true
            end
            hamilt_params[k] = v[timestep]
        end
    end

    # the magnetic field gradient enters U through its time integral, which keeps growing
    # while the gradient is non-zero even when the gradient value itself repeats, so the
    # rebuild has to be keyed on the accumulated integral rather than on the ramped value
    # (and the integral, not the gradient, is what long_range_scaling needs). make_tevo_params
    # precomputes it for the whole run, so this is a lookup rather than an O(n) sum per call
    magnetic_integral::Float64 = 0.0
    if haskey(t_evo_params,"magnetic_gradient_time")
        magnetic_integral = t_evo_params["magnetic_gradient_integral_time"][timestep]
        magnetic_integral !== get(hamilt_params,"magnetic_gradient_integral",nothing) && (if_rebuild_ulr = true)
        hamilt_params["magnetic_gradient_integral"] = magnetic_integral
    end

    if if_rebuild_ulr
        # applyHam reads the interaction from hamilt_params["U"] only, so the ramped value
        # has to be propagated into U or the interaction stays frozen at its starting value;
        # the fallbacks cover hamilt_params dicts built before "lr_dist" and
        # "interaction_strength" were stored at setup (U[1] equals the overall strength
        # for every scaling type in long_range_scaling)
        lr_dist::Int = get(hamilt_params,"lr_dist",length(hamilt_params["U"])-1)
        stren::Float64 = get(hamilt_params,"interaction_strength",hamilt_params["U"][1])
        hamilt_params["U"] = long_range_scaling(lr_dist,lattice_params["Ly"],stren;
            scaling=hamilt_params["scaling_type"],corr_length=hamilt_params["corr_length"],
            sigma=hamilt_params["sigma"],blockade_radius=hamilt_params["blockade_radius"],
            magnetic_spacing=hamilt_params["magnetic_spacing"],magnetic_gradient_integral=magnetic_integral)
    end

    # build the Hamiltonian, from the saved undressed matrices when if_reading is set
    if get(kwargs,:if_reading,false)
        ht = getHamiltonian(lattice_params,localtime_hamilt_params; kwargs...,output_level=0)
    else
        ht = buildHam(lattice_params,localtime_hamilt_params; kwargs...,output_level=0)
    end

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

    #filename_dict["dt"] = timeevo_dict["dt"]

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
        if k != "dt" && k != "nsteps" && k != "tmax" && k != "when_dt_ends" && k != "other_dt" && k != "when_change_dt" && k != "trap_frequency"
            filename_dict["rampparam"] = k
            if string(v[1]) == "linear_ramp"
                filename_dict["ramptype"] = "linear"
            elseif string(v[1]) == "pulse_ramp"
                filename_dict["ramptype"] = "pulse"
            else
                error("Unknown ramp type for time evolution parameter $k")
            end
        end
    end
    dataloc = get(kwargs, :dataloc, dataloc)

    # two gradient runs that differ only in the trap frequency agree on every other field
    # here, so w has to enter the name or their saved data collides; the default 1.0 is left
    # out to keep the names of runs made before the trap frequency was a parameter
    if timeevo_dict["trap_frequency"] != 1.0
        filename_dict["trap_frequency"] = timeevo_dict["trap_frequency"]
    end

    filename_dict["ramptime"] = timeevo_dict["tmax"]

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
		elseif hamilt_dict["scaling_type"] == "dd"
			filename_dict["magnetic_spacing"] = hamilt_dict["magnetic_spacing"]
		elseif hamilt_dict["scaling_type"] == "magnetic_gradient"
			# magnetic_spacing is the initial spacing a0 the gradient displaces from
			filename_dict["magnetic_spacing"] = hamilt_dict["magnetic_spacing"]
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

    opl > 0 && display(t_evo_params)

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
    if_instant_gs::Bool = get(kwargs, :if_instant_gs, false)
    if_instant_exact::Bool = get(kwargs, :if_instant_exact, if_instant_gs)
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
    tevo_nrg = [zeros(Float64,Int(1+(nsteps+1)/2)) for i in 1:length(wavefunc)]

    # if continuous saving, save the initial state
    if_save_data && if_continuous_saving && (actual_filename = save_tevo_data(tevo_wavefunc,metadata))

    if if_instant_gs
        nev = get(kwargs, :nev, 10)
        instant_spec = Dict{String,SparseMatrixCSC}()
        instant_nrgs = Dict{String,Vector{Float64}}()
        for i in 1:nev
            instant_spec[string(i)] = spzeros(ComplexF64,length(wavefunc[1]),Int(1+(nsteps-1)/2))
            instant_nrgs[string(i)] = zeros(Float64,Int(1+(nsteps-1)/2))
        end
        # exact (dense) diagonalization correctly resolves degenerate eigenvalues that a
        # single-vector Lanczos run structurally cannot (a Krylov chain from one starting
        # vector only ever yields one Ritz pair per degenerate eigenvalue, however large
        # krylovdim is) -- only worth it for small Hilbert spaces since it rebuilds and
        # fully diagonalizes the dense Hamiltonian every RK4 step
        running_args = get_quick_running_args(nev; if_exact=if_instant_exact)
    end

    opl > 0 && println("Starting time evolution")

    opl > 0 && display(t_evo_params)

    # perform the time evolution
    for timestep in 1:2:nsteps

        ht_start = ht_prev
        for i in 1:length(wavefunc)
            wavefunc[i],ht_new = runge_kutta_step(wavefunc[i],ht_start,timestep,t_evo_params,lattice_params,hamilt_params; kwargs...)
            ht_prev = ht_new
        end
        
        for i in 1:length(wavefunc)
            normalize!(wavefunc[i])
        end

        if if_instant_gs
            
            #fulloverlap = 0.0

            states,nrgs,rhos,hh = find_eigenstates(running_args.nev,lattice_params,hamilt_params; running_args...)
            for i in 1:nev
                local_state = states[i]
                #fulloverlap += abs2(adjoint(wavefunc[1]) * local_state)
                instant_spec[string(i)][:,Int((timestep+1)/2)] = local_state
                instant_nrgs[string(i)][Int((timestep+1)/2)] = nrgs[i]
            end
            opl > 1 && println("Found instantaneous eigenstates at step $timestep")

            #fulloverlap < 1e-6 && error("State is Lost! Overlap with instantaneous eigenstates: $fulloverlap")
        end
        
        for i in 1:length(wavefunc)
            tevo_wavefunc[i][:,Int((timestep+1)/2)] = wavefunc[i]
            tevo_nrg[i][Int((timestep+1)/2)] = real(adjoint(wavefunc[i]) * (ht_prev * wavefunc[i]))
        end

        # save data if continuous saving is enabled
        if_save_data && if_continuous_saving && save_tevo_data_local(wavefunc,timestep; kwargs...,filename=actual_filename)
        
        opl > 0 && println("Finished $(round(timestep / nsteps * 100, digits = 2))% of time evolution")
    end

    # save data if not continuous saving
    if_save_data && !if_continuous_saving && save_tevo_data(tevo_wavefunc,metadata; kwargs...)

    opl > 0 && println("Time evolution completed.")

    if if_instant_gs
        return [tevo_wavefunc,tevo_nrg],[instant_spec,instant_nrgs]
    else
        return [tevo_wavefunc,tevo_nrg], [nothing, nothing]
    end
end

function make_tevo_params(given_parameters::Dict)
    # initialize the time evolution parameters
    t_evo_params::Dict = Dict{String,Any}()

    t_evo_params["dt"] = [given_parameters["dt"], given_parameters["other_dt"]]
    t_evo_params["when_dt_ends"] = [given_parameters["when_change_dt"],given_parameters["nsteps"]]
    
    t_evo_params["nsteps"] = 2*given_parameters["nsteps"] - 1
    t_evo_params["tmax"] = t_evo_params["nsteps"] * t_evo_params["dt"]

    # the trap frequency is a fixed scalar of the gradient kernel rather than a ramped
    # control, so it is carried over as-is and left out of the ramp expansion below. Set on
    # every run, so get_magnetic_gradient_integral can read it straight out of t_evo_params;
    # callers that name no frequency get the 1.0 the retarded kernel was first written with
    t_evo_params["trap_frequency"] = get(given_parameters,"trap_frequency",1.0)

    for (k,v) in given_parameters
        if k != "dt" && k != "nsteps" && k != "tmax" && k != "when_change_dt" && k != "other_dt" && k != "trap_frequency"
            t_evo_params[k] = v[1](t_evo_params["nsteps"],t_evo_params["dt"][1]; v[2]...)
        end
    end

    # the retarded gradient integral timeham needs at every half-step, computed once for the
    # whole run, with the sample spacings get_magnetic_gradient_integral uses (two-dt aware)
    if haskey(t_evo_params,"magnetic_gradient_time")
        n = length(t_evo_params["magnetic_gradient_time"])
        dt_early, dt_late = t_evo_params["dt"] ./ 2
        steps = ifelse.((1:n-1) .> t_evo_params["when_dt_ends"][1], dt_late, dt_early)
        t_evo_params["magnetic_gradient_integral_time"] = magnetic_gradient_response(t_evo_params["magnetic_gradient_time"],steps,t_evo_params["trap_frequency"])[1]
    end

    return t_evo_params
end

function linear_ramp(nsteps::Int,dt::Float64; kwargs...)

    starting_value::Float64 = kwargs[:starting_value]
    ending_value::Float64 = kwargs[:ending_value]

    starting_time::Float64 = 0.0#get(kwargs, :starting_time, 0.0)

    # nsteps counts raw RK4 half-steps (timeham is indexed directly by this array),
    # each spaced dt/2 apart in real time, so ending_time must be converted using
    # the half-step spacing rather than the full step dt
    raw_dt::Float64 = dt / 2
    ending_time::Float64 = get(kwargs, :ending_time, nsteps * raw_dt)

    #steps_until_start::Int = Int(ceil(starting_time / raw_dt))
    steps_until_end::Int = Int(ceil(ending_time / raw_dt))

    return vcat(range(starting_value, ending_value, length = steps_until_end + 1), ending_value .* ones(nsteps - steps_until_end + 1))
end

# Applies an externally supplied pulse (e.g. from a QuOCS optimization) as the shape
# of a ramp control, on the same half-step grid as linear_ramp. Must be named
# pulse_ramp: get_tevo_filename below whitelists ramp functions by literal name.
# pulse_ramp's required length is ceil(ending_time / (dt/2)) + 1.
function pulse_ramp(nsteps::Int,dt::Float64; kwargs...)

    ending_time::Float64 = kwargs[:ending_time]
    pulse_values = kwargs[:pulse_ramp][:]

    raw_dt::Float64 = dt / 2
    steps_until_end::Int = Int(ceil(ending_time / raw_dt))

    @assert length(pulse_values) == steps_until_end + 1 "Pulse length $(length(pulse_values)) does not match the expected number of ramp half-steps $(steps_until_end + 1) for ending_time=$(ending_time), dt=$(dt)"

    return vcat(pulse_values, pulse_values[end] .* ones(nsteps - steps_until_end + 1))
end

# the sample times of that grid for a pulse of length ending_time: the ceil(ending_time/(dt/2)) + 1
# half-steps pulse_ramp expects
halfstep_times(ending_time::Float64,dt::Float64) = (0:Int(ceil(ending_time/(dt/2)))) .* (dt/2)

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

# Largest interaction that actually appears in the Hamiltonian. U is indexed by synthetic
# separation as U[dist+1], and applyHam only pairs distinct particles that already share
# the other coordinate, so dist = 0 would mean two particles on one site -- unreachable in
# a hardcore basis. U[1] therefore never enters the energy (it only feeds the lr_dist reach
# count) and must be left out of any estimate of the spectral scale. abs is needed because
# the tail alone can be attractive, where a bare maximum would return a negative scale.
function get_interaction_scale(hamilt_dict::Dict)
    us = hamilt_dict["U"]
    return length(us) > 1 ? maximum(abs,@view us[2:end]) : 0.0
end

# Largest RK4 step that stays inside the stability region for the current interaction
# scale: the spectral radius is estimated as (number of interacting pairs) * the largest
# coupling that enters H, and RK4 on the imaginary axis is stable up to |E|*dt ~ 2.
# Exposed so callers that have to build a control pulse on the half-step grid before
# starting the run (the grid spacing is dt/2) can size that grid with the same step the
# evolution will use.
# The second branch keeps the step from being so coarse that a short tmax gets only a
# handful of samples, by falling back to the step a U = 300 interaction would demand. It
# also catches a vanishing interaction scale, where the first expression is Inf: the
# hopping, not U, sets the limit there and this fallback stands in for it.
function get_critical_dt(tmax::Float64,lattice_dict::Dict,hamilt_dict::Dict)
    npairs = lattice_dict["N"] * (lattice_dict["N"]-1) / 2
    dt_crit = 2.0 / (npairs * get_interaction_scale(hamilt_dict))
    dt_crit > 0.01*tmax && (dt_crit = 2.0 / (npairs * 300.0))
    return dt_crit
end

# The same limit for the dd profile at an arbitrary spacing, for spacing / gradient pulses: the
# step from the t = 0 Hamiltonian only bounds a run while a(t) stays at or above its starting
# value, and U ~ 1/a^3 grows quickly as a pulse pulls the states together
function get_critical_dt_at_spacing(spacing::Float64,tmax::Float64,lattice_dict::Dict,hamilt_dict::Dict)
    hp = copy(hamilt_dict)
    hp["U"] = long_range_scaling(get(hamilt_dict,"lr_dist",length(hamilt_dict["U"])-1),lattice_dict["Ly"],
                                 hamilt_dict["interaction_strength"]; scaling="dd",magnetic_spacing=spacing)
    return get_critical_dt(tmax,lattice_dict,hp)
end

function run_timeevo(starting_gs::Vector,time_params::Dict,lattice_dict::Dict,hamilt_dict::Dict; kwargs...)
    opl::Int = get(kwargs, :output_level, 1)
    
    #tmax_global = 25.0
    #dt_global = 0.05

    tmax = time_params["tmax"]    

    dt_crit = get_critical_dt(tmax,lattice_dict,hamilt_dict)

    # use caller-supplied dt if provided, otherwise fall back to interaction-derived critical
    # step. A supplied step above dt_crit is outside the RK4 stability region: the top of the
    # spectrum is amplified every step and, because time_evolution renormalizes the state each
    # step, the divergence never shows as a growing norm -- it silently turns the transported
    # state into the highest-energy eigenvector instead, so this has to be rejected rather than
    # warned about
    dt::Float64 = haskey(time_params, "dt") ? time_params["dt"] : dt_crit
    if dt > dt_crit
        error("Time step dt = $dt is above the RK4 stability limit dt_crit = $dt_crit for " *
              "N = $(lattice_dict["N"]) and max interacting U = $(get_interaction_scale(hamilt_dict)): " *
              "the evolution would diverge and be hidden by the per-step renormalization. " *
              "Pass a smaller dt, " *
              "omit \"dt\" from time_params to use dt_crit, or lower the interaction scale.")
    end
    max_nsteps::Int = Int(ceil(tmax / dt))
    when_change_dt::Int = max_nsteps + 1
    
    tevo_pdict::Dict{String,Any} = Dict([("dt",dt),("tmax",tmax),("nsteps",max_nsteps),("when_change_dt",when_change_dt),("other_dt",dt)])

    # the trap frequency of the magnetic-gradient kernel is a plain scalar, not a ramp
    # triple/quadruple, so it bypasses the control structuring and reaches make_tevo_params
    # (and from there get_magnetic_gradient_integral) unchanged. Always set, so every dict
    # downstream of here carries it whether or not the caller named one: an explicit entry in
    # time_params wins, otherwise the frequency the state was set up with in
    # get_normal_model_params_ed, otherwise the 1.0 the retarded kernel was first written with
    tevo_pdict["trap_frequency"] = get(time_params,"trap_frequency",get(hamilt_dict,"trap_frequency",1.0))

    # structure the control parameter values
    for (k,v) in time_params
        if k != "dt" && k != "tmax" && k != "trap_frequency"
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

    # build the t=0 Hamiltonian from the ramp's actual starting values rather than reusing
    # hamilt_dict["H"], which may still carry a one-off pinning/barrier term baked in by whatever
    # produced the starting state (e.g. position_state's corner-confinement potential)
    starting_ham = timeham(1,tevo_dict,lattice_dict,hamilt_dict; kwargs...)

    tevo_data,instant_data = time_evolution(starting_gs,starting_ham,tevo_dict,lattice_dict,hamilt_dict; saving_args...,kwargs...) #output_level=1, nev=speccount

    return tevo_data,tevo_dict,instant_data,saving_args
end

# Runs a chain of sequential ramp stages, feeding each stage's final states in as the
# next stage's starting states. Each stage is (control_key,pulse_values,duration): the
# Hamiltonian parameter named control_key is ramped via pulse_ramp over [0,duration]
# while every other parameter stays fixed at whatever hamilt_params currently holds.
# Used by control-functions.jl's figure-of-merit functions to turn one or more
# QuOCS-optimized pulses into a final state for fidelity comparison.
function run_ramp_stages(states::Vector{Vector{ComplexF64}},stages,lattice_params::Dict,hamilt_params::Dict,dt::Float64; kwargs...)
    # timeham writes each ramped parameter's current value back into hamilt_params as a
    # side effect; work on a copy so stages still chain correctly (each stage starts from
    # the previous stage's final control values) without corrupting the caller's dict
    hamilt_params = copy(hamilt_params)
    for (control_key,pulse_values,duration) in stages
        tevo_params = Dict([(control_key,(pulse_ramp,duration,pulse_values)),("tmax",duration),("dt",dt)])
        tevo_data,_,_,_ = run_timeevo(states,tevo_params,lattice_params,hamilt_params; kwargs...)
        # end-1 skips the final save point which lands at tmax rather than the last full Trotter step
        states = [Vector{ComplexF64}(tevo_data[1][i][:,end-1]) for i in 1:length(states)]
    end
    return states
end






































"fin"
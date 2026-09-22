#####################################################
#=

Setup for the QuOCS magnetic-gradient-pulse optimization (optimal-control/config_maggradPulse.py).

The figure of merit there is purely classical -- it only needs the driven response
a(t) = a0 + int_0^t B(t') sin(w(t-t')) dt' / w of the displaced trap mode -- so it is
computed in Python and this file is not called per evaluation. The one thing the Python
side cannot work out for itself is the time step: the pulse has to be sampled on the same
RK4 half-step grid the eventual run_timeevo uses, and that step follows from the
interaction scale of the t = 0 Hamiltonian, which needs the ED. Run once at setup.

Depends on:
    execute-ed.jl, time-evolution.jl

=#
######################################################


# The t = 0 Hamiltonian of the gradient run: the accumulated integral is zero there, so the
# "magnetic_gradient" scaling reduces to the ordinary "dd" profile at the initial spacing a0
# and this is also the ED whose groundstate the pulse starts from.
function setup_maggrad_pulse(params_dict::Dict)
    pdict = Dict{String,Any}(
        "output_level" => 0,
        "Lx" => params_dict["Lx"],
        "Ly" => params_dict["Ly"],
        "N" => params_dict["N"],
        "lr" => params_dict["lr"],
        "if_periodic_x" => params_dict["if_periodic_x"],
        "if_periodic_y" => params_dict["if_periodic_y"],
        "hopping_anisotropy" => 1.0,
        "scaling_type" => "dd",
        "trap_frequency" => params_dict["trap_frequency"],
        "magnetic_spacing" => params_dict["a0"],
        "interaction_strength" => params_dict["interaction_strength"],
        "filling" => 0.5,
        "nev" => params_dict["speccount"],
        "if_find_data" => false,
        "if_save_data" => false)

    _,_,_,_,_,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=0)

    return (lattice_params=lattice_params, hamilt_params=hamilt_params,
            dt=get_critical_dt(params_dict["tmax"],lattice_params,hamilt_params))
end


# The RK4 stability limit at an arbitrary spacing, for checking a finished pulse. The step
# setup_maggrad_pulse returns is the one for a0, which is the tightest spacing only while
# a(t) >= a0; B is free to go negative here, so a pulse that pulls the states together
# raises max(U) ~ 1/a^3 and demands a smaller step than the grid was built with.
function maggrad_critical_dt_at_spacing(spacing::Float64,tmax::Float64,
                                        lattice_params::Dict,hamilt_params::Dict)
    hp = copy(hamilt_params)
    hp["U"] = long_range_scaling(hamilt_params["lr_dist"],lattice_params["Ly"],
                                 hamilt_params["interaction_strength"];
                                 scaling="dd",magnetic_spacing=spacing)
    return get_critical_dt(tmax,lattice_params,hp)
end

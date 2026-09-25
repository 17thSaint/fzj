#####################################################
#=

Setup for the QuOCS magnetic-gradient-pulse optimization (optimal-control/config_maggradPulse.py).

The figure of merit there is purely classical -- it only needs the driven response
a(t) = a0 + int_0^t B(t') sin(w(t-t')) dt' / w of the displaced trap mode -- so it is
computed in Python and this file is not called per evaluation. The one thing the Python
side cannot work out for itself is the time step: the pulse has to be sampled on the same
RK4 half-step grid the eventual run_timeevo uses, and that step follows from the
interaction scale of the t = 0 Hamiltonian, which needs the ED. Run once at setup.
A finished pulse that dips below a0 is checked with get_critical_dt_at_spacing
(time-evolution.jl).

Depends on:
    control-functions.jl (optimization_ed_params), time-evolution.jl

=#
######################################################


# The t = 0 Hamiltonian of the gradient run: the accumulated integral is zero there, so the
# "magnetic_gradient" scaling reduces to the ordinary "dd" profile at the initial spacing a0
# and this is also the ED whose groundstate the pulse starts from.
function setup_maggrad_pulse(params_dict)
    pdict = optimization_ed_params(params_dict,"scaling_type"=>"dd","magnetic_spacing"=>params_dict["a0"])
    _,_,_,_,_,lattice_params,hamilt_params = run_normal_ed(pdict; output_level=0)

    return (lattice_params=lattice_params, hamilt_params=hamilt_params,
            dt=get_critical_dt(params_dict["tmax"],lattice_params,hamilt_params))
end

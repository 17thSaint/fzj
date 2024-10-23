#####################################################
#=

This file contains observer types for 1D effective MPS results

Depends on:

=#
######################################################

mutable struct NRGObserver <: AbstractObserver
    nrg_tol::Float64
    local_ham
    nrg::Vector{Float64}
 
    NRGObserver(nrg_tol=0.0,local_ham=10.0) = new(nrg_tol,local_ham,[1000.0])
end

function ITensors.checkdone!(o::NRGObserver;kwargs...)
	if abs(o.nrg[end] - o.nrg[end-1]) < o.nrg_tol && length(o.nrg) > 10
	  #println("Stopping DMRG after sweep $sw")
	  return true
	end
	return false
end
  
function ITensors.measure!(o::NRGObserver; kwargs...)
	  #display(kwargs)
	  half_sweep = kwargs[:half_sweep]
	  bond = kwargs[:bond]
	  outputlevel = kwargs[:outputlevel]

	  psi = kwargs[:psi]
	  ham = o.local_ham
	
	  if bond == 1 && half_sweep == 2
		# Update last_energy and keep going
		append!(o.nrg,[calculate_energy(psi,ham)])
		if outputlevel > 0
			println("The energy change is $(round(abs(o.nrg[end] - o.nrg[end-1]),digits=7)) for tolerance $(o.nrg_tol)")
		end
	  end
end

mutable struct LinkDimObserver <: AbstractObserver
    mlds::Vector{Int64}
 
    LinkDimObserver() = new([10000])
end

function ITensors.checkdone!(o::LinkDimObserver;kwargs...)
	if length(o.mlds) > 10
		if all(o.mlds[end-5:end] .== o.mlds[end])
			return true
		end
	end
	return false
end
  
function ITensors.measure!(o::LinkDimObserver; kwargs...)
	  #display(kwargs)
	  half_sweep = kwargs[:half_sweep]
	  bond = kwargs[:bond]
	  outputlevel = kwargs[:outputlevel]

	  psi = kwargs[:psi]
	
	  if bond == 1 && half_sweep == 2
		# Update max link dim
		append!(o.mlds,[maxlinkdim(psi)])
		if outputlevel > 0 && length(o.mlds) > 5
			println("Current Percent Change = ",round(100*std(o.mlds[end-5:end])/mean(o.mlds[end-5:end]),digits=2))
		end
	  end
end

mutable struct NRGErrorObserver <: AbstractObserver

	var_tol::Float64
	local_ham
	nrgs::Vector{Float64}
    nrg_var::Vector{Float64}
 
    NRGErrorObserver(var_tol=0.0,local_ham=10.0) = new(var_tol,local_ham,[10000.0,1000.0],[0.0])
end

function ITensors.checkdone!(o::NRGErrorObserver;kwargs...)
	sw = kwargs[:sweep]
	#psi = kwargs[:psi]
	#ham = o.local_ham
	if o.nrg_var[end] < o.var_tol && abs(o.nrgs[end] - o.nrgs[end-1]) < o.nrg_var[end]
	  return true
	end
	#o.nrgs[1] = o.nrgs[2]
	return false
end
  
function ITensors.measure!(o::NRGErrorObserver; kwargs...)
	  half_sweep = kwargs[:half_sweep]
	  bond = kwargs[:bond]
	  outputlevel = kwargs[:outputlevel]
	
	  if bond == 1 && half_sweep == 2 && outputlevel > 0
		psi = kwargs[:psi]
	    ham = o.local_ham
	    append!(o.nrg_var,[energy_variance(psi,ham)])
	    append!(o.nrgs,[real(calculate_energy(psi,ham))])
		if outputlevel > 0
			println("The energy variance is $(round(o.nrg_var[end],digits=10)) with energy change $(round(abs(o.nrgs[end] - o.nrgs[end-1]),digits=10))")
		end
	end
end
































"fin"
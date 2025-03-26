#####################################################
#=

This file is for reading/writing data types in HDF5 format.

=#
######################################################

using HDF5

#include("ttn.jl")

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, observer::Type{<:TTN.ITensorMPS.AbstractObserver})
    group = open_group(parent, name)
    
    observer_type = read(attributes(group)["type"])

    if observer_type == "SavingNRGVarObserver"
        return read_SNVO(group, SavingNRGVarObserver)
    elseif observer_type == "SavingMeasurementsObserver"
        return read_SMO(group, SavingMeasurementsObserver)
    elseif observer_type == "SavingExcitedNRGVarObserver"
        return read_SENVO(group, SavingExcitedNRGVarObserver)
    elseif observer_type == "NRGVarObserver"
        return read_NVO(group, NRGVarObserver)
    else
        error("Observer type $observer_type not recognized")
    end
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, observer::SavingNRGVarObserver)
	group = create_group(parent, name)
	
	file_path = observer.file_path
	var_tol = observer.var_tol
	nrg = observer.nrg

	write(group, "file_path", file_path)
	write(group, "var_tol", var_tol)
	write(group, "nrg", nrg)

    attributes(group)["type"] = "SavingNRGVarObserver"
end

function read_SNVO(group::HDF5.Group, observer::Type{<:SavingNRGVarObserver})
	#group = open_group(parent, name)
	
	file_path = read(group, "file_path")
	var_tol = read(group, "var_tol")
	nrg = read(group, "nrg")
	
	return SavingNRGVarObserver(file_path, var_tol, nrg)
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, observer::NRGVarObserver)
	group = create_group(parent, name)
	
	var_tol = observer.var_tol
	nrg = observer.nrg

	write(group, "var_tol", var_tol)
	write(group, "nrg", nrg)

    attributes(group)["type"] = "NRGVarObserver"
end

function read_NVO(group::HDF5.Group, observer::Type{<:NRGVarObserver})
	#group = open_group(parent, name)
	
	var_tol = read(group, "var_tol")
	nrg = read(group, "nrg")
	
	return NRGVarObserver(var_tol, nrg)
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, observer::SavingMeasurementsObserver)
	group = create_group(parent, name)
	
	file_path = observer.file_path
	var_tol = observer.var_tol
	nrg = observer.nrg
	measurement_functions = observer.measurement_functions
	measurements = observer.measurements

	write(group, "file_path", file_path)
	write(group, "var_tol", var_tol)
	write(group, "nrg", nrg)
	write(group, "measurement_functions", measurement_functions)
	write(group, "measurements", measurements)

    attributes(g)["type"] = "SavingMeasurementsObserver"
end

function read_SMO(group::HDF5.Group, observer::Type{<:SavingMeasurementsObserver})
	#group = open_group(parent, name)
	
	file_path = read(group, "file_path")
	var_tol = read(group, "var_tol")
	nrg = read(group, "nrg")
	measurement_functions = read(group, "measurement_functions")
	measurements = read(group, "measurements")
	
	return SavingMeasurementsObserver(measurement_functions, measurements, file_path, var_tol, nrg)
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, observer::SavingExcitedNRGVarObserver)
    group = create_group(parent, name)
    
    file_path = observer.file_path
    var_tol = observer.var_tol
    nrg_level = observer.nrg_level
    nrg = observer.nrg

    write(group, "file_path", file_path)
    write(group, "var_tol", var_tol)
    write(group, "nrg_level", nrg_level)
    write(group, "nrg", nrg)

    attributes(group)["type"] = "SavingExcitedNRGVarObserver"
end

function read_SENVO(group::HDF5.Group, observer::Type{<:SavingExcitedNRGVarObserver})
    #group = open_group(parent, name)
    
    file_path = read(group, "file_path")
    var_tol = read(group, "var_tol")
    nrg_level = read(group, "nrg_level")
    nrg = read(group, "nrg")
    
    return SavingExcitedNRGVarObserver(file_path, var_tol, nrg_level, nrg)
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, ham_op::TTN.ITensorMPS.Sum{TTN.ITensorMPS.Scaled{ComplexF64, TTN.ITensorMPS.Prod{TTN.ITensorMPS.Op}}})
    group = create_group(parent, name)

    for (i, val) in enumerate(ham_op)
        local_subgroup = create_group(group, "op_$i")

        coeff = coefficient(val)
        for (j, t) in enumerate(TTN.terms(val))

            # find each local operator and sites
            local_op = TTN.which_op(t)
            local_sites = TTN.site(t)

            # write data to local_subgroup
            write(local_subgroup, "op_$j", local_op)
            write(local_subgroup, "sites_$j", collect(local_sites))

        end

        write(local_subgroup, "coeff", coeff)
    end

end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, ham::Type{<:TTN.ITensorMPS.Sum{TTN.ITensorMPS.Scaled{ComplexF64, TTN.ITensorMPS.Prod{TTN.ITensorMPS.Op}}}})
    group = open_group(parent, name)

    ham_op = TTN.OpSum()

    for i in 1:length(group)
        local_subgroup = open_group(group, "op_$i")

        coeff = read(local_subgroup, "coeff")
        
        if length(local_subgroup) == 5
            op1 = read(local_subgroup, "op_1")
            op2 = read(local_subgroup, "op_2")

            sites1 = read(local_subgroup, "sites_1")
            sites2 = read(local_subgroup, "sites_2")

            ham_op += (coeff,op1,sites1,op2,sites2)
        elseif length(local_subgroup) == 3
            op1 = read(local_subgroup, "op_1")
            sites1 = read(local_subgroup, "sites_1")

            ham_op += (coeff,op1,sites1)
        else
            display(local_subgroup)
            error("Invalid number of operators in the subgroup")
        end
    end

    return ham_op

end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, expander::TTN.DefaultExpander)
    group = create_group(parent, name)

    write(group, "p", expander.p)
    write(group, "min", expander.min)
end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, expander::Type{<:TTN.DefaultExpander})
    group = open_group(parent, name)

    p = read(group, "p")
    p == 0.0 && return HDF5.read(parent, name, TTN.NoExpander)
    min = read(group, "min")

    return TTN.DefaultExpander(p; min=min)
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, expander::TTN.NoExpander)
    group = create_group(parent, name)

    write(group, "p", 0.0)
end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, expander::Type{<:TTN.NoExpander})
    group = open_group(parent, name)

    return TTN.NoExpander()
end

function TTN._padding(j::TTN.Index{Vector{Pair{TTN.QN, Int}}}, jp::TTN.Index{Vector{Pair{TTN.QN, Int}}}, p::Int, min::Int; tags = "Padded", kwargs...)
    return TTN._padding(j,jp,p; tags = "Padded", kwargs...)
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, measurement_functions::Vector{NamedTuple})
    group = create_group(parent, name)

    for (i, val) in enumerate(measurement_functions)
        write(group, "measurement_$i", named_tuple_to_dict(val))
    end
end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, measurement_functions::Type{<:Vector{NamedTuple}})
    group = open_group(parent, name)

    measurement_functions = []

    for i in 1:length(group)
        local_subgroup = open_group(group, "measurement_$i")

        measurement = Dict{String,Any}()
        for key in keys(local_subgroup)
            measurement[key] = read(local_subgroup, key)
        end

        append!(measurement_functions, [dict_to_symbols(measurement)])
    end

    return measurement_functions
end





























"fin"
#####################################################
#=

This file is for reading/writing data types in HDF5 format.

=#
######################################################

using HDF5

include("ttn.jl")

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, observer::SavingNRGVarObserver)
	group = create_group(parent, name)
	
	file_path = observer.file_path
	var_tol = observer.var_tol
	nrg = observer.nrg

	write(group, "file_path", file_path)
	write(group, "var_tol", var_tol)
	write(group, "nrg", nrg)
end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, observer::Type{<:SavingNRGVarObserver})
	group = get_group(parent, name)
	
	file_path = read(group, "file_path")
	var_tol = read(group, "var_tol")
	nrg = read(group, "nrg")
	
	SavingNRGVarObserver(file_path, var_tol, nrg)
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
end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, observer::Type{<:SavingMeasurementsObserver})
	group = get_group(parent, name)
	
	file_path = read(group, "file_path")
	var_tol = read(group, "var_tol")
	nrg = read(group, "nrg")
	measurement_functions = read(group, "measurement_functions")
	measurements = read(group, "measurements")
	
	SavingMeasurementsObserver(measurement_functions, measurements, file_path, var_tol, nrg)
end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, ham_op::Sum{Scaled{ComplexF64, Prod{Op}}})
    group = create_group(parent, name)

    for (i, val) in enumerate(ham_op)
        local_subgroup = create_group(group, "op_$i")

        coeff = coefficient(val)
        for (j, t) in enumerate(TTNKit.terms(val))

            # find each local operator and sites
            local_op = TTNKit.which_op(t)
            local_sites = TTNKit.site(t)

            # write data to local_subgroup
            write(local_subgroup, "op_$j", local_op)
            write(local_subgroup, "sites_$j", collect(local_sites))

        end

        write(local_subgroup, "coeff", coeff)
    end

end

function HDF5.read(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, ham_op::Type{<:Sum{Scaled{ComplexF64, Prod{Op}}}})
    group = get_group(parent, name)

    ham_op = TTNKit.OpSum()

    for i in 1:length(group)
        local_subgroup = get_group(group, "op_$i")

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

    ham_op

end

function HDF5.write(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, expander::TTNKit.DefaultExpander)
    group = create_group(parent, name)

    error("Not implemented yet")
end





























"fin"
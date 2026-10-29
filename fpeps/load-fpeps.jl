#####################################################
#=

This file loads everything for fPEPS: packages, other-funcs and the fpeps library files.
Drivers (execute-fpeps.jl, daily-things.jl) include other-funcs/include-other-files.jl and then this;
the library files themselves include nothing and only list what they use

Depends on:
	other-funcs/include-other-files.jl (include_other_files, paths from pwd, so pwd inside fzj/ or local/<user>/)
	QuantumNaturalGradient, QuantumNaturalfPEPS, ITensors, ITensorMPS, JLD2 (fpeps env)

=#
######################################################

using QuantumNaturalGradient, QuantumNaturalfPEPS, ITensors, ITensorMPS, JLD2

include_other_files(["other-funcs/basic-2d-stuff.jl","other-funcs/basic-2d-observables.jl","other-funcs/data-storage-funcs.jl",
	"fpeps/build-ham.jl","fpeps/savings.jl","fpeps/observables.jl"])

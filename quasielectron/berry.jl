using Statistics
include("reverse-flux.jl")
include("laughlin-wavefunc.jl")
include("analysis-functions.jl")


function slice_matrix(data,corrlength)
	if length(size(data)) > 1
		sliced_data = [data[:,j*corrlength] for j in 1:Int(floor(size(data)[2]/corrlength))]
	else
		sliced_data = [data[j*corrlength] for j in 1:Int(floor(size(data)[1]/corrlength))]
	end
	return sliced_data
end

function get_berry_phase(allconfigs,allwavefuncs,qe_loc,clockcount,m=3; kwargs...)
	if_check = get(kwargs, :if_check, false)
	wavefunc_type = get(kwargs, :vers, "P")
	wave_function = wavefunc_type == "P" ? laughlin_wavefunction_girvinjach : reverse_flux_wavefunction
	corrlength = get(kwargs, :corrlength, 1)
	
	configurations = slice_matrix(allconfigs,corrlength)
	wavefunctions = slice_matrix(allwavefuncs,corrlength)
	
	num_parts = length(configurations[1])
	qe_cutoff = get(kwargs, :qe_cutoff, num_parts)	
	rm = sqrt(2*m*num_parts)
	
	dtheta = 2*pi/clockcount
	
	calced_vals = [0.0 for i in 1:length(configurations)]
	for i in 1:length(configurations)
		local_calced_vals = [0.0*im for i in 1:clockcount+1]
		previous_wavefunc = wavefunctions[i]
		for j in 1:clockcount
			new_qe_loc = rm*qe_loc*exp(im*j*dtheta)
			new_wavefunc = wave_function(configurations[i],m; qe_cutoff = qe_cutoff, qe_loc = new_qe_loc)[1]
			#=
			if if_check
				old_wavefunc = wave_function(configurations[i],m; qe_cutoff = qe_cutoff, qe_loc = rm*qe_loc)[1]
				if abs(real(old_wavefunc) - real(wavefunctions[i])) / abs(real(old_wavefunc)) > 0.001
					println("Old Wavefunction calculation off: Recalced = ",old_wavefunc,", Saved = ",wavefunctions[i])
				end
			end
			=#
			#ratio_wavefunc = exp(new_wavefunc - wavefunctions[i])
			ratio_wavefunc = exp(new_wavefunc - previous_wavefunc)
			#local_calced_vals[j] = imag(ratio_wavefunc)#/dtheta
			local_calced_vals[j] = ratio_wavefunc
			previous_wavefunc = new_wavefunc
		end
		calced_vals[i] = -angle(prod(local_calced_vals))
		#println(calced_vals[i])
	end
	
	final_val = mean(calced_vals)
	val_error = std(calced_vals)
	
	return final_val,val_error,calced_vals
end

#
parts = 8
allfiles = find_data_file(Dict([("particles",parts),("qe_cutoff",parts)]),"laugh","jld2","../cluster-data/quasielectron/")
#=
alldata = []
allmetadata = []
for f in [allfiles[1]]
	data = read_data_jld2(f,"../cluster-data/quasielectron/")
	append!(alldata,[data[1]])
	append!(allmetadata,[data[2]])
end
=#

#for cl in [1000]
cl = 100
phases = [0.0 for i in 1:length(allfiles)]
phase_errors = [0.0 for i in 1:length(allfiles)]
qe_radii = [0.0 for i in 1:length(allfiles)]
for (i,f) in enumerate(allfiles)
	println(i/length(allfiles))
	alldata,allmetadata = read_data_jld2(f,"../cluster-data/quasielectron/")
	configs = alldata["configs"]
	wavefuncs = alldata["wavefuncs"]
	berry_data = get_berry_phase(configs,wavefuncs,allmetadata["qe_loc"],10; corrlength = cl,vers="P")
	phases[i] = berry_data[1]
	phase_errors[i] = berry_data[2]
	qe_radii[i] = allmetadata["qe_loc"]
	#scatter([j for j in 1:length(berry_data[3])],berry_data[3])
end
#
println(mean(phase_errors))
plot(qe_radii,phases .+ 1,"-p",label="$cl")
plot(qe_radii,[1/6 for i in 1:length(phases)]) 
xlabel("QE Radius / rm")
ylabel("Aharanov-Bohm Phase")
legend()
#











































"fin"

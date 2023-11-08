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
	throw_counts = [0.0 for i in 1:length(configurations)]
	for i in 1:length(configurations)
		local_calced_vals = [0.0*im for k in 1:clockcount]
		previous_wavefunc = wavefunctions[i]
		throwout_count = 0
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
			#ratio_wavefunc = exp(new_wavefunc - previous_wavefunc)
			ratio_wavefunc = new_wavefunc - previous_wavefunc
			#ratio_wavefunc_exp = exp(new_wavefunc - previous_wavefunc)
			#local_calced_vals[j] = imag(ratio_wavefunc)#/dtheta
			#local_calced_vals_exp[j] = ratio_wavefunc_exp
			if abs(imag(ratio_wavefunc)) < 1.0
				local_calced_vals[j] = ratio_wavefunc
			else
				throwout_count += 1
			end
			previous_wavefunc = new_wavefunc
		end
		#angle(exp(sum(log.(Complex.(local_calced_vals_exp)))))
		#println(angle(prod(local_calced_vals_exp)))
		#println(local_calced_vals_exp[end],", ",local_calced_vals[end])
		#println(angle(prod(exp.(local_calced_vals))))
		#-angle(exp(sum(local_calced_vals)))
		#calced_vals[i] = angle(prod(local_calced_vals))
		#plot([i+k/(2*pi*clockcount) for k in 1:clockcount],imag.(local_calced_vals),"-p",label="$clockcount")
		#plot([i+k/(2*pi*clockcount) for k in 1:clockcount],[0.0 for k in 1:clockcount])
		#plot([i+clockcount/100 for k in 1:1],[sum(imag.(local_calced_vals))],"-p",label="$clockcount")
		#legend()
		throw_counts[i] = throwout_count/clockcount
		calced_vals[i] = imag(sum(local_calced_vals))
		#println(imag(sum(log.(Complex.(local_calced_vals)))))
	end
	
	final_val = mean(calced_vals)
	val_error = std(calced_vals)
	
	return final_val,val_error,calced_vals,throw_counts
end

#
whichtype = "laugh"
m = 1

if whichtype == "rfa"
	parts = 10
	vers = "R"
else
	parts = 8
	vers = "P"
end
allfiles = find_data_file(Dict([("m",m),("particles",parts),("qe_cutoff",parts)]),whichtype,"jld2","../cluster-data/quasielectron/")
#=
alldata = []
allmetadata = []
for f in [allfiles[1]]
	data = read_data_jld2(f,"../cluster-data/quasielectron/")
	append!(alldata,[data[1]])
	append!(allmetadata,[data[2]])
end
=#

#fig2 = figure()
#fig1 = figure()
quadr(x,p) = p[1] .* (x.^2)
for cc in [100,200]
cl = 1000
#cc = 50
phases = [0.0 for i in 1:length(allfiles)]
phase_errors = [0.0 for i in 1:length(allfiles)]
qe_radii = [0.0 for i in 1:length(allfiles)]
throws = [0.0 for i in 1:length(allfiles)]
for (i,f) in enumerate(allfiles)
	println(i/length(allfiles))
	alldata,allmetadata = read_data_jld2(f,"../cluster-data/quasielectron/")
	configs = alldata["configs"]
	wavefuncs = alldata["wavefuncs"]
	berry_data = get_berry_phase(configs,wavefuncs,allmetadata["qe_loc"],cc; corrlength = cl,vers=vers)
	phases[i] = berry_data[1]
	phase_errors[i] = berry_data[2]
	qe_radii[i] = allmetadata["qe_loc"]
	throws[i] = mean(berry_data[4])
	#scatter([j for j in 1:length(berry_data[3])],berry_data[3],label="$cc,$(qe_radii[i])")
	#legend()
end
#
println(mean(phase_errors))
display(throws)
plot(qe_radii,phases,"-p",label="$cc")
#plot(qe_radii,[1/6 for i in 1:length(phases)]) 
xlabel("QE Radius / rm")
ylabel("Aharanov-Bohm Phase")
#
thisfit = LsqFit.curve_fit(quadr,qe_radii[1:5],phases[1:5],[pi/3])
plot(qe_radii,quadr(qe_radii,thisfit.param),label="$(round(thisfit.param[1]/pi,digits=3))")
#
legend()
#
end










































"fin"

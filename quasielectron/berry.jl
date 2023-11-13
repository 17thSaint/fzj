using Statistics
include("reverse-flux.jl")
include("laughlin-wavefunc.jl")
include("analysis-functions.jl")

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

function find_actual_qeloc(configs,expected_loc,rm; kwargs...)
	findqeloc(x,p) = p[1] .* ((x .+ p[2]).^2) .+ p[3]
	raddens = radial_density_full(configs,rm; rend=expected_loc*rm*1.5, points = 300, if_plot=false,titlestring="$expected_loc")
	if expected_loc > 0.4
		shift = 100
		limit = maximum(raddens[2]) * 0.1
		left = findfirst(x -> raddens[2][x] > limit && raddens[2][x+1] < limit,10:length(raddens[2])-1) + 10
		right = findfirst(x -> raddens[2][x+left+shift] < limit && raddens[2][x+left+shift+1] > limit,1:length(raddens[2])-left-shift-1) + left + shift
	else
		left = 2
		right = length(raddens[1])
	end
	findlocfit = LsqFit.curve_fit(findqeloc,raddens[1][left:right],raddens[2][left:right],[1.0,-expected_loc,0.0])
	#plot(raddens[1][left:right],findqeloc(raddens[1][left:right],findlocfit.param),label="$(round(-findlocfit.param[2],digits=3))")
	#legend()
	qe_actual_loc = -findlocfit.param[2]#-findlocfit.param[2] < 0.0 ? expected_loc : -findlocfit.param[2]
	if abs(qe_actual_loc) > 1.2
		return expected_loc
	else
		return qe_actual_loc
	end
end

#
parts = 10
ms = [5]
phases = [[0.0 for i in 1:10] for j in 1:length(ms)]
qe_radii = [[0.0 for i in 1:10] for j in 1:length(ms)]
#part = [8,9,10,11]
for (j,m) in enumerate(ms)
whichtype = "rfa"
#m = 3
#parts = 9
rm = sqrt(2*parts*m)

if whichtype == "rfa"
	#parts = 10
	vers = "R"
else
	#parts = 8
	vers = "P"
end
allfiles = find_data_file(Dict([("m",m),("particles",parts),("qe_cutoff",parts)]),whichtype,"jld2","../cluster-data/quasielectron/")
loc = findfirst(x -> get_params_dict_from_filename(allfiles[x])["qe_loc"] == 0.0,1:length(allfiles))
!isnothing(loc) ? deleteat!(allfiles,loc) : nothing
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
#for cl in [1000]
cl = 100
cc = 100
#phases = [0.0 for i in 1:length(allfiles)]
phase_errors = [0.0 for i in 1:length(allfiles)]
#qe_radii = [0.0 for i in 1:length(allfiles)]
throws = [0.0 for i in 1:length(allfiles)]
for (i,f) in enumerate(allfiles)
	println(i/length(allfiles))
	alldata,allmetadata = read_data_jld2(f,"../cluster-data/quasielectron/")
	therm_time = allmetadata["therm_time"]
	samp_freq = allmetadata["samp_freq"]
	configs = alldata["configs"][:,Int((1000 - therm_time)/samp_freq)+1:end]
	wavefuncs = alldata["wavefuncs"][Int((1000 - therm_time)/samp_freq)+1:end]
	berry_data = get_berry_phase(configs,wavefuncs,allmetadata["qe_loc"],cc; corrlength = Int(cl / samp_freq),vers=vers)
	phases[j][i] = berry_data[1]
	phase_errors[i] = berry_data[2]
	#qe_radii[parts == 9 ? 1 : 2][i] = allmetadata["qe_loc"]
	throws[i] = mean(berry_data[4])
	
	qe_radii[j][i] = find_actual_qeloc(configs,allmetadata["qe_loc"],rm)
	#plot(raddens[1][left:right],findqeloc(raddens[1][left:right],findlocfit.param),label="$(round(-findlocfit.param[2],digits=3))")
	#legend()
	#=
	occrez,bw = get_occupancy(configs,sqrt(2*parts*m),100; if_plot=false,title_string="QE Loc = $(allmetadata["qe_loc"]), m = $m")
	quartersize = Int(size(occrez)[1]/4)
	scatter([allmetadata["qe_loc"]],[bw*(findfirst(x->occrez[2*quartersize,x] ==minimum(occrez[Int(2*quartersize),quartersize:Int((80/25)*quartersize)]),collect(1:4*quartersize))-50)],label="$(allmetadata["qe_loc"]),$m")
	=#
	
	#scatter([j for j in 1:length(berry_data[3])],berry_data[3],label="$cc,$(qe_radii[i])")
	#legend()
end
#
println(mean(abs.(phase_errors ./ phases[j])))
display(throws)
plot(qe_radii[j],phases[j],"-p",label="$parts")
#
#shift = 1.0
#errorbar(qe_radii,phases ./ ((pi*shift) .* (qe_radii.^2)),yerr=[phase_errors ./ ((pi*shift) .* (qe_radii.^2)),phase_errors ./ ((pi*shift) .* (qe_radii.^2))],label="$cl")
#=if cc == 50
plot(qe_radii,[1/3 for j in 1:length(qe_radii)],label="1/3")
end
=#
#plot(qe_radii,[1/6 for i in 1:length(phases)]) 
#xlabel("QE Radius / rm")
#ylabel("Charge of Excitation")
#=
edge = 5
thisfit = LsqFit.curve_fit(quadr,qe_radii[j][1:edge],phases[j][1:edge],[pi/3])
plot(qe_radii[j],quadr(qe_radii[j],thisfit.param),label="$(round(thisfit.param[1]/pi,digits=3))")
#
legend()
=#
end
#=

quadr(x,p) = p[1] .* (x.^2)
slopes = [0.0 for i in 1:3]
edges = [5,5,8]
givenslope = [36,12,36/5]
for j in 1:3
plot(qe_radii[j],phases[j],"-p",label="1/$(ms[j])")
edge = edges[j]
thisfit = LsqFit.curve_fit(quadr,qe_radii[j][1:edge],phases[j][1:edge],[pi/3])
if j == 1
plot(qe_radii[j],quadr(qe_radii[j],[givenslope[j]*pi]),c="k",label="TH")#"$(round(thisfit.param[1]/pi,digits=3))")
else
plot(qe_radii[j],quadr(qe_radii[j],[givenslope[j]*pi]),c="k")#"$(round(thisfit.param[1]/pi,digits=3))")
end
slopes[j] = thisfit.param[1]/pi
end
xlabel("QE Radius / rm")
ylabel("Berry Phase")
legend()
ylim(0,60)
ratio13 = slopes[1] / slopes[2]#mean(phases[1][1:10] ./ phases[2][1:10])
ratio15 = slopes[1] / slopes[3]#mean(phases[1][1:10] ./ phases[3][1:10])
title("Laughlin v=1,1/3,1/5")#: 1/3 Ratio = $(round(ratio13,digits=4)), 1/5 Ratio = $(round(ratio15,digits=4))")
=#
#=
fig = figure()
plot(qe_radii[1],phases[1] ./ phases[2],"-p")
xlabel("QE Radius / rm")
=#
#=
quadr(x,p) = p[1] .* (x.^2)
for parts in [9,10]
fig = figure()
plot(qe_radii[parts == 9 ? 1 : 2],phases[parts == 9 ? 1 : 2],"-p",label="$parts")
xlabel("QE Radius / rm")
ylabel("Charge of Excitation")


thisfit = LsqFit.curve_fit(quadr,qe_radii[parts == 9 ? 1 : 2],phases[parts == 9 ? 1 : 2],[pi/3])
plot(qe_radii[parts == 9 ? 1 : 2],quadr(qe_radii[parts == 9 ? 1 : 2],thisfit.param),label="$(round(thisfit.param[1]/pi,digits=3))")
legend()
end
=#






































"fin"

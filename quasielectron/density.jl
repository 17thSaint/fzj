include("reverse-flux.jl")
include("laughlin-wavefunc.jl")
include("analysis-functions.jl")

dataloc = "../cluster-data/quasielectron/"
m = 3
parts = 10
rm = sqrt(2*parts*m)
whichtype = "laugh"
bins = 100
rs = [-1.2 + (i-1)*2.4/(bins-1) for i in 1:bins]

allfiles = find_data_file(Dict([("m",m),("particles",parts)]),whichtype,"jld2",dataloc)
display(allfiles)

gsfile = allfiles[end]
gsdata,gsmetadata = read_data_jld2(gsfile,dataloc)

gs_density = get_occupancy(gsdata["configs"],rm,bins; max_x = 1.2, max_y = 1.2, if_plot = false)
nn = integrate(rs,gs_density[1][Int(bins/2),:])
plot(rs,gs_density[1][Int(bins/2),:] ./ nn)
title("GS")

all_excessdensity = [0.0 for i in 1:length(allfiles)-1]
qe_radii = [0.0 for i in 1:length(allfiles)-1]
for i in length(allfiles)-1:length(allfiles)-1
	f = allfiles[i]
	alldata,allmetadata = read_data_jld2(f,dataloc)
	qe1_density = get_occupancy(alldata["configs"],rm,bins; max_x = 1, max_y = 1, if_plot = false)
	#fig = figure()
	nn2 = integrate(rs,qe1_density[1][Int(bins/2),:])
	plot(rs,qe1_density[1][Int(bins/2),:] ./ nn2)
	title("Single QE")

	rig = figure()
	plot(rs,(qe1_density[1][Int(bins/2),:] ./ nn2) .- (gs_density[1][Int(bins/2),:] ./ nn),label="$(allmetadata["qe_loc"])")
	xlabel("Radius / rm")
	legend()
	#=
	intexdens = integrate([(j-1)*rm/(0.5*bins-1) for j in 1:Int(bins/2)],qe1_density[1][Int(bins/2),Int(bins/2)+1:end] .- gs_density[1][Int(bins/2),Int(bins/2)+1:end])
	all_excessdensity[i] = intexdens
	qe_radii[i] = allmetadata["qe_loc"]
	=#
end

#plot(qe_radii,all_excessdensity,"-p")


































"fin"

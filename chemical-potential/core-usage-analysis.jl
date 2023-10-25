
dataloc = "../cluster-data/chemical-potential/data-2023-10-25_13-59/"

#basename = "ttn-maxocc-1-chem-1.0-alpha-0.25-t-0.5-layers-4"
#metadatas = read_data_jld2(basename,dataloc)[2]

files = readdir(dataloc)

alllines = []
for f in files
	fo = open(dataloc*f)
	lines = readlines(fo)
	append!(alllines,[lines])
end

cores = [Int(parse(Float64,split(alllines[i][5],">")[end])) for i in 1:length(alllines)]
times = [mean([parse(Float64,split(split(alllines[j][14 + 3*i]," ")[end],"s")[1]) for i in 1:4]) for j in 1:length(alllines)]
scatter(cores,times)
xlabel("Cores Being Used")
ylabel("Average Time")
title("Avg Time for single sweep of BondDim 100 4 Layer TTN vs Cores used")
#=
times = []
for m in metadatas
	append!(times,[m["runtime"]])
	append!(cores,[m["open_cores"]])
end
end
=#




























"fin"

using DelimitedFiles

which_file = "../cluster-data/synth-dims/output-layers-4-filling-.5-mdim-20-hopping_anisotropy-0.5-lr-3-onsite_strength-5.0.txt"

all_lines = open(which_file) do file
    readlines(file)
end

profile_start = findfirst(i -> occursin("Running time",all_lines[i]), 1:length(all_lines)) + 3
full_profile = all_lines[profile_start:end-1]

data_dict = Dict{Int64,String}()
for line in full_profile
    split_profile = split(line,"╎")
    profile_counts = tryparse(Int64,split_profile[1])
    if profile_counts != nothing && profile_counts >= 100
        data_dict[profile_counts] = split(split_profile[end],";")[2]
    end
end

for (key,value) in data_dict
    println(key," ",value)
end
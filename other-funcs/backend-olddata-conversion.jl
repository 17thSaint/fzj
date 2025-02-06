#####################################################
#=

This file is for converting the old Backend TTN data type to the main_new version without the Backend variable

    Need to convert:
        TreeTensorNetwork
        Lattice
        Network
        TPO
        MPO

=#
######################################################

include("data-storage-funcs.jl")
include("../review-practice-codes/ttn.jl")
include("../review-practice-codes/hdf5-highlevel-saving.jl")

# convert all jld2 files to hdf5
if true
    starting_location = pwd()

    root_dir = get_folder_location("cluster-data/synth-dims")

    for (root,dirs,files) in walkdir(root_dir)

        println("Current directory: $root")

        if occursin("torus",root)
            println("SKIPPING TORUS FOLDER")
            continue
        end

        for f in files

            if occursin(".jld2",f) && occursin("ttn",f)

                println("Working on file: $(joinpath(root,f))")

                try

                    d,m = read_data_jld2(joinpath(root,f))
                    write_data_hdf5(joinpath(root,f),d,m)

                catch

                    println("File $f broke")
                    
                    h5open(starting_location*"/files-that-broke.h5","r+") do binaryfile
                        i = length(keys(binaryfile)) + 1
                        write(binaryfile, "f-$i", joinpath(root,f))
                    end

                end
                
            end

        end
    end
end

# do just ttn alone by making new data file
if false
    all_files = readdir()
    filter!(x -> occursin(".jld2",x),all_files)
    filter!(x -> occursin("ttn",x),all_files)
    f = all_files[1]

    d,m = read_data_jld2(f);

    write_data_hdf5("ttn-testing-true",d,m)

    #=h5open("ttn-testing-true.h5","w") do f
		g_alldata = create_group(f,"all_data")
		for (datum_key,datum) in d
			isnothing(datum) && continue
			write(g_alldata, datum_key, datum)
		end
	end=#

end








































"fin"
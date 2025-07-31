#####################################################
#=

This file contains the tests and functions for multithreading with GPUs

Depends on:
    

=#
######################################################

using CUDA,TTN,LinearAlgebra,Statistics

function indexed_contractall(tnlist::Vector{TTN.ITensor})
    result::Float64 = 0.0
    for i in 1:1:Int(length(tnlist)/2)
        result += TTN.scalar(tnlist[2i-1] * tnlist[2i])
    end
    return result
end

function threaded_contractall(tnlist::Vector{TTN.ITensor})
    result::Float64 = 0.0
    Threads.@spawn for i in 1:Int(length(tnlist)/2)
        result += TTN.scalar(tnlist[2i-1] * tnlist[2i])
    end
    return result
end

avg_counts = 50
list_length = 2

allsizes = collect(50:10:150)
alltimes_cpu = zeros(Float64, length(allsizes))
error_cpu = zeros(Float64, length(allsizes))
alltimes_gpu = zeros(Float64, length(allsizes))
error_gpu = zeros(Float64, length(allsizes))

for (idx,dimsize) in enumerate(allsizes)
    println("Testing with size: ", dimsize)

    local_times_cpu = zeros(Float64, avg_counts)
    local_times_gpu = zeros(Float64, avg_counts)

    i,j = TTN.Index(dimsize,"i"), TTN.Index(dimsize*10,"k")

    for x in 1:avg_counts

        # Create random tensors
        tnlist = [TTN.randomITensor(i,j) for m in 1:list_length]

        # run scalar indexing contraction on CPU
        cpu_time = Base.@elapsed indexed_contractall(tnlist)
        local_times_cpu[x] = cpu_time

        # move tensors to GPU
        gpu_tnlist = [TTN.adapt(CuArray,tn) for tn in tnlist]

        # run threaded contraction on GPU
        gpu_time = CUDA.@elapsed threaded_contractall(gpu_tnlist)
        local_times_gpu[x] = gpu_time

    end
    alltimes_cpu[idx] = mean(local_times_cpu)
    error_cpu[idx] = std(local_times_cpu)
    alltimes_gpu[idx] = mean(local_times_gpu)
    error_gpu[idx] = std(local_times_gpu)
end

allspeedups = alltimes_cpu ./ alltimes_gpu
display(allspeedups)
speeduperror = allspeedups .* sqrt.((error_cpu ./ alltimes_cpu.^2 + error_gpu ./ alltimes_gpu.^2))
display(speeduperror)

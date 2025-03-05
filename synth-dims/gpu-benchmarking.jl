#####################################################
#=

This file contains GPU testing and benchmarking on TTNs

=#
######################################################

#using Pkg
#Pkg.activate(".")

include("../other-funcs/include-other-files.jl")
include_other_files(["synth-dims/long-range-ttn.jl"])

#using CUDA
using LinearAlgebra

mutable struct BenchmarkObserver <: TTN.AbstractObserver
	sweep_count::Int
    sweep_times::Vector{Float64}
 
    BenchmarkObserver(sweep_count=5) = new(sweep_count,[])
end

function TTN.measure!(o::BenchmarkObserver; kwargs...)
    append!(o.sweep_times,[kwargs[:dt]])
end

function TTN.checkdone!(o::BenchmarkObserver;kwargs...)
	sh = kwargs[:sweep_handler]
    if (kwargs[:sweep_handler]).current_sweep == o.sweep_count
        return true
    else
        return false
    end
end

function initialize_benchmark(pu_type::String,net::TTN.AbstractNetwork; kwargs...)

    # build TTN
    psi = RandomTreeTensorNetwork(net; maxdim=25)
    if pu_type == "gpu"
        psi = gpu(psi)
    end
    lat = physical_lattice(net)

    # make Ham
    params_dict = Dict([("hopping_anisotropy",1.0),("if_check_fluxes",false),("particles",10),("layers",number_of_layers(psi)),("filling",0.5),("onsite_strength",0.0),("lr","all"),("if_periodic_phys",false),("if_periodic_synth",false)])
    model_paras = get_normal_model_params(params_dict)
    ham = long_range_HH_ham(net,model_paras[:ts],model_paras[:alpha]; model_paras...)
    tpo = TTN.TPO(ham,lat)

    # do DMRG
    sp = TTN.dmrg(psi,tpo; number_of_sweeps=5, maxdims=25, cutoff=0, eigsolve_krylovdim=5, eigsolve_verbosity=0)

    println("Finished $pu_type initialization")
end

function run_benchmark(pu_type::String,mdim::Int,net::TTN.AbstractNetwork,filepath::String; kwargs...)
    # get kwargs
    opl::Int = get(kwargs,:opl,1)
    n::Int = get(kwargs,:n,8)
    stren::Float64 = get(kwargs,:stren,0.0)
    if_periodic::Bool = get(kwargs,:if_periodic,false)
    num_sweeps::Int = get(kwargs,:num_sweeps,5)

    # build TTN
    psi = RandomTreeTensorNetwork(net; maxdim=mdim)
    if pu_type == "gpu"
        psi = gpu(psi)
    end
    lat = physical_lattice(net)

    opl > 1 && println("Starting linkdim = $(maxlinkdim(psi))")

    # make Ham
    params_dict = Dict([("hopping_anisotropy",1.0),("if_check_fluxes",false),("particles",n),("layers",number_of_layers(psi)),("filling",0.5),("onsite_strength",stren),("lr","all"),("if_periodic_phys",if_periodic),("if_periodic_synth",if_periodic)])
    model_paras = get_normal_model_params(params_dict)
    ham = long_range_HH_ham(net,model_paras[:ts],model_paras[:alpha]; model_paras...)
    tpo = TTN.TPO(ham,lat)

    # make observer
    obs = BenchmarkObserver(num_sweeps)

    # do DMRG
    time_start = time()
    sp = TTN.dmrg(psi,tpo; observer=obs, number_of_sweeps=num_sweeps, maxdims=mdim, cutoff=0, eigsolve_krylovdim=5, eigsolve_verbosity=0)
    final_runtime = (time() - time_start) / num_sweeps
    
    append!(obs.sweep_times,[final_runtime])

    final_mdim = maxlinkdim(sp.ttn)

    new_data = Dict([("$final_mdim",obs.sweep_times)])
    modify_data_hdf5(new_data,filepath,"all_data")

    opl > 0 && println("GPU: mdim = $final_mdim, time = $final_runtime")

    return obs,final_mdim
end

function get_pu_info(pu_type::String)
    if pu_type == "gpu"
        pu_info = CUDA.versioninfo()
    elseif pu_type == "cpu"
        pu_info = read(`lscpu`, String)
    end
    return pu_info
end

#Dict([("benchmark_type","cpu"),("min_mdim",50),("max_mdim",50),("count_mdim",1)])
args_dict = make_args_dict(ARGS)

# set benchmarking parameters
benchmark_type::String = args_dict["benchmark_type"]
min_mdim::Int = args_dict["min_mdim"]
max_mdim::Int = args_dict["max_mdim"]
count_mdim::Int = args_dict["count_mdim"]
diff_mdim::Int = Int(floor((max_mdim - min_mdim) / count_mdim))
thread_count::Int = get(args_dict,"thread_count",1)

# set system and Ham parameters
lx = get(args_dict,"lx",4)
ly = get(args_dict,"ly",lx)
n = Int(ceil(0.5*lx*ly))
layers = Int(log(2,lx*ly))
stren = get(args_dict,"stren",0.0)
if_periodic = get(args_dict,"if_periodic",false)
max_occ = 2
etol = 1e-5
conserve_qns = true
dataloc = get_folder_location("cluster-data/gpu-benchmarking")

pu_info = get_pu_info(benchmark_type)
metadata = Dict([("pu_info",pu_info),("thread_count",thread_count),("lx",lx),("ly",ly),("layers",layers),("stren",stren),("if_periodic",if_periodic),("max_occ",max_occ),("etol",etol),("conserve_qns",conserve_qns)])
filename = "$(benchmark_type)-benchmarking-data-layers-$layers-startmdim-$min_mdim-endmdim-$max_mdim-countmdim-$count_mdim.h5"

filepath = joinpath(dataloc,filename)
filename = write_data_hdf5(filepath,Dict(),metadata)

# create network and lattice
net = BinaryRectangularNetwork((lx,ly), "Boson"; conserve_qns=conserve_qns, dim = max_occ+1)

# initialization step
initialize_benchmark(benchmark_type,net)

# run benchmark
if benchmark_type == "cpu"
    BLAS.set_num_threads(thread_count)
end
mdims = collect(min_mdim:diff_mdim:max_mdim)
alltimes = zeros(Float64,length(mdims))
allmdims = zeros(Float64,length(mdims))
for (idx,mdim) in enumerate(mdims)
    alltimes[idx],allmdims[idx] = run_benchmark(benchmark_type,mdim,net,filepath)
end


#= plot benchmarking results
if false
    using PyPlot
    dataloc = get_folder_location("cluster-data/gpu-benchmarking")
    all_files = readdir(dataloc)
    filter!(x -> occursin("h5",x),all_files)
    filter!(x -> occursin("layers-6",x),all_files)

    cpu_files = filter(x -> occursin("cpu",x),all_files)
    for (idx,f) in enumerate(cpu_files)
        d,m = read_data(joinpath(dataloc,f))
        alltimes = []
        allmdims = []
        for (k,v) in d
            push!(alltimes,v)
            push!(allmdims,parse(Float64,k))
        end
        scatter(allmdims,alltimes,label="CPU",color="r")
    end

    gpu_files = filter(x -> occursin("gpu",x),all_files)
    for (idx,f) in enumerate(gpu_files)
        d,m = read_data(joinpath(dataloc,f))
        alltimes = []
        allmdims = []
        for (k,v) in d
            push!(alltimes,v)
            push!(allmdims,parse(Float64,k))
        end
        scatter(allmdims,alltimes,label="GPU",color="g")
    end

    xlabel("Max Bond Dimension")
    ylabel("Time (s)")
    title("AVG Time per DMRG Sweep vs Max Bond Dimension")
    legend()
    xscale("log")
    yscale("log")
end=#



































"fin"
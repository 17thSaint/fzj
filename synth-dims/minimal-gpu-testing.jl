#####################################################
#=

This file contains testing of GPU functions

=#
######################################################




include("../other-funcs/include-other-files.jl")
include_other_files(["synth-dims/long-range-ttn.jl","other-funcs/basic-2d-observables.jl"])

# testing tensor times mpo which has DenseVector type issues
if true
    #=
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))
    dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    pdict = Dict([("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    filter!(x -> !occursin("if_synth_rectangle",x),all_files)
    f = all_files[1]
    f = all_files[1]
    d,m = read_data(dataloc * "/" * f; output_level=0)
    psi = gpu(d["ttn"])
    lat = TTN.physical_lattice(psi.net)
    mapss = zigzag_curve(lx,ly)
    fourpt = four_point_mpo(psi; momentum1 = [0.0,0], momentum2 = [0.0,0.0], mapping = mapss)
    =#

    tn = TTN.ITensors.tensor(psi[1,1])
    mp = TTN.ITensors.tensor(fourpt[1])
    la,lb = TTN.ITensors.compute_contraction_labels(TTN.inds(tn),TTN.inds(mp))
    lr = TTN.NDTensors.contract_labels(la,lb)

    tr = TTN.NDTensors.contraction_output_type(typeof(tn),typeof(mp),lr)
    indsR = TTN.NDTensors.contract_inds(TTN.inds(tn), la, TTN.inds(mp), lb, lr)

    bosR, cplan = TTN.NDTensors.contract_blockoffsets(
        TTN.NDTensors.blockoffsets(tn),
    TTN.inds(tn),
    la,
    TTN.NDTensors.blockoffsets(mp),
    TTN.inds(mp),
    lb,
    indsR,
    lr,
  )
    R = TTN.NDTensors.similar(tr, bosR, indsR)

end







































"fin"
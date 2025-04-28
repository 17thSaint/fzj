#####################################################
#=

This file contains testing of GPU functions

=#
######################################################




include("../other-funcs/include-other-files.jl")
include_other_files(["synth-dims/long-range-ttn.jl","other-funcs/basic-2d-observables.jl"])

# testing tensor times mpo which has DenseVector type issues
if true
    #
    lx,ly,n = 8,4,4
    layers = Int(log(2,lx*ly))
    #=dataloc = get_folder_location("cluster-data/synth-dims/excited-states")
    pdict = Dict([("layers",layers),("particles",n),("if_periodic_phys",true),("if_periodic_synth",true),("hopping_anisotropy",1.0)])
    all_files = find_data_file(pdict,"ttn",dataloc)
    filter!(x -> !occursin("if_synth_rectangle",x),all_files)
    f = all_files[1]
    f = all_files[1]
    d,m = read_data(dataloc * "/" * f; output_level=0)
    psi = gpu(d["ttn"])
    =#

    #lat = TTN.physical_lattice(psi.net)
    #mapss = zigzag_curve(lx,ly)
    #fourpt = four_point_mpo(psi; momentum1 = [0.0,0], momentum2 = [0.0,0.0], mapping = mapss, coeff_kwargs=(Lx=lx,Ly=ly,))
    #fourpt_wrapped = easy_mpowrapper(fourpt, lat; mapping=mapss)
    #CUDA.@allowscalar rez = calculate_mpo_expectation(psi, fourpt_wrapped)
    #tlist = [psi[1,1], TTN.adapt(CuArray,fourpt[1]), TTN.dag(TTN.prime(psi[1,1]))]

    net = BinaryRectangularNetwork((lx,ly), "Boson"; conserve_qns=true, dim = 2)
    psi = RandomTreeTensorNetwork(net; maxdim=25)
    psi = gpu(psi)

    tn1 = psi[1,1]

    # the issue seems to be independent of the MPO in the middle, it's the inner product that can't be done


    #=psi_phys_index1 = TTN.Index([TTN.QN("Number",0)=>1,TTN.QN("Number",1)=>1]; tags="Site 1")
    psi_phys_index2 = TTN.Index([TTN.QN("Number",0)=>1,TTN.QN("Number",1)=>1]; tags="Site 2")
    psi_up_index = TTN.Index([TTN.QN("Number",0)=>1,TTN.QN("Number",1)=>2,TTN.QN("Number",2)=>1]; tags="Upsite")
    ttn_tensor = TTN.ITensors.randomITensor(ComplexF64,[psi_phys_index1,psi_phys_index2,TTN.dag(psi_up_index)])
    ttn_tensor = TTN.adapt(CuArray,ttn_tensor)=#

    #=phys_ind = psi_phys_index1
    virt_index = TTN.Index([TTN.QN("Number",0)=>1,TTN.QN("Number",1)=>1]; tags="virtual")
    blahmpo = TTN.ITensors.randomITensor(ComplexF64,[TTN.dag(phys_ind),TTN.prime(phys_ind),virt_index])
    blahmpo = TTN.adapt(CuArray,blahmpo)=#

    #tlist = [ttn_tensor, blahmpo, TTN.dag(TTN.prime(ttn_tensor))]
    #CUDA.@allowscalar rez = contract(tlist)

    #=pos = (2,1)
    net = psi.net
    # getting the child tensor
    tn_child = psi[pos]
    # getting the parent node
    pos_parent = TTN.parent_node(net, pos)
    tn_parent = psi[pos_parent]
    # the left index for the splitting is simply the not commoninds of tn_child
    idx_r = TTN.commonind(tn_child, tn_parent)
    #idx_l = uniqueinds(tn_child, tn_parent)
    idx_l = TTN.uniqueinds(tn_child, idx_r)
    #Q,R = qr(tn_child, idx_l; tags = tags(idx_r))
    Q,R = TTN.factorize(tn_child, idx_l; tags = TTN.tags(idx_r))
    # handles large normed TTN's. Specially for random initialization
    if true
        R .= R./TTN.norm(R)
    end
    res = R*tn_parent=#
    # R is not on GPU memory, tn_child is GPU but the resulting Q,R are not on the GPU, because Factorize uses scalar indexing and happens on the GPU
    # however when R is put on GPU afterwards it still gets an illegal memory error



end







































"fin"
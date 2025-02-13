using ITensorMPS
using LinearAlgebra
using LaTeXStrings
using KrylovKit
using Printf, Random

create_wavefunction(sz::NTuple{D, Int}) where{D} = create_wavefunction(ComplexF64, sz)
create_wavefunction(elT, sz::NTuple{D,Int}) where{D} = normalize!(randn(elT, sz))

function patron_application!(ttn::TTNKit.TreeTensorNetwork, wf_coefs::Array, op_ins::String; maxdim::Int = maxlinkdim(ttn), normalize::Bool = true)
    net = TTNKit.network(ttn)
    lat = TTNKit.physical_lattice(net)

    all(size(lat) .== size(wf_coefs)) || error("Trying to apply a patron wavefunction of dimensionality $(size(wf_coefs)) to a TTN defined on a lattice $(size(lat))")


    patron_mpo = wf_mpo(wf_coefs, net, op_ins)
    
    for p in eachindex(TTNKit.lattice(net, 1))
        ttn = TTNKit.move_ortho!(ttn, (1,p))
        pc = TTNKit.child_nodes(net, (1, p))
        Tpat = map(j -> patron_mpo[j], last.(pc))
        
        pr = TTNKit.parent_node(net, (1,p))
        #println("Site $p has parent ",pr)
        #println("Child nodes with lattice coords are $(pc[1][2]) at $(TTNKit.coordinate(lat,pc[1][2])) and $(pc[2][2]) at $(TTNKit.coordinate(lat,pc[2][2]))")
        Tt  = ttn[(1, p)]

        rind = TTNKit.commonind(Tt, ttn[pr])
        linds = TTNKit.uniqueinds(Tt, rind)
        Tn = TTNKit.noprime(TTNKit.contract(Tt, Tpat...))

        #println("Multiplication worked")
        #if p == 2
        #    return Tn,linds,rind
        #end
        A, R = TTNKit.factorize(Tn, linds, maxdim = maxdim, tags = TTNKit.tags(rind))
        
        #=
        println("Physical Site $p")
        println("Sizes: ",", Tn",length(TTNKit.inds(Tn)),", R",length(TTNKit.inds(R)),", Pr",length(TTNKit.inds(ttn[pr])),", A",length(TTNKit.inds(A)))
        println(", Tn",TTNKit.tags.(TTNKit.inds(Tn)),", R",TTNKit.tags.(TTNKit.inds(R)),", Pr",TTNKit.tags.(TTNKit.inds(ttn[pr])),", A",TTNKit.tags.(TTNKit.inds(A)))
        =#
        ttn[(1,p)] = A
        ttn[pr] = ttn[pr] * R    
        
        # moving the orhto center
        ttn.ortho_center[1] = pr[1]
        ttn.ortho_center[2] = pr[2]
        # this has to be deleted in the future... dont need this anymore
        ttn.ortho_direction[1][p] = TTNKit.number_of_child_nodes(net, (1,p)) + 1
        ttn.ortho_direction[pr[1]][pr[2]] = -1
    end
    # now move to the higher layers layer by layer, excluding the top node
    TTNKit.ITensors.set_warn_order(20)
    for ll in 2:TTNKit.number_of_layers(net)-1
        for p in 1:TTNKit.number_of_tensors(net, ll)
            pp = (ll, p)
            ttnc = TTNKit.move_ortho!(ttn, pp)
            
            Tt = ttn[pp]

            pr = TTNKit.parent_node(net, pp)

            pc = TTNKit.child_nodes(net, pp)
            
            linds = map(p -> TTNKit.commonind(Tt, ttn[p]), pc)
            ptags = "Link,nl=$(ll),np=$(p)"#tags(commonind(Tt, ttn[pr]))
        
            A, R = factorize(Tt, linds; maxdim = maxdim, tags = ptags)
            #=
            println(ll,", ",p)
            println("With Parent $pr")
            println("Sizes: ",", R",length(TTNKit.inds(R)),", Tt",length(TTNKit.inds(Tt)),", A",length(TTNKit.inds(A)))
            println(", R",TTNKit.tags.(TTNKit.inds(R)),", Tt",TTNKit.tags.(TTNKit.inds(Tt)),", A",TTNKit.tags.(TTNKit.inds(A)))
            =#
            ttn[pp] = A
            ttn[pr] = ttn[pr] * R
        
            # moving the orhto center
            ttn.ortho_center[1] = pr[1]
            ttn.ortho_center[2] = pr[2]
            # this has to be deleted in the future... dont need this anymore
            ttn.ortho_direction[ll][p] = TTNKit.number_of_child_nodes(net, (ll,p)) + 1
            ttn.ortho_direction[pr[1]][pr[2]] = -1
        end
    end

    if normalize
        tpnd = (TTNKit.number_of_layers(net), 1) 
        T = ttn[tpnd]
        ttn[tpnd] = T/norm(T)
    end

    return ttn
end

function get_2dttn_path(size) # written by ChapGPT 29.05.2023
    path = []

    function traverse_lattice(x_start, y_start, x_end, y_end)
        if x_start == x_end && y_start == y_end
            push!(path, (x_start, y_start))
        else
            x_mid = (x_start + x_end) ÷ 2
            y_mid = (y_start + y_end) ÷ 2
            traverse_lattice(x_start, y_start, x_mid, y_mid)
            traverse_lattice(x_mid + 1, y_start, x_end, y_mid)
            traverse_lattice(x_start, y_mid + 1, x_mid, y_end)
            traverse_lattice(x_mid + 1, y_mid + 1, x_end, y_end)
        end
    end

    traverse_lattice(1, 1, size, size)
    return path
end

function get_site_number_square(site_label, side_length)
    x,y = site_label
    site_number = (y - 1) * side_length + x
    if site_number > side_length^2
    	println("ERROR: Outside Square")
    	return Int(site_number)
    end
    return Int(site_number)
end

function get_site_number(x,y,a,b) # written by ChatGPT 12.06.2023
    # Check if the coordinates are within the lattice range
    if x < 1 || x > a || y < 1 || y > b
        println("Outside Lattice")
        return
    end

    # Calculate the number of complete rows and columns
    complete_rows = div(x - 1, a)
    complete_cols = div(y - 1, b)

    # Calculate the number of remaining elements in the incomplete row and column
    remaining_row = rem(x - 1, a)
    remaining_col = rem(y - 1, b)

    # Calculate the number of the lattice site
    if complete_rows % 2 == 0
        site_number = complete_rows * a + complete_cols * b + 1 + remaining_row * b + remaining_col
    else
        site_number = complete_rows * a + (a - remaining_row) * b + 1 + remaining_col
    end

    return site_number
end

function add_rect_block(path,edge_size)
    second_block = []
    for i in 1:length(path)
        push!(second_block,(path[i][1]+edge_size,path[i][2]))
    end
    return [path; second_block]
end

function ttn_2d_mapping(size)
    path = get_2dttn_path(minimum(size))
    if size[1] != size[2]
        path = add_rect_block(path,minimum(size))
    end
    map2d::Vector{Int64} = []
    for i in 1:length(path)
        
        site_idx = Int(get_site_number(path[i][2],path[i][1],size[2],size[1]))
        append!(map2d,[site_idx])
    end
    return map2d
end

function construct_first_layer(mpo,mapping,net)
    bEnvironment = Vector{ITensor}(undef, 1) 

    bEnvironment = map(eachindex(net,1)) do pp
        chdnds = TTNKit.child_nodes(net, (1,pp))
        map(1:TTNKit.number_of_child_nodes(net, (1,pp))) do nn
          mpo[TTNKit.inverse_mapping(mapping)[chdnds[nn][2]]]
        end
    end
    
    return bEnvironment
end

function wf_mpo(wf, net, op_ins)
    ampo = TTNKit.OpSum()
    lat = TTNKit.physical_lattice(net)
    for pci in keys(wf)
        pt = to_indices(wf, (pci,))

        plin = TTNKit.linear_ind(lat, pt)
        ampo += (wf[pci], op_ins, plin)
    end
    
    mapping = ttn_2d_mapping(size(lat))
    
    # build MPO out of OpSum
    idx_lat = map(mapping) do pos
        TTNKit.hilbertspace(lat[pos])
    end
    mpo = TTNKit.MPO(ampo,idx_lat)
    
    rez_data = construct_first_layer(mpo,mapping,net)
    
    remapped_mpo::Vector{ITensor} = []
    for i in 1:size(rez_data,1)
        for j in 1:size(rez_data[i],1)
            append!(remapped_mpo,[rez_data[i][j]])
        end
    end
    return remapped_mpo
    
end

























module MajoranaLadder

using Combinatorics
using Printf, LinearAlgebra, Arpack

struct LocalSpace
    L::Int64 # physical length of the system
    N::Int64 # number of electrons in the system

    dimHilb :: Int64 # total dimension of hilbertspace
    totalSites :: Int64 # Size with enrolled unit cell

    basis_states :: Array{String,2} # basis states in string rep
    basis_ints   :: Array{Int,2}    # basis states in integer rep
    basis_dict   :: Dict{String,Int} # fast access to the position of basis in array

    t   :: Float64
    w1  :: Float64
    w2  :: Float64
    u   :: Float64

    name :: String
end


function LocalSpace_init(L::Int64,N::Int64, paraDict :: Dict{String, Float64}; use_z2::Bool = true, z_sector::Int64 = 0)

    bit_strings, bit_ints, basis_dict = generate_basis(L,N, use_z2 = use_z2, z_sector = z_sector)

    dimHilb = size(bit_strings,2)

    t  = get(paraDict, "t", 1.0)
    w1 = get(paraDict, "w1", 0.0)
    w2 = get(paraDict, "w2", w1)
    u  = get(paraDict, "u", 0.0)

    majLad = LocalSpace(L,N,dimHilb, L*2, bit_strings, bit_ints, basis_dict, t, w1, w2, u, "MajoranaLadder")
    return majLad
end

precompile(LocalSpace_init,(Int,Int,Bool, Int))
function generate_basis(L::Int64,N::Int64; use_z2::Bool = true, z_sector::Int64 = 0)

    @printf "Strart to calulate basis states for L=%.0f and N=%.0f\n" L N
    @printf "Use reduction of hilbertspace by Z2 symmetry: %s\n" use_z2

    z_sector_ = mod(z_sector,2)

    n_unit = 2 # unit cell size
    totalSites = L*n_unit
    n_states = binomial(totalSites,N)
    @printf "Numbers of basis states: %.0f\n" n_states

    bit_strings = Array{String,2}(undef,totalSites,n_states)
    bit_ints    = Array{Int,2}(undef,   totalSites,n_states)

    basis_dict  = Dict{String,Int}()

    subsets = combinations([1:1:totalSites;],N)

    jj = 1
    for subs in subsets
        bit_strings[:,jj]  = fill("0",totalSites)
        bit_ints[:,jj]     = fill(0,  totalSites)
        for en in subs
            bit_strings[en,jj]  = "1"
            bit_ints[en,jj]     =  1
        end
        basis_dict[join(bit_strings[:,jj])] = jj
        jj = jj + 1
    end

    if(use_z2)
        @printf "\nReducing the basis state to the Z2 sector: %.0f\n" z_sector_
        #bits_new = Array{String,2}(undef,0,L)
        # first count the number of surviving states
        n_states_new = 0
        match_ind    = Array{Int,1}(undef,0)
        for jj = 1:n_states
            sum = 0
            for kk = 1:2:totalSites
                sum = sum + bit_ints[kk,jj]
            end
            if(mod(sum,2) == z_sector_)
                n_states_new += 1
                # remember the matching index
                push!(match_ind, jj)
            end
        end
        @printf "\nDimension of total Hilbertspace after reduction: %.0f\n" n_states_new

        # create the new Arrays
        str_new  = Array{String,2}(undef, totalSites, n_states_new)
        int_new  = Array{Int,2}(undef, totalSites, n_states_new)
        new_dict = Dict{String,Int}()
        # and fill it the the matching numbers
        for jj = 1:n_states_new
            str_new[:,jj] = bit_strings[:,match_ind[jj]]
            int_new[:,jj] = bit_ints[:,match_ind[jj]]
            new_dict[join(str_new[:,jj])] = jj
        end
        bit_strings = str_new
        bit_ints    = int_new
        basis_dict  = new_dict
    end

    return bit_strings, bit_ints, basis_dict
end

precompile(generate_basis,(Int,Int,Bool, Int))

function applyHam(n::Int, model::LocalSpace)
    # n :: active basis state to act on
    # model :: holding all informationa about the system

    # load the bit representation
    bit_str = model.basis_states[:,n]
    bit_int = model.basis_ints[:,n]

    # creating the output basis states and weights
    output_n      = Array{Int,1}(undef, 0)
    output_weight = Array{Float64,1}(undef, 0)

    # perform the diagonla action of the Hamiltonian
    # Interaction -> we need to measure the number of neigbouring a&b in the string
    # the string is organized as: a1 b1 a2 b2 a3 b3 ...
    # so we need to multiply the addition of one pack by the next pack and add
    # all these contributions together
    elem = 0
    for jj = 1:2:model.totalSites-2
        tmp1 = bit_int[jj] + bit_int[jj+1]
        tmp2 = bit_int[jj+2] + bit_int[jj+3]
        elem = elem + tmp1*tmp2
    end
    if(elem>0)
        push!(output_n, n)
        push!(output_weight, elem*model.u)
    end

    # inner unit hubbard interaction
    n_doubl =  0
    for jj = 1:model.totalSites -1
        if((bit_int[jj] & bit_int[jj+1]) == 1)
            n_doubl = n_doubl + 1
        end
        if(n_doubl > 0)
            push!(output_n, model.basis_dict[join(bit_str)])
            push!(output_weight, n_doubl*model.w2)
        end
    end

    # offdiagonal part.
    for jj = 1:model.totalSites - 2
        # hopping part
        # see if hopping will happens

        if((bit_int[jj] ⊻ bit_int[jj+2]) == 1)
            work_str = copy(bit_str)
            # calcultate the sign wich will be aquired after hopping
            # sign by annihilation
            # for the hopping, we only need to to take an possible
            # additional a(b) particle which may sit in between into account
            # this will fix the sign completly
            #sign1 = mod(sum(bit_int[1:jj-1]),2) == 1 ? -1 : 1
            #sign2 = bit_int[jj+1] == 0 ? sign1 : -1*sign1
            sign = bit_int[jj+1] == 0 ? 1 : -1

            # now we need to perform the hopping
            work_str[jj]   = bit_str[jj+2]
            work_str[jj+2] = bit_str[jj]
            push!(output_n, model.basis_dict[join(work_str)])
            push!(output_weight,-1*sign*model.t)
        end


        # only select the terms where jj represents a 'a' particle
        if((jj < model.totalSites - 2) && mod(jj,2) == 1)
            # this selects all possible configurations for the pair hopping terms
            # like ... 1 0 1 0 .... or ... 1 0 0 1 ...
            if(((bit_int[jj] ⊻ bit_int[jj+1]) & (bit_int[jj+2] ⊻ bit_int[jj+3])) == 1)
                work_str = copy(bit_str)
                # now decide wether this is a W1 or W2 configuration
                # W1 configuration have an alternative pattern 1 0 1 0 or 0 1 0 1
                # W2 a cluster 1 0 0 1 or 0 1 1 0
                if(bit_int[jj+1] ⊻ bit_int[jj+2] == 1)
                    work_str[jj]   = bit_str[jj+1]
                    work_str[jj+1] = bit_str[jj]
                    work_str[jj+2] = bit_str[jj+3]
                    work_str[jj+3] = bit_str[jj+2]

                    # we dont need a sign here
                    push!(output_weight, -1*model.w1)
                    push!(output_n, model.basis_dict[join(work_str)])
                else
                    # W2 configuration#
                    #work_str[jj]   = bit_str[jj+2]
                    #work_str[jj+2] = bit_str[jj]
                    #work_str[jj+1] = bit_str[jj+3]
                    #work_str[jj+3] = bit_str[jj+1]

                    # we dont need a sign here
                    #push!(output_weight, -1*model.w2)
                end
            end
        end
    end

    return output_n, output_weight
end

precompile(applyHam,(Int, LocalSpace, ))

function buildHam(model::LocalSpace)
        @printf "Start building the Hamiltonian for L=%.0f and N=%.0f\n" model.L model.N
        @printf "Parameter for %s model: t=%.1f, W1=%.1f, W2=%.1f, U=%.1f\n" model.name model.t model.w1 model.w2 model.u

        ham = zeros(Float64, model.dimHilb, model.dimHilb)

        for jj = 1:model.dimHilb
            out_n, out_weight = applyHam(jj, model)
            for kk = 1:size(out_n,1)
                ham[jj,out_n[kk]] = out_weight[kk]
            end
        end

        return ham
end

precompile(buildHam, (LocalSpace,))

function measure_n(vec::Array{Float64,1}, model::LocalSpace; n::Int = -1)

    out = zeros(Float64, model.L)
    if(n == -1)
        n_   = 0
        mult = 2
    elseif(n == 1)
        n_ = 0
        mult = 1
    elseif(n == 2)
        n_ = 1
        mult = 1
    else
        return out
    end
    for jj = 1:model.L
        for uu = 1:mult
            for kk = 1:model.dimHilb
                if model.bit_ints[(jj-1)*2 + uu + n_,kk] == 1
                    out[jj] = out[jj] + abs(vec[kk])^2
                end
            end
        end
    end

    return out
end

precompile(measure_n, (Array{Float64,1}, LocalSpace, Int, ))

function measure_double(vec::Array{Float64,1}, model::LocalSpace)
    out = zeros(Float64, model.L)

    for jj = 1:model.L
        for kk = 1:model.dimHilb
            if ((model.basis_ints[(jj-1)*2 + 1, kk] & model.basis_ints[(jj-1)*2 + 2, kk]) == 1)
                out[jj] = out[jj] + abs(vec[kk])^2
            end
        end
    end

    return out
end
precompile(measure_n, (Array{Float64,1}, LocalSpace,))

end # module

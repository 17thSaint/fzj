module MajoranaLadder

using Combinatorics
using Printf, LinearAlgebra, Arpack
using SparseArrays
using ProgressMeter

struct LocalSpace
    L::Int64 # physical length of the system
    N::Int64 # number of electrons in the system

    dimHilb :: Int64 # total dimension of hilbertspace
    totalSites :: Int64 # Size with enrolled unit cell

    basis_states :: Array{String,2} # basis states in string rep
    basis_ints   :: Array{Int,2}    # basis states in integer rep
    basis_dict   :: Dict{String,Int} # fast access to the position of basis in array

    t_a    :: Real
    t_b    :: Real
    mu     :: Real
    w1     :: Real
    w2     :: Real
    u      :: Real
    ualpha :: Real

    name :: String
end


function LocalSpace_init(L::Int64,N::Int64, paraDict :: Dict{String, T} where T<:Real; use_z2::Bool = true, z_sector::Int64 = 0)

    bit_strings, bit_ints, basis_dict = generate_basis(L,N, use_z2 = use_z2, z_sector = z_sector)

    dimHilb = size(bit_strings,2)

    t      = get(paraDict, "t", 1.0)
    t_a    = get(paraDict, "ta", t)
    t_b    = get(paraDict, "tb", t)
    mu     = get(paraDict, "mu", 0.0)
    w1     = get(paraDict, "w1", 0.0)
    w2     = get(paraDict, "w2", w1)
    u      = get(paraDict, "u", 0.0)
    ualpha = get(paraDict, "ualpha", 0.0)

    majLad = LocalSpace(L,N,dimHilb, L*2, bit_strings, bit_ints, basis_dict, t_a, t_b, mu, w1, w2, u, ualpha,"MajoranaLadder")
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
    @showprogress for subs in subsets
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
        @showprogress for jj = 1:n_states
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

    # perform the diagonal action of the Hamiltonian
    # Interaction -> we need to measure the number of neigbouring a&b in the string
    # the string is organized as: a1 b1 a2 b2 a3 b3 ...
    # so we need to multiply the addition of one pack by the next pack and add
    # all these contributions together
    elem1 = 0
    elem2 = 0 # intra wire interaction, assumed to be the same for a and b
    chm_pot = 0 # chemical potential, energy offest
    for jj = 1:2:model.totalSites-2
        tmp1 = bit_int[jj] + bit_int[jj+1]
        tmp2 = bit_int[jj+2] + bit_int[jj+3]
        elem1 = elem1 + tmp1*tmp2

        tmp1 = bit_int[jj] * bit_int[jj+2]
        tmp2 = bit_int[jj+1] * bit_int[jj+3]
        elem2 = elem2 + tmp1 + tmp2
    end
    # chemical potential is just the sum over all integers
    chm_pot = sum(bit_int)*model.mu
    elem = model.u*elem1 + model.ualpha*elem2 + chm_pot
    if(elem>0)
        push!(output_n, n)
        push!(output_weight, elem)
    end


    # offdiagonal part.
    for jj = 1:model.totalSites - 2
        # hopping part
        # see if hopping will happens
        if((bit_int[jj] ⊻ bit_int[jj+2]) == 1)
            # choosing the correct hopping amplitude
            t = mod(jj,2) == 1 ? model.t_a : model.t_b
            # odd particles describing the upper chain (A), even
            # one the lower chain(B)
            work_str = copy(bit_str)
            # calcultate the sign wich will be aquired after hopping
            # sign by annihilation
            # for the hopping, we only need to to take an possible
            # additional a(b) particle which may sit in between into account
            # this will fix the sign completly
            #sign1 = mod(sum(bit_int[1:jj-1]),2) == 1 ? -1 : 1
            #sign2 = bit_int[jj+1] == 0 ? sign1 : -1*sign1
            sign = iszero(bit_int[jj+1]) ? 1 : -1

            # now we need to perform the hopping
            work_str[jj]   = bit_str[jj+2]
            work_str[jj+2] = bit_str[jj]
            push!(output_n, model.basis_dict[join(work_str)])
            push!(output_weight,-1*sign*t)
        end


        # only select the terms where jj represents a 'a' particle
        if((jj < model.totalSites - 2) && isodd(jj))
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
                else
                    # W2 configuration#
                    work_str[jj]   = bit_str[jj+2]
                    work_str[jj+2] = bit_str[jj]
                    work_str[jj+1] = bit_str[jj+3]
                    work_str[jj+3] = bit_str[jj+1]

                    # we dont need a sign here
                    push!(output_weight, -1*model.w2)
                end
                push!(output_n, model.basis_dict[join(work_str)])
            end

        end
    end

    return output_n, output_weight
end

precompile(applyHam,(Int, LocalSpace, ))

function buildHam(model::LocalSpace)
        @printf "Start building the Hamiltonian for L=%.0f and N=%.0f\n" model.L model.N
        @printf "Parameter for %s model: t_a=%.1f, W1=%.1f, W2=%.1f, U=%.1f, Ualpha = %.1f\n" model.name model.t_a model.w1 model.w2 model.u model.ualpha

        ham = zeros(Float64, model.dimHilb, model.dimHilb)
        p = Progress(model.dimHilb,1) #progress meter

        for jj = 1:model.dimHilb
            out_n, out_weight = applyHam(jj, model)
            for kk = 1:size(out_n,1)
                ham[jj,out_n[kk]] += out_weight[kk]
            end
            next!(p) # update meter
        end

        return ham
end

precompile(buildHam, (LocalSpace,))

function buildHam_sp(model::LocalSpace)
    @printf "Start building the Hamiltonian for L=%.0f and N=%.0f\n" model.L model.N
    @printf "Parameter for %s model: t=%.1f, W1=%.1f, W2=%.1f, U=%.1f\n" model.name model.ta model.w1 model.w2 model.u

    ham = spzeros(model.dimHilb,model.dimHilb)

    p = Progress(model.dimHilb,1) # progress meter
    for jj = 1:model.dimHilb
        out_n, out_weight = applyHam(jj, model)

        for kk = 1:size(out_n,1)
            ham[jj,out_n[kk]] += out_weight[kk]
        end
        next!(p) # update meter
    end

    return ham
end
precompile(buildHam_sp, (LocalSpace,))

function buildHam_sp_pa(model::LocalSpace)
    @printf "Start building the Hamiltonian for L=%.0f and N=%.0f\n" model.L model.N
    @printf "Parameter for %s model: ta=%.1f, tb = %.1f, W1=%.1f, W2=%.1f, U=%.1f\n" model.name model.t_a model.t_b model.w1 model.w2 model.u

    ham = spzeros(model.dimHilb,model.dimHilb)

    nThreads = Threads.nthreads() # local save of number of threads

    # generating a sparse array for every thread
    sparse_arrays = Array{SparseMatrixCSC{Float64,Int64}}(undef,nThreads)

    # initialize them to be zero
    for jj = 1:nThreads
        sparse_arrays[jj] = spzeros(model.dimHilb,model.dimHilb)
    end

    p = Progress(model.dimHilb,1) # progress meter
    Threads.@threads for jj = 1:model.dimHilb
        out_n, out_weight = applyHam(jj, model)

        for kk = 1:size(out_n,1)
            sparse_arrays[Threads.threadid()][jj,out_n[kk]] += out_weight[kk]
            #ham[jj,out_n[kk]] += out_weight[kk]
        end
        next!(p) # update meter
    end

    # no collapse the sparse arrays to the final Hamiltonian
    for jj = 1:nThreads
        ham += sparse_arrays[jj]
    end

    return ham
end
precompile(buildHam_sp_pa, (LocalSpace,))

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

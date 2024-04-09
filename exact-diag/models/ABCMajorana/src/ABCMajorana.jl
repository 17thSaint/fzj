module ABCMajorana

using Combinatorics
using Printf, LinearAlgebra, Arpack
using SparseArrays
using ProgressMeter


struct LocalSpace
    L::Int64 # physical length of the system
    N::Int64 # number of electrons in the system
    use_z2::Bool # z2 symmetry was used to reduce the hilbertspace
    z_sector::Int64 # Z2 symmetry sector

    dimHilb :: Int64 # total dimension of hilbertspace
    totalSites :: Int64 # Size with enrolled unit cell

    basis_states :: Array{String,2} # basis states in string rep
    basis_ints   :: Array{Int,2}    # basis states in integer rep
    basis_dict   :: Dict{String,Int} # fast access to the position of basis in array

    t      :: Float64
    J      :: Float64
    Vp     :: Float64
    Vm     :: Float64
    mu     :: Float64
    ualpha :: Float64

    name :: String
end


function LocalSpace_init(L::Int64,N::Int64, paraDict :: Dict{String, Float64}; use_z2::Bool = true, z_sector::Int64 = 0)

    bit_strings, bit_ints, basis_dict = generate_basis(L,N, use_z2 = use_z2, z_sector = z_sector)

    dimHilb = size(bit_strings,2)

    t      = get(paraDict, "t", 1.0)
    J      = get(paraDict, "J", 1.0)
    V1     = get(paraDict, "V1", 0.0)
    V2     = get(paraDict, "V2", 0.0)
    ualpha = get(paraDict, "ualpha", 0.0)
    mu     = get(paraDict, "mu", 0.0)

    Vp = (V1 + V2)/2
    Vm = (V1 - V2)/2

    abcMaj = LocalSpace(L,N,use_z2, z_sector,dimHilb, L*4, bit_strings, bit_ints, basis_dict, t, sqrt(2)*J, Vp, Vm, mu, ualpha, "ABCMajorana")
    return abcMaj
end

precompile(LocalSpace_init,(Int,Int,Bool, Int))
function generate_basis(L::Int64,N::Int64; use_z2::Bool = true, z_sector::Int64 = 0)

    @printf "Strart to calulate basis states for L=%.0f and N=%.0f\n" L N
    @printf "Use reduction of hilbertspace by Z2 symmetry: %s\n" use_z2

    z_sector_ = mod(z_sector,2)

    n_unit = 4 # unit cell size
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
        # we have the convention  .... a m p b .... and the Z2 symmetry
        # is related to the parity of a+m
        n_states_new = 0
        match_ind    = Array{Int,1}(undef,0)
        @showprogress for jj = 1:n_states
            sum = 0
            for kk = 1:4:totalSites
                sum = sum + bit_ints[kk,jj]
                sum = sum + bit_ints[kk + 1, jj]
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
    # the string is organized as: a1 m1 p1 b1 a2 m2 p2 b2 a3 m3 p3 b3 ...
    # so we need to multiply the addition of one pack by the next pack and add
    # all these contributions together
    elem1 = 0
    elem2 = 0
    elem3 = 0
    for jj = 1:4:model.totalSites
        # chemical potential part
        elem1 += bit_int[jj+1] + bit_int[jj+2]
        # interaction part
        if(jj < model.totalSites - 3)
            # a-a and b-b interaction
            tmp1 = bit_int[jj] * bit_int[jj+4] # a part.
            tmp2 = bit_int[jj+3] * bit_int[jj+7] # b part.
            elem2 += tmp1 + tmp2

            # total unit cell p-m interaction
            tmp1 = bit_int[jj+1] + bit_int[jj+2]
            tmp2 = bit_int[jj+5] + bit_int[jj+6]
            elem3 += tmp1*tmp2
        end

    end
    elem = model.mu * elem1 + model.ualpha * elem2 + model.Vp*elem3
    if(abs(elem)>0)
        push!(output_n,n)
        push!(output_weight, elem)
    end

    # ABC Clock part
    for jj = 1:model.totalSites
        if(mod(jj,4) == 1 || mod(jj,4) == 3)
            if((bit_int[jj] ⊻ bit_int[jj+1]) == 1)
                work_str = copy(bit_str)
                work_str[jj]   = bit_str[jj+1]
                work_str[jj+1] = bit_str[jj]

                push!(output_n, model.basis_dict[join(work_str)])
                push!(output_weight,-1*model.J)
            end
        end
    end

    # offdiagonal part.
    for jj = 1:model.totalSites - 4
        # hopping part
        # see if hopping will happens

        # only a and b particles may hop, so pick the corresponding jj
        # jj mod 4 == 1 -> a particle
        # jj mod 4 == 0 -> b particle
        if(mod(jj,4) == 1 || mod(jj,4) == 0)
            if((bit_int[jj] ⊻ bit_int[jj+4]) == 1)
                work_str = copy(bit_str)
                # calcultate the sign wich will be aquired after hopping
                # sign by annihilation
                # for the hopping, we only need to to take an possible
                # additional a(b) particle which may sit in between into account
                # this will fix the sign completly
                #sign1 = mod(sum(bit_int[1:jj-1]),2) == 1 ? -1 : 1
                #sign2 = bit_int[jj+1] == 0 ? sign1 : -1*sign1

                # JW string for hopping, counting all particles in between
                sign = mod(sum(bit_int[jj+1:jj+3]),2) == 0 ? 1 : -1

                # now we need to perform the hopping
                work_str[jj]   = bit_str[jj+4]
                work_str[jj+4] = bit_str[jj]
                push!(output_n, model.basis_dict[join(work_str)])
                push!(output_weight,-1*sign*model.t)
            end
        end



        # now we need the pairhopping part, for this select the terms where
        # jj represents a 'm' particle, i.e. mod(jj,4) == 2
        if(mod(jj,4) == 2)
            # this selects all possible configurations for the pair hopping terms
            # like ... na 1 0 nb na 1 0 nb .... or ... na 1 0 nb na 0 1 nb ...
            if(((bit_int[jj] ⊻ bit_int[jj+1]) & (bit_int[jj+4] ⊻ bit_int[jj+5])) == 1)
                work_str = copy(bit_str)
                # now decide wether this is a W1 or W2 configuration
                # W1 configuration have an alternative pattern na 1 0 nb na 1 0 nb or na 0 1 nb na 0 1 nb
                # W2 a cluster na 1 0 nb na 0 1 nb or na 0 1 nb na 1 0 nb
                if(bit_int[jj+1] ⊻ bit_int[jj+4] == 1)
                    work_str[jj]   = bit_str[jj+1]
                    work_str[jj+1] = bit_str[jj]
                    work_str[jj+4] = bit_str[jj+5]
                    work_str[jj+5] = bit_str[jj+4]

                    # we dont need a sign here
                    push!(output_weight, model.Vm)
                else
                    # W2 configuration#
                    work_str[jj]   = bit_str[jj+4]
                    work_str[jj+4] = bit_str[jj]
                    work_str[jj+1] = bit_str[jj+5]
                    work_str[jj+5] = bit_str[jj+1]

                    # we dont need a sign here
                    push!(output_weight, model.Vm)
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
        #@printf "Parameter for %s model: t=%.1f, W1=%.1f, W2=%.1f, U=%.1f, Ualpha = %.1f\n" model.name model.t model.w1 model.w2 model.u model.ualpha

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
    #@printf "Parameter for %s model: t=%.1f, W1=%.1f, W2=%.1f, U=%.1f\n" model.name model.t model.w1 model.w2 model.u

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
    #@printf "Parameter for %s model: t=%.1f, W1=%.1f, W2=%.1f, U=%.1f\n" model.name model.t model.w1 model.w2 model.u

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

function measure_n(vec::Array{T,1} where T<:Number, model::LocalSpace; n::String = "tot")

    out = zeros(Float64, model.L)
    if(n == "tot")
        n_   = 0
        mult = 4
    elseif(n == "a")
        n_ = 0
        mult = 1
    elseif(n == "m")
        n_ = 1
        mult = 1
    elseif(n == "p")
        n_ = 2
        mult = 1
    elseif(n == "b")
        n_ = 3
        mult = 1
    else
        return out
    end
    for jj = 1:model.L
        for uu = 1:mult
            for kk = 1:model.dimHilb
                if model.basis_ints[(jj-1)*4 + uu + n_,kk] == 1
                    out[jj] = out[jj] + abs(vec[kk])^2
                end
            end
        end
    end

    return out
end

precompile(measure_n, (Array{Float64,1}, LocalSpace, Int, ))
precompile(measure_n, (Array{Complex{Float64},1}, LocalSpace, Int, ))

# Measures the pair correlation between u_1 electron at position j1 and u2 electron at position j2. I.e:
# C(u1,j1; u2, j2) := <c_{u1,j1} c_{u2,j2}^\dagger>
function measure_pcor_single(vec::Array{T,1} where T <: Number, model::LocalSpace, u1 :: String, j1 :: Int64, u2 :: String, j2 :: Int64)

    @assert j2 <= model.L "Position j2 has to be smaller than the system size."
    @assert j1 > 0 "Position j1 has to be larger than 0."

    n_unitcell = 4

    out = 0

    is_complex = typeof(vec[1]) <: Complex

    # mapping of strings to integers representing the internal counting of the unitcell
    u1_int = 0
    if(u1 == "a")
        u1_int = 1
    elseif(u1 == "m")
        u1_int = 2
    elseif(u1 == "p")
        u1_int = 3
    elseif(u1 == "b")
        u1_int = 4
    else
        return out
    end

    u2_int = 0
    if(u2 == "a")
        u2_int = 1
    elseif(u2 == "m")
        u2_int = 2
    elseif(u2 == "p")
        u2_int = 3
    elseif(u2 == "b")
        u2_int = 4
    else
        return out
    end

    pos1 = (j1 - 1)*n_unitcell + u1_int
    pos2 = (j2 - 1)*n_unitcell + u2_int

    onpoint = pos1 == pos2

    for (jj,vk) in enumerate(vec)

        # load the bit representation
        bit_int = model.basis_ints[:,jj]
        if(onpoint)
            if(bit_int[pos1] == 1)
                out += abs(vk)^2
            end
            continue
        end
        if(bit_int[pos1] == 1 && bit_int[pos2] == 0)
            # load the string representation
            bit_str = model.basis_states[:,jj]
            # performing the hopping
            bit_str[pos1] = "0"
            bit_str[pos2] = "1"

            # get the position inside the basis states of this string
            # if not present in the basis state, simply go to the next state
            # since in this case the state is not available in the symmetry reduced
            # basis. -> TODO better catching of this case....
            hh = 0
            try
                hh = model.basis_dict[join(bit_str)]
            catch
                continue
            end
            # the sign is the parity of all particles sitting between pos1 and pos2
            sign = mod(sum(bit_int[pos1 + 1:pos2 - 1]),2) == 0 ? 1 : -1
            if(is_complex)
                out += sign*conjugate(vec[hh]) * vk
            else
                out += sign*vec[hh] * vk
            end
        end
    end
    return out
end

precompile(measure_pcor_single, (Array{Float64,1}, LocalSpace, String, Int64, String, Int64, ))
precompile(measure_pcor_single, (Array{Complex{Float64},1}, LocalSpace, String, Int64, String, Int64, ))

# measures the pair correlation between j1 and k for k between j1 and j2
function measure_pair_correlation(vec::Array{T,1} where T <: Number, model::LocalSpace, u1 :: String, j1 :: Int64, u2 :: String, j2 :: Int64)
    @assert j2 >= j1 "Position j2 has to be larger than position j1."
    out = zeros(Float64, j2 - j1 + 1)

    for jj = j1:j2
        out[jj - j1 + 1] = measure_pcor_single(vec, model, u1 , j1 , u2, jj)
    end
    return out
end

precompile(measure_pair_correlation, (Array{Float64,1}, LocalSpace, String, Int64, String, Int64, ))
precompile(measure_pair_correlation, (Array{Complex{Float64},1}, LocalSpace, String, Int64, String, Int64, ))

function initialize_product_wf(model::LocalSpace,state_arr::Array{String,1})
    @assert model.L == size(state_arr,1) "Size of given state array not compatible with system length!"

    # converting the state array into the bits

    # getting an empty integer bit array
    bit_array = Array{String,1}(undef, model.totalSites)
    fill!(bit_array,"0")

    # now identify the positions in the bit_array with the given states
    # '0' -> empty state
    # '1' - '4' -> single state
    # '12', '13', .. -> double occupied state
    # '123', '124', ... -> three electron state
    # '1234' -> fully occupied

    # initialize the symmetry observers
    n_tot = 0
    z_sec = 0
    for jj = 1:model.L
        @assert 0 < length(state_arr[jj]) < 5 "Unknown state given: "*state_arr[jj]

        # splitting the state into its components
        split_state = split(state_arr[jj],"")
        # handle empty state specially
        if (size(split_state,1) == 1 && split_state[1] == "0")
            continue
        end

        # adding the number of particles == size(split_state,1) to
        # the observer
        n_tot += size(split_state,1)
        # now fill the bit_string prepared with the correct elements.
        # for this, first calculate the offset comming from the position
        # in the chain
        offset = (jj-1)*4

        # now loop through all character bins,
        # try to convert them to integers,
        # add this integer to the offset and overwrite the 0 in the
        # bit_arr to '1', on the fly, update the parity observer
        # '1' and '2' representing 'a' and 'm' states -> they contribute
        # to the parity
        for jj in 1:size(split_state,1)
            bit = split_state[jj]
            @assert !in(bit, split_state[1:jj-1]) "Bit appearces twice."
            if (bit == "1" || bit == "2")
                z_sec += 1
            end
            bit_int = 0
            try
                bit_int = parse(Int64,bit)
            catch
                throw(DomainError(bit, "Bits have to be encoded numbers."))
            end
            @assert 0<bit_int<5 "Bits have to be in the range from 1 to 4."
            bit_array[bit_int + offset] = "1"
        end
    end

    # now checkt the overall symmetry sectors
    @assert n_tot == model.N "Given state has a different number of particles as the given model!"
    if (model.use_z2)
        @assert mod(z_sec,2) == model.z_sector "Given state has a different Z2 as the given model!"
    end

    # now initialize an empty weight vector
    weights = Array{Float64,1}(undef, model.dimHilb)
    fill!(weights,0)

    # and set the position corresponding to the given state to 1
    weights[model.basis_dict[join(bit_array)]] = 1
    return weights
end

precompile(initialize_product_wf, (LocalSpace, Array{string,1}, ))

end # module

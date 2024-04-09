module XXZ

using Combinatorics
using Printf, LinearAlgebra, Arpack

function generate_basis(L::Int64,N::Int64)

    n_states = binomial(L,N)

    bit_strings = Array{String,2}(undef,n_states,L)
    bit_ints    = Array{Int,2}(undef,n_states,L)

    basis_dict  = Dict()

    subsets = combinations([1:1:L;],N)

    j = 1
    for subs in subsets
        bit_strings[j,:] = fill("0",L)
        bit_ints[j,:]     = fill(0,L)
        for en in subs
            bit_strings[j,en] = "1"
            bit_ints[j,en]     = 1
        end
        basis_dict[join(bit_strings[j,:])] = j
        j = j + 1
    end
    return bit_strings, bit_ints, basis_dict
end

precompile(generate_basis,(Int,Int,))

function applyHam(n::Int64, basis::Array{String,2}, int_basis::Array{Int,2},
                  basis_dict::Dict, J::Float64, Delta::Float64, H::Float64)
    #J     = para[1]
    #Delta = para[2]
    #H     = para[3] #-> total magnetic field, precomputed, since sum_{j=1}^L S^Z is fixed
    L     = size(basis,2)
    #show(stdout, "text/plain", basis)
    bit_str = basis[n,:]
    bit_int = int_basis[n,:]

    # creating the output arrays holding the final basis numbers
    # and the corresponding weights
    output_n     = Array{Int,1}(undef,0)
    output_weight = Array{Float64,1}(undef,0)

    # performing the diagonal part of the hamiltonian
    # i.e. the ZZ compontent
    elem = 0.0
    for jj = 1:L-1
        elem = elem + (bit_int[jj] - 0.5)*(bit_int[jj+1] - 0.5)
    end
    push!(output_n,n)
    push!(output_weight,elem*Delta + H)

    # off diagonal part
    for jj = 1:L-1
        if((bit_int[jj] ⊻ bit_int[jj+1]) == 1 ) # performing a xor to see if flip happens
            work_str = copy(bit_str)
            work_str[jj]   = bit_str[jj+1]
            work_str[jj+1] = bit_str[jj]
            push!(output_n, basis_dict[join(work_str)])
            push!(output_weight,J/2)
        end
    end

    return output_n, output_weight
end
precompile(applyHam,(Array{String,2}, Array{Int,2}, Dict, Float64, Float64, Float64, ))

function buildHam(L::Int64, N::Int64, para::Array{Float64,1})
    J     = para[1]
    Delta = para[2]
    h     = para[3]

    @printf "Start building the Hamiltonian for L=%.0f and N=%.0f\n" L N
    @printf "Parameter for XXZ model: J=%.1f, Delta=%.1F, h=%.1f\n" J Delta h

    bit_strings, bit_ints, bit_dict = generate_basis(L,N)

    dim = size(bit_strings,1)
    @printf "Size of local many body hilbertspace: %.0f\n" dim
    # the total magnetic energy is a constant throughout the basis sets, precompute it
    H = h*sum(bit_ints[1,:] - fill(0.5,L))

    # now prepare the hamilton matrix
    ham = zeros(Float64, dim, dim)

    for jj = 1:dim
        out_n, out_weight = applyHam(jj, bit_strings, bit_ints, bit_dict,
                                     J, Delta, H)
        for kk = 1:size(out_n,1)
            ham[jj,out_n[kk]] = out_weight[kk]
        end
    end

    return ham, bit_strings, bit_ints, bit_dict
end
precompile(buildHam, (Int, Int, Array{Float64,1}))

function measure_sZ(vec::Array{Float64,1},bit_ints::Array{Int,2})
    L = size(bit_ints,2)
    dim = size(bit_ints,1)
    out = zeros(Float64,L)

    for jj = 1:dim
        for kk = 1:L
            out[kk] = out[kk] + abs(vec[jj])^2 * (bit_ints[jj,kk] - 0.5)
        end
    end
    return out
end
precompile(measure_sZ, (Array{Float64,1}, Array{Int,2},))

function measure_sZ(vec::Array{Complex,1},bit_ints::Array{Int,2})
    L = size(bit_ints,2)
    dim = size(bit_ints,1)
    out = zeros(Float64,L)

    for jj = 1:dim
        for kk = 1:L
            out[kk] = out[kk] + abs(vec[jj])^2 * (bit_ints[jj,kk] - 0.5)
        end
    end
    return out
end
precompile(measure_sZ, (Array{Complex,1}, Array{Int,2},))

function quench_sp(vec::Array{Float64,1}, pos::Int64,
                   basis::Array{String,2}, int_basis::Array{Int,2})

    dim_old = size(vec,1)
    L       = size(basis,2)
    N       = sum(int_basis[1,:])

    # if the system is completly up, return a 0 dim. vector
    if(N == L)
        return []
    end
    # build the basis of the hilbertspace with one up added
    bit_strings_en, bit_ints_en, bit_dict_en = generate_basis(L,N+1)
    dim_new = size(bit_strings_en,1)

    out_vec = zeros(Float64,dim_new)

    for jj =1:dim_old
        if(int_basis[jj,pos] == 0)
            work_arr = copy(basis[jj,:])
            work_arr[pos] = "1"
            pos_new = bit_dict_en[join(work_arr)]
            out_vec[pos_new] = vec[jj]
        end
    end
    return out_vec/norm(out_vec)
end

function get_weights(vec::Array{Float64,1}, eig_vec::Array{Float64,2})
    weights = transpose(eig_vec)*vec
    return weights
end
precompile(get_weights, (Array{Float64,1}, Array{Float64,2}, ))

function get_real_weigths(weights::Array{Complex,1},eig_vec::Array{Float64,2})
    real_weights = eig_vec*weights
    return real_weights
end
precompile(get_real_weigths, (Array{Complex,1}, Array{Float64,2},))

function time_measurement(L::Int64, N::Int64, para::Array{Float64,1}, pos::Array{Int64,1},
                          dt::Float64, timesteps::Int64)

    t1 = time()
    # first build the Hamiltonian without quenches
    ham, basis, int_basis, basis_dict = buildHam(L,N, para)
    t2 = time()

    @printf "Finished building hamiltonian. Time needed: %f\n" t2 - t1


    # calculate the eigenstates/energies
    @printf "Calculating the eigensystem.\n"
    Eigsys = eigen(ham)
    t3 = time()

    @printf "Finished. Time needed: %f\n" t3 - t2

    vals = Eigsys.values
    vecs = Eigsys.vectors

    @printf "\nEnergy of GS: %.0f\n" vals[1]
    if(size(vals,1)>1)
        @printf "First energy gap: %.0f\n" vals[2] - vals[1]
    end

    # getting the groundstate in this sector:
    gs = vecs[:,1]



    @printf "\n\nPerforming quenches on gs.\n"
    quenched_gs = copy(gs)
    for jj = 1:size(pos,1)
        quenched_gs = quench_sp(quenched_gs, pos[jj], basis, int_basis)
        basis, int_basis, basis_dict = generate_basis(L,N+jj)
    end
    t4 = time()

    @printf "Finished. Time needed: %f\n" t4 - t3
    ham, basis, int_basis, basis_dict = buildHam(L,N+size(pos,1), para)
    t5 = time()
    @printf "Finished building hamiltonian. Time needed: %f\n" t5 - t4


    @printf "Calculating the eigensystem.\n"
    Eigsys = eigen(ham)
    t6 = time()
    @printf "Finished. Time needed: %f\n" t6 - t5

    vals = Eigsys.values
    vecs = Eigsys.vectors


    # getting the weights on the eigenvectors for the quenched state
    weights = get_weights(quenched_gs, vecs)


    # opening outstream
    out_stream = open("local_observables.dat", "w")


    #s_z = zeros(Complex,L)
    # printing the first line
    write(out_stream, "timestep\tsite\tReal(s_z)\tIm(s_z)\n")

    # printing the first measurement of sz
    s_z = measure_sZ(quenched_gs,int_basis)
    for ll = 1:L
        write(out_stream,"0\t",string(ll),"\t",string(real(s_z[ll])),"\t",string(imag(s_z[ll])),"\n")
    end

    # now multiply the phases at the different stages and rotate back to measure the
    # sz:
    @printf "Start Time measurments.\n"
    for tt = 1:timesteps
        # getting the complex phases
        phases = map(E->exp(-dt*tt*E*im), vals)
        # now get the rotated weights
        rot_weights = zeros(Complex, size(phases,1))
        for jj = 1:size(phases,1)
            rot_weights[jj] = phases[jj]*weights[jj]
        end
        # get the time evolved vector
        vec_t = vecs*rot_weights
        # measure the time evolved local magnetization
        s_z = measure_sZ(vec_t,int_basis)
        for ll = 1:L
            write(out_stream,string(tt),"\t",string(ll),"\t",string(real(s_z[ll])),"\t",string(imag(s_z[ll])),"\n")
        end
    end
    t7 = time()
    @printf "Finished. Time Needed: %f\n" t7-t6

    close(out_stream)
end
precompile(time_measurement, (Int, Int, Array{Float64,1}, Array{Int,1}, Float64, Int, ))


end # module

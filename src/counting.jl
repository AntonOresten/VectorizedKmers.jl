"""
    count_kmers!(kmer_array, sequence; reset=true)
"""
function count_kmers! end

# primarily meant to be used for testing
function count_kmers!(
    kmer_array::KmerArray{N, K}, sequence::Vector{Int};
    reset::Bool = true,
) where {N, K}
    @assert 0 <= maximum(sequence) < N
    reset && zeros!(kmer_array)
    mask = N^K
    kmer = zero(Int)
    for (i, m) in enumerate(sequence)
        kmer = kmer * N % mask + m
        K <= i && (kmer_array[kmer] += 1)
    end
    return kmer_array
end

"""
    count_kmers(sequence, N, K, T=Int, zeros=zeros)
"""
function count_kmers end

function count_kmers(sequence, N::Integer, K::Integer, T::Type{<:Real}=Int, zeros=zeros)
    return count_kmers!(KmerArray(N, K, T, zeros), sequence; reset=false)
end

alphabet_size(T::Type) = error("$(T) does not have a defined alphabet size. Please define `alphabet_size(::Type{<:$(T)})` or insert the alphabet size as a second argument in the `count_kmers` function call.")
alphabet_size(::T) where T = alphabet_size(T)

function count_kmers(sequence, K::Integer, T::Type{<:Real}=Int, zeros=zeros)
    return count_kmers(sequence, alphabet_size(sequence), K, T, zeros)
end
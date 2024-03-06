"""
    count_kmers!(kmer_array, seq; reset=true)
"""
function count_kmers! end

# primarily meant for testing purposes
function count_kmers!(
    kmer_array::KmerArray{N, K}, seq::Vector{Int};
    reset::Bool = true,
) where {N, K}
    0 <= maximum(seq) < N || throw(ArgumentError("seq contains elements outside the range 0:N-1"))
    reset && zeros!(kmer_array)
    mask = N^K
    kmer = zero(Int)
    for (i, m) in enumerate(seq)
        kmer = kmer * N % mask + m
        K <= i && (kmer_array[kmer] += 1)
    end
    return kmer_array
end

"""
    count_kmers(seq, N, K, T=Int, zeros=zeros)
"""
function count_kmers end

count_kmers(seq, N::Integer, K::Integer, T::Type{<:Real}=Int, zeros=zeros) = count_kmers!(KmerArray(N, K, T, zeros), seq; reset=false)

alphabet_size(T::Type) = error("$(T) does not have a defined alphabet size. Please define `alphabet_size(::Type{<:$(T)})` or insert the alphabet size as a second argument in the `count_kmers` function call.")
alphabet_size(::T) where T = alphabet_size(T)

count_kmers(seq, K::Integer, T::Type{<:Real}=Int, zeros=zeros) = count_kmers(seq, alphabet_size(seq), K, T, zeros)
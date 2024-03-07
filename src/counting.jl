"""
    count_kmers!(kmer_array, seq; reset=true)
"""
function count_kmers! end

# primarily meant for testing purposes
function count_kmers!(
    kmer_array::KmerArray{N, K}, seq;
    reset::Bool = true,
) where {N, K}
    reset && zeros!(kmer_array)
    mask = N^K
    kmer = zero(UInt)
    for (i, m) in enumerate(seq)
        kmer = kmer * N % mask + axis_index(kmer_array, m)
        kmer_array[kmer] += K <= i
    end
    return kmer_array
end

"""
    count_kmers(seq, N, K, T=Int, zeros=zeros)
"""
function count_kmers end

count_kmers(seq, N::Integer, K::Integer, T::Type{<:Real}=Int, zeros=zeros) = count_kmers!(KmerArray(N, K, T, zeros), seq; reset=false)
count_kmers(seq, K::Integer, T::Type{<:Real}=Int, zeros=zeros) = count_kmers(seq, default_alphabet_size(eltype(seq)), K, T, zeros)
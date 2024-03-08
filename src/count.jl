"""
    count_kmers!(kmer_array, sequence; reset=true)
"""
function count_kmers! end

"""
    count_kmers!(kmer_array, sequence; reset=true)
"""
function count_kmers!(kmer_array::KmerArray{N, K}, sequence; reset::Bool = true) where {N, K}
    reset && zeros!(kmer_array)
    mask = N^K
    kmer = zero(UInt)
    for (i, m) in enumerate(sequence)
        kmer = kmer * N % mask + axis_index(kmer_array, m)
        kmer_array[kmer] += K <= i
    end
    return kmer_array
end

"""
    count_kmers(sequence, K, T=Int, zeros=zeros; N=default_alphabet_size(eltype(sequence)))
"""
function count_kmers end

function count_kmers(
    sequence, K::Integer, T::Type{<:Real}=Int, zeros=zeros;
    N::Integer = default_alphabet_size(eltype(sequence)),
)
    return count_kmers!(KmerArray(N, K, T, zeros), sequence; reset=false)
end
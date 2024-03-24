"""
    count_kmers!(kmer_array, sequence; reset=true)

Requires method `axis_index(::KmerArray{N}, ::eltype(sequence)) where N` to be defined
"""
function count_kmers!(kmer_array::KmerArray{N, K}, sequence; reset::Bool = true) where {N, K}
    reset && zeros!(kmer_array)
    kmer_array_values = kmer_array.values
    mask = N^K
    kmer = zero(UInt)
    for (i, m) in enumerate(sequence)
        kmer = (kmer * N + axis_index(kmer_array, m)) % mask
        kmer_array_values[kmer + 1] += K <= i
    end
    return kmer_array
end

"""
    count_kmers(sequence, K, T=Int, zeros=zeros; N=default_alphabet_size(eltype(sequence)))
"""
function count_kmers(sequence, K::Integer, T::Type{<:Real}=Int, zeros=zeros; N::Integer = default_alphabet_size(eltype(sequence)))
    return count_kmers!(KmerArray(N, K, T, zeros), sequence; reset=false)
end
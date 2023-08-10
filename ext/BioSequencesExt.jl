module BioSequencesExt
 
using VectorizedKmers, BioSequences
import VectorizedKmers: KmerCount, count_kmers!, KmerCountColumns, KmerCountRows

function count_kmers!(kmer_count::KmerCount{4, k, T}, seq::LongDNA{4}; reset::Bool=true) where {k, T}
    reset && fill!(kmer_count.data, zero(T))
    len = length(seq)
    mask = UInt(1) << 2k - 1
    kmer = UInt(0)
    i = 0
    @inbounds for data_int in seq.data
        for j in 0:4:63 # could maybe do some SIMD shit on middle k-mers
            i += 1
            i > len && break # only necessary for the last element in seq.data
            kmer = kmer << 2 & mask + trailing_zeros(data_int >> j) & 0b11
            kmer_count[kmer + 1] += k <= i
        end
    end
    kmer_count
end

end
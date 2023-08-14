module BioSequencesExt
 
using VectorizedKmers, BioSequences

function VectorizedKmers.count_kmers!(
    kmer_count::KmerCount{4, K, T},
    seq::LongDNA{4};
    reset::Bool = true,
) where {K, T}
    counts = kmer_count.counts
    reset && reset!(kmer_count)
    len = length(seq)
    mask = UInt(1) << 2K - 1
    kmer = UInt(0)
    i = 0
    @inbounds for data_int in seq.data
        for j in 0:4:63 # could maybe do some SIMD shit on middle K-mers
            i += 1
            i > len && break # only necessary for the last element in seq.data
            kmer = kmer << 2 & mask + trailing_zeros(data_int >> j) & 0b11
            counts[kmer + 1] += K <= i
        end
    end
    kmer_count
end

function VectorizedKmers.count_kmers(
    ::Type{KmerCount{4, K, T}},
    seq::LongDNA{4};
    zeros_func::Function = zeros,
) where {K, T}
    kmer_count = KmerCount{4, K, T}(zeros_func)
    count_kmers!(kmer_count, seq, reset=false)
    kmer_count
end

end
module BioSequencesExt
 
using VectorizedKmers, BioSequences

function VectorizedKmers.count_kmers!(
    kmer_count::KmerCountVector{4, k},
    seq::LongDNA{2};
    reset::Bool = true,
) where k
    counts = kmer_count.counts
    reset && reset!(kmer_count)
    len = length(seq)
    mask = UInt(1) << 2k - 1
    kmer = UInt(0)
    i = 0
    @inbounds for data_int in seq.data
        for j in 0:2:63
            i += 1
            i > len && break
            kmer = ((kmer << 2) & mask) | ((data_int >> j) & 0b11)
            counts[kmer + 1] += k <= i
        end
    end
    kmer_count
end

function VectorizedKmers.count_kmers!(
    kmer_count::KmerCountVector{4, k},
    seq::LongDNA{4};
    reset::Bool = true,
) where k
    counts = kmer_count.counts
    reset && reset!(kmer_count)
    len = length(seq)
    mask = UInt(1) << 2k - 1
    kmer = UInt(0)
    i = 0
    @inbounds for data_int in seq.data
        for j in 0:4:63 # could maybe do some SIMD shit on middle k-mers
            i += 1
            i > len && break # only necessary for the last element in seq.data
            kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
            counts[kmer + 1] += k <= i
        end
    end
    kmer_count
end

function VectorizedKmers.count_kmers(
    ::Type{KmerCountVector{4, k, T}},
    seq::LongDNA;
    zeros_func::Function = zeros,
) where {k, T}
    kmer_count = KmerCountVector{4, k, T}(zeros_func)
    count_kmers!(kmer_count, seq, reset=false)
    kmer_count
end

function VectorizedKmers.count_kmers(seq::LongDNA, K::Integer, T::DataType=UInt32)
    count_kmers(KmerCountVector{4, K, T}, seq)
end

end
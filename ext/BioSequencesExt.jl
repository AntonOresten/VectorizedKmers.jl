module BioSequencesExt
 
using VectorizedKmers, BioSequences

const LongNucOrView{N} = Union{LongSequence{<:NucleicAcidAlphabet{N}}, LongSubSeq{<:NucleicAcidAlphabet{N}}}

VectorizedKmers.alphabet_size(::Type{<:LongNucOrView}) = 4

@inline function VectorizedKmers.count_kmers!(
    kmer_count_vector::KmerCountVector{4, k},
    sequence::LongNucOrView{2};
    reset::Bool = true,
) where k
    reset && VectorizedKmers.zeros!(kmer_count_vector)
    counts = kmer_count_vector.counts
    len = length(sequence)
    mask = UInt(1) << 2k - 1
    kmer = UInt(0)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, len)
    first_count_index = k + start - 1
    i = 0
    @inbounds for data_int in sequence.data
        if i + 32 < start
            i += 32
            continue
        end
        for j in 0:2:63
            i += 1
            i < start && continue
            i > stop && break
            kmer = ((kmer << 2) & mask) | ((data_int >> j) & 0b11)
            counts[kmer + 1] += first_count_index <= i
        end
    end
    kmer_count_vector
end

@inline function VectorizedKmers.count_kmers!(
    kmer_count_vector::KmerCountVector{4, k},
    sequence::LongNucOrView{4};
    reset::Bool = true,
) where k
    reset && VectorizedKmers.zeros!(kmer_count_vector)
    counts = kmer_count_vector.counts
    len = length(sequence)
    mask = UInt(1) << 2k - 1
    kmer = UInt(0)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, len)
    first_count_index = k + start - 1
    i = 0
    @inbounds for data_int in sequence.data
        if i + 16 < start
            i += 16
            continue
        end
        for j in 0:4:63 # could maybe do some SIMD shit on middle k-mers
            i += 1
            i < start && continue
            i > stop && break
            kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
            counts[kmer + 1] += first_count_index <= i
        end
    end
    kmer_count_vector
end

end
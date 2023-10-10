module BioSequencesExt
 
using VectorizedKmers, BioSequences

import BioSequences: SeqOrView

VectorizedKmers.alphabet_size(::Type{<:SeqOrView{<:NucleicAcidAlphabet}}) = 4
VectorizedKmers.alphabet_size(::Type{<:SeqOrView{<:AminoAcidAlphabet}}) = 20

@inline function VectorizedKmers.count_kmers!(
    kmer_vector::KmerVector{4, k, T, A},
    sequence::SeqOrView{<:NucleicAcidAlphabet{2}};
    reset::Bool = true,
) where {k, T, A <: AbstractVector{T}}
    reset && VectorizedKmers.zeros!(kmer_vector)
    length(sequence) < k && return kmer_vector
    mask = UInt(4^k - 1)
    kmer = UInt(0)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, length(sequence))
    data_start, data_stop = cld(start, 32), cld(stop, 32)
    first_count_index = k + start - 1
    i = 32 * (data_start - 1)
    @inbounds for data_int in @view sequence.data[data_start:data_stop]
        for j in 0:2:63
            i += 1
            i > stop && break
            kmer = ((kmer << 2) & mask) | ((data_int >> j) & 0b11)
            kmer_vector.values[kmer + 1] += first_count_index <= i
        end
    end
    return kmer_vector
end

@inline function VectorizedKmers.count_kmers!(
    kmer_vector::KmerVector{4, k, T, A},
    sequence::SeqOrView{<:NucleicAcidAlphabet{4}};
    reset::Bool = true,
) where {k, T, A <: AbstractVector{T}}
    reset && VectorizedKmers.zeros!(kmer_vector)
    length(sequence) < k && return kmer_vector
    mask = UInt(4^k - 1)
    kmer = UInt(0)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, length(sequence))
    data_start, data_stop = cld(start, 16), cld(stop, 16)
    first_count_index = k + start - 1
    i = 16 * (data_start - 1)
    @inbounds for data_int in @view sequence.data[data_start:data_stop]
        for j in 0:4:63
            i += 1
            i > stop && break
            kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
            kmer_vector.values[kmer + 1] += first_count_index <= i
        end
    end
    return kmer_vector
end

function VectorizedKmers.count_kmers!(
    kmer_vector::KmerVector{20, k, T, A},
    sequence::SeqOrView{<:AminoAcidAlphabet};
    reset::Bool = true,
) where {k, T, A}
    reset && VectorizedKmers.zeros!(kmer_vector)
    length(sequence) < k && return kmer_vector
    mask = UInt(20^k)
    kmer = UInt(0)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, length(sequence))
    data_start, data_stop = cld(start, 8), cld(stop, 8)
    first_count_index = k + start - 1
    i = 8 * (data_start - 1)
    @inbounds for data_int in @view sequence.data[data_start:data_stop]
        for j in 0:8:63
            i += 1
            i > stop && break
            kmer = (kmer * 20 + ((data_int >> j) & 0xff) % 20) % mask
            kmer_vector.values[kmer + 1] += first_count_index <= i
        end
    end
    return kmer_vector
end

@inline function VectorizedKmers.count_kmers!(
    kmer_vector::KmerVector{4, k, Bool, BitVector},
    sequence::SeqOrView{<:NucleicAcidAlphabet{4}};
    reset::Bool = true,
) where k
    chunks = kmer_vector.values.chunks
    reset && fill!(chunks, zero(UInt))
    length(sequence) < k && return kmer_vector
    mask = UInt(4^k - 1)
    kmer = UInt(0)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, length(sequence))
    data_start, data_stop = cld(start, 16), cld(stop, 16)
    first_count_index = k + start - 1
    i = 16 * (data_start - 1)
    @inbounds for data_int in @view sequence.data[data_start:data_stop]
        for j in 0:4:63
            i += 1
            i > stop && break
            kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
            chunk_index = kmer >> 6 + one(UInt)
            index_in_chunk = kmer & 0x3f
            first_count_index <= i && (chunks[chunk_index] ⊻= one(UInt) << index_in_chunk)
        end
    end
    return kmer_vector
end

end
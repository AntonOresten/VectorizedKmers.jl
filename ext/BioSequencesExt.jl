module BioSequencesExt
 
using VectorizedKmers, BioSequences

import BioSequences: SeqOrView

VectorizedKmers.axis_index(::KmerArray{4}, m::NucleicAcid) = trailing_zeros(reinterpret(Int8, m))
VectorizedKmers.axis_index(::KmerArray{20}, m::AminoAcid) = (reinterpret(Int8, m) - 1) % 20 + 1

VectorizedKmers.default_alphabet_size(::Type{<:NucleicAcid}) = 4
VectorizedKmers.default_alphabet_size(::Type{AminoAcid}) = 20

function VectorizedKmers.count_kmers!(
    kmer_array::KmerArray{4, K, T, A},
    sequence::SeqOrView{<:NucleicAcidAlphabet{2}};
    reset::Bool = true,
) where {K, T, A}
    reset && VectorizedKmers.zeros!(kmer_array)
    mask = one(UInt) << 2K - 1
    kmer = zero(UInt)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, length(sequence))
    data_start, data_stop = cld(start, 32), cld(stop, 32)
    first_count_index = K + start - 1
    i = 32 * (data_start - 1)
    @inbounds for data_int in @view sequence.data[data_start:data_stop]
        for j in 0:2:63
            i += 1
            i > stop && break
            kmer = ((kmer << 2) & mask) | ((data_int >> j) & 0b11)
            kmer_array[kmer] += first_count_index <= i
        end
    end
    return kmer_array
end

function VectorizedKmers.count_kmers!(
    kmer_array::KmerArray{4, K, T, A},
    sequence::SeqOrView{<:NucleicAcidAlphabet{4}};
    reset::Bool = true,
) where {K, T, A}
    reset && VectorizedKmers.zeros!(kmer_array)
    mask = one(UInt) << 2K - 1
    kmer = zero(UInt)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, length(sequence))
    data_start, data_stop = cld(start, 16), cld(stop, 16)
    first_count_index = K + start - 1
    i = 16 * (data_start - 1)
    @inbounds for data_int in @view sequence.data[data_start:data_stop]
        for j in 0:4:63
            i += 1
            i > stop && break
            kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
            kmer_array.offset_values.parent[kmer + 1] += first_count_index <= i
        end
    end
    return kmer_array
end

function VectorizedKmers.count_kmers!(
    kmer_array::KmerArray{20, K, T, A},
    sequence::SeqOrView{AminoAcidAlphabet};
    reset::Bool = true,
) where {K, T, A}
    reset && VectorizedKmers.zeros!(kmer_array)
    mask = UInt(20^K)
    kmer = UInt(0)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, length(sequence))
    data_start, data_stop = cld(start, 8), cld(stop, 8)
    first_count_index = K + start - 1
    i = 8 * (data_start - 1)
    @inbounds for data_int in @view sequence.data[data_start:data_stop]
        for j in 0:8:63
            i += 1
            i > stop && break
            kmer = (kmer * 20 + ((data_int >> j) & 0xff) % 20) % mask
            kmer_array[kmer] += first_count_index <= i
        end
    end
    return kmer_array
end

function VectorizedKmers.count_kmers!(
    kmer_array::KmerArray{4, K, Bool, BitArray{K}},
    sequence::SeqOrView{<:NucleicAcidAlphabet{4}};
    reset::Bool = true,
) where K
    chunks = kmer_array.offset_values.parent.chunks
    reset && fill!(chunks, zero(UInt))
    mask = one(UInt) << 2K - 1
    kmer = zero(UInt)
    start, stop = sequence isa LongSubSeq ? (sequence.part.start, sequence.part.stop) : (1, length(sequence))
    data_start, data_stop = cld(start, 16), cld(stop, 16)
    first_count_index = K + start - 1
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
    return kmer_array
end

end
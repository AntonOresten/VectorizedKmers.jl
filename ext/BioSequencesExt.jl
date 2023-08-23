module BioSequencesExt
 
using VectorizedKmers, BioSequences

function VectorizedKmers.count_kmers!(
    kmer_count_vector::KmerCountVector{4, k},
    seq::LongDNA{2};
    reset::Bool = true,
) where k
    reset && zeros!(kmer_count_vector)
    counts = kmer_count_vector.counts
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
    kmer_count_vector
end

function VectorizedKmers.count_kmers!(
    kmer_count_vector::KmerCountVector{4, k},
    seq::LongDNA{4};
    reset::Bool = true,
) where k
    reset && zeros!(kmer_count_vector)
    counts = kmer_count_vector.counts
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
    kmer_count_vector
end

function VectorizedKmers.count_kmers(
    seq::LongDNA,
    k::Integer,
    T::Type{<:Real} = Int,
    zeros_func::Function = zeros,
)
    kmer_count_vector = KmerCountVector{4, k}(T, zeros_func)
    count_kmers!(kmer_count_vector, seq, reset=false)
    kmer_count_vector
end


function VectorizedKmers.count_kmers!(
    kmer_count_columns::KmerCountColumns{4, k},
    sequences::Vector{LongDNA{N}};
    column_offset::Integer = 0,
    reset::Bool = true
) where {k, N}
    kcv_gen = Iterators.drop(eachvec(kmer_count_columns), column_offset)
    for (kcv, seq) in zip(kcv_gen, sequences)
        count_kmers!(kcv, seq, reset=reset)
    end
    kmer_count_columns
end

function VectorizedKmers.count_kmers(
    sequences::Vector{LongDNA{4}},
    k::Integer,
    T::Type{<:Real} = Int,
    zeros_func::Function = zeros,
)
    kmer_count_columns = KmerCountColumns{4, k}(length(sequences), T, zeros_func)
    count_kmers!(kmer_count_columns, sequences, reset=false)
    kmer_count_columns
end

end

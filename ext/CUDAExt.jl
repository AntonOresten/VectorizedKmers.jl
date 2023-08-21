module CUDAExt

using VectorizedKmers, CUDA
using BioSequences

# transpose doesn't copy?
# store kmer count vectors in columns, and transpose when taking matmul
# store sequences in columns 

"""
    count_kmers!(kcc::KmerCountColumns{4, k, T, CuMatrix{T}}, sequences::CuMatrix{UInt8})

Turbocharged k-mer counting powered by CUDA. DNA sequences (must have same length) are stored in the columns of `sequences`.
Values in `sequences` must be between 0 and 3.

Chars 'A', 'C', 'G', and 'T' can be converted to 0, 1, 3, and 2 respectively using the function:
`char -> UInt8(char) >> 1 & 0x03`, or `byte -> byte >> 1 & 0x03`,
which is easy to broadcast to an array of bytes.

Only defined for kmer_count_columns because so that each k-mer count is stored contiguously in memory.

This is a good bare-bones method for when performance is critical, and I/O speeds are a bottleneck.
"""
function VectorizedKmers.count_kmers!(
    kmer_count_columns::KmerCountColumns{4, k, T, M},
    sequences::CuMatrix{UInt8},
    seq_lengths::Union{Missing, CuVector{Int}} = missing;
    column_offset::Int = 0,
    reset::Bool = true,
) where {k, T, M <: CuMatrix{T}}
    counts = kmer_count_columns.counts
    reset && CUDA.fill!(counts, 0)

    max_seq_len, num_seqs = size(sequences)
    seq_lengths = ismissing(seq_lengths) ? CUDA.fill(max_seq_len, num_seqs) : seq_lengths

    @assert length(kmer_count_columns) >= num_seqs "The k-mer counts of $num_seqs sequences would not fit in $(length(kmer_count_columns)) columns."
    @assert length(kmer_count_columns) - column_offset >= num_seqs "The column offset is too high. The sequences don't fit."
    @assert length(seq_lengths) == num_seqs "The length of the seq_lengths vector ($(length(seq_lengths))) does not match the number of sequences ($num_seqs)."
    @assert maximum(seq_lengths) <= max_seq_len "The maximum sequence length provided ($(maximum(seq_lengths))) is higher than the size of the sequences matrix ($max_seq_len)."
    max_seq_len > k-1 || @warn "All sequences are shorter than k. No k-mers will be counted."

    mask = unsigned(4^k - 1)

    function count_kmers_column!(counts, sequences, seq_lengths, num_seqs, k, mask)
        seq_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x

        @inbounds if seq_idx <= num_seqs
            count_vector = view(counts, :, seq_idx + column_offset)
            sequence = view(sequences, :, seq_idx)

            kmer = zero(UInt64)
            for i in 1:k-1
                base = sequence[i]
                kmer = (kmer << 2) | base
            end

            for i in k:seq_lengths[seq_idx]
                base = sequence[i]
                kmer = ((kmer << 2) & mask) | base
                count_vector[kmer + 1] += one(T)
            end
        end
    end

    threads = 256
    blocks = ceil(Int, num_seqs / threads)

    @cuda threads=threads blocks=blocks count_kmers_column!(counts, sequences, seq_lengths, num_seqs, k, mask)

    kmer_count_columns
end

# maybe have a method that takes vector of sequences (as CuVectors) and converts to CuMatrix with sequences padded to max seq len and calls the method above

@inline data_length(seq::LongDNA{4}) = length(seq.data) 

function VectorizedKmers.count_kmers!(
    kmer_count_columns::KmerCountColumns{4, k, T, M},
    sequences::Vector{LongDNA{4}};
    column_offset::Int = 0,
    reset::Bool = true,
) where {k, T, M <: CuMatrix{T}}
    counts = kmer_count_columns.counts
    reset && CUDA.fill!(counts, 0)

    num_seqs = length(sequences)

    seq_lengths = CuVector(length.(sequences))
    data_lengths = CuVector(data_length.(sequences))
    #max_seq_len = maximum(seq_lengths)
    max_data_length = maximum(data_lengths)

    data_matrix_h = Matrix{UInt64}(undef, (max_data_length, num_seqs))
    for (i, seq) in enumerate(sequences)
        data_matrix_h[1:data_length(seq), i] = seq.data
    end

    data_matrix = CuMatrix(data_matrix_h)

    function count_kmers_column!(
        counts,
        data_matrix,
        seq_lengths,
        num_seqs,
        k,
        mask,
        data_lengths,
        column_offset,
    )
        seq_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x

        @inbounds if seq_idx <= num_seqs
            seq_length = seq_lengths[seq_idx]
            data_length = data_lengths[seq_idx]

            count_vector = view(counts, :, seq_idx + column_offset)
            data_vector = view(data_matrix, 1:data_length-1, seq_idx)

            kmer = UInt(0)
            i = 0
            for data_int in data_vector
                for j in 0:4:63
                    i += 1
                    kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
                    count_vector[kmer + 1] += k <= i
                end
            end

            data_int = data_vector[data_length]
            for j in 0:4:(4 * ((seq_length - 1) % 16 + 1) - 1)
                i += 1
                kmer = ((kmer << 2) & mask) | (trailing_zeros(data_int >> j) & 0b11)
                count_vector[kmer + 1] += k <= i
            end
        end

        nothing
    end

    mask = unsigned(4^k - 1)
    threads = 256
    blocks = ceil(Int, num_seqs / threads)

    @cuda threads=threads blocks=blocks count_kmers_column!(
        counts, data_matrix, seq_lengths, num_seqs, k, mask, data_lengths, column_offset)

    kmer_count_columns
end

end
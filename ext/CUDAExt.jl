module CUDAExt

using VectorizedKmers, CUDA

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
"""
function VectorizedKmers.count_kmers!(
    kmer_count_columns::KmerCountColumns{4, k, T, M},
    sequences::CuMatrix{UInt8},
    sequence_lengths::Union{Missing, CuVector{Int}} = missing;
    column_offset::Int = 0,
    reset::Bool = true,
) where {k, T, M <: CuMatrix{T}}
    count_matrix = kmer_count_columns.counts
    max_seq_len, num_sequences = size(sequences)
    sequence_lengths = ismissing(sequence_lengths) ? CUDA.fill(max_seq_len, num_sequences) : sequence_lengths
    reset && CUDA.fill!(count_matrix, 0)

    @assert length(kmer_count_columns) >= num_sequences "The k-mer counts of $num_sequences sequences would not fit in $(length(kmer_count_columns)) columns."
    @assert length(kmer_count_columns) - column_offset >= num_sequences "The column offset is too high. Sequences can't fit."
    @assert length(sequence_lengths) == num_sequences "The length of the sequence_lengths vector ($(length(sequence_lengths))) does not match the number of sequences ($num_sequences)."
    @assert maximum(sequence_lengths) <= max_seq_len "The maximum sequence length provided ($(maximum(sequence_lengths))) is higher than the size of the sequences matrix ($max_seq_len)."
    max_seq_len > k-1 && @warn "All sequences are shorter than k. No k-mers will be counted."

    function count_kmers_column!(count_matrix, sequences, sequence_lengths, k, mask)
        seq_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x

        @inbounds if seq_idx <= num_sequences
            count_vector = view(count_matrix, :, seq_idx + column_offset)
            sequence = view(sequences, :, seq_idx)

            kmer = zero(UInt64)
            for i in 1:k-1
                base = sequence[i]
                kmer = (kmer << 2) | base
            end

            for i in k:sequence_lengths[seq_idx]
                base = sequence[i]
                kmer = ((kmer << 2) & mask) | base
                count_vector[kmer + 1] += one(T)
            end
        end
    end

    mask = unsigned(4^k - 1)
    threads = 256
    blocks = ceil(Int, num_sequences / threads)

    @cuda threads=threads blocks=blocks count_kmers_column!(count_matrix, sequences, sequence_lengths, k, mask)

    kmer_count_columns
end

end
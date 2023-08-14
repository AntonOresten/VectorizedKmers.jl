module CUDAExt

export
    count_kmers_cuda!

using VectorizedKmers, CUDA
import VectorizedKmers: KmerCount, count_kmers!, KmerCountColumns, KmerCountRows

# transpose doesn't copy?
# store kmer count vectors in columns, and transpose when taking matmul
# store sequences in columns 

"""
    count_kmers!(kcc::KmerCountColumns{4, K, T, CuMatrix{T}}, sequences::CuMatrix{UInt8})

Turbocharged K-mer counting powered by CUDA. DNA sequences (must have same length) are stored in the columns of `sequences`.
Values in `sequences` must be between 0 and 3.

Chars 'A', 'C', 'G', and 'T' can be converted to 0, 1, 3, and 2 respectively using the function:
`char -> UInt8(char) >> 1 & 0x03`, or `byte -> byte >> 1 & 0x03`,
which is easy to broadcast to an array of bytes.
"""
function VectorizedKmers.count_kmers!(
    kmer_count_columns::KmerCountColumns{4, K, T, M},
    sequences::CuMatrix{UInt8};
    reset::Bool = true,
) where {K, T, M <: CuMatrix{T}}
    count_matrix = kmer_count_columns.counts
    seq_len, num_sequences = size(sequences)
    reset && CUDA.fill!(count_matrix, 0)

    function count_kmers_column!(count_matrix, sequences, K, seq_len, num_sequences, mask)
        seq_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x

        @inbounds if seq_idx <= num_sequences
            count_vector = view(count_matrix, :, seq_idx)
            sequence = view(sequences, :, seq_idx)

            kmer = zero(UInt64)
            for i in 1:K-1
                base = sequence[i]
                kmer = (kmer << 2) + base
            end

            for i in K:seq_len
                base = sequence[i]
                kmer = ((kmer << 2) & mask) + base
                count_vector[kmer + 1] += one(T)
            end
        end

        return nothing
    end

    mask = unsigned(4^K - 1)
    threads = 256
    blocks = ceil(Int, num_sequences / threads)

    @cuda threads=threads blocks=blocks count_kmers_column!(count_matrix, sequences, K, seq_len, num_sequences, mask)

    kmer_count_columns
end

end
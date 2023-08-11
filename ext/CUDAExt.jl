module CUDAExt

export
    count_kmers_cuda!

using VectorizedKmers, CUDA
import VectorizedKmers: KmerCount, count_kmers!, KmerCountColumns, KmerCountRows

# transpose doesn't copy?
# store kmer count vectors in columns, and transpose when taking matmul
# store sequences in columns 

"""
    count_kmers!(kcc::KmerCountColumns{4, k, T, CuMatrix{T}}, sequences::CuMatrix{UInt8})

Turbocharged k-mer counting using CUDA. DNA sequences (all with same length) are stored in the columns of `sequences`.
Values in `sequences` must be between 0 and 3.

Chars 'A', 'C', 'G', and 'T' can be converted to 0, 1, 3, and 2 respectively using the function:
`char -> UInt8(char) >> 1 & 0x03`
"""
function count_kmers!(kmer_count_columns::KmerCountColumns{4, k, T, M}, sequences::CuMatrix{UInt8}) where {k, T, M <: CuMatrix{T}}
    data = kmer_count_columns.data
    seq_len, num_sequences = size(sequences)
    CUDA.fill!(data, zero(T))

    function kernel(data, sequences, k, seq_len, num_sequences, mask)
        seq_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if seq_idx <= num_sequences
            kmer = unsigned(0)
            for i in 1:k-1
                base = sequences[i, seq_idx]
                kmer = kmer << 2 + base
            end
            for i in k:seq_len
                base = sequences[i, seq_idx]
                kmer = kmer << 2 & mask + base
                CUDA.@atomic data[kmer + 1, seq_idx] += one(T)
            end
        end
        return
    end

    mask = unsigned(4^k - 1)
    threads = 256
    blocks = ceil(Int, num_sequences / threads)

    @cuda threads=threads blocks=blocks kernel(data, sequences, k, seq_len, num_sequences, mask)

    kmer_count_columns
end

end
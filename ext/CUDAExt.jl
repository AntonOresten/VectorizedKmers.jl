module CUDAExt

export
    count_kmers_cuda!

using VectorizedKmers, CUDA
import VectorizedKmers: KmerCount, count_kmers!, KmerCountColumns, KmerCountRows

function count_kmers_cuda!(kmer_count::KmerCount{4, k, T}, seq::CuVector{UInt8}; reset::Bool=true) where {k, T}
    data = kmer_count.data
    seq_len = length(seq)
    reset && CUDA.fill!(data, zero(T))

    mask = unsigned(4^k - 1)
    kmer = unsigned(0)
    for i in 1:k-1
        base = seq[i]
        kmer = (kmer << 2) + base
    end
    for i in k:seq_len
        base = seq[i]
        kmer = ((kmer << 2) & mask) + base
        CUDA.@atomic data[kmer + 1] += one(T)
    end
end

function count_kmers_cuda!(kmer_count_vector::AbstractKmerCountVector{4, k, T, CuMatrix{T}}, sequences::CuMatrix{UInt8}; reset::Bool=true) where {k, T}
    num_sequences = kmer_count_vector isa KmerCountColumns ? size(sequences, 2) : size(sequences, 1)
    @assert length(kmer_count_vector) >= num_sequences
    reset && fill!(kmer_count_vector.data, zero(T))

    function kernel(kmer_count_vector::AbstractKmerCountVector{4, k, T, CuMatrix{T}}, sequences::CuMatrix{UInt8})
        seq_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
        if seq_idx <= num_sequences
            count_kmers_cuda!(kmer_count_vector[seq_idx], view(sequences, :, seq_idx), reset=reset)
        end
        return
    end

    threads = 256
    blocks = ceil(Int, num_sequences / threads)
    @cuda threads=threads blocks=blocks kernel(kmer_count_vector, sequences)

    kmer_count_vector
end

end
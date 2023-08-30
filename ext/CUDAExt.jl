module CUDAExt

using VectorizedKmers, CUDA

# transpose doesn't copy?
# store kmer count vectors in columns, and transpose for matmul

# maybe not needed?
#@inline zeros!(ka::AbstractKmerArray{N, S, k, T, A <: CuArray}) where {N, S, k, T, A} = CUDA.fill!(ka.values, zero(eltype(ka)))

"""
    count_kmers!(kc::KmerColumns{4, k, T, CuMatrix{T}}, sequences::CuMatrix{UInt8})

Turbocharged k-mer counting powered by CUDA. DNA sequences (must have same length) are stored in the columns of `sequences`.
Values in `sequences` must be between 0 and 3.

Chars 'A', 'C', 'G', and 'T' can be converted to 0, 1, 3, and 2 respectively using the function:
`char -> UInt8(char) >> 1 & 0x03`, or `byte -> byte >> 1 & 0x03`,
which is easy to broadcast to an array of bytes.

Only defined for kmer_columns because so that each k-mer count is stored contiguously in memory.

This is a good bare-bones method for when performance is critical, and I/O speeds are a bottleneck.
"""
function VectorizedKmers.count_kmers!(
    kmer_columns::KmerColumns{4, k, T, M},
    sequences::CuMatrix{UInt8},
    seq_lengths::Union{Missing, CuVector{Int}} = missing;
    offset::Int = 0,
    reset::Bool = true,
) where {k, T, M <: CuMatrix{T}}
    values = kmer_columns.values
    reset && CUDA.fill!(values, 0)

    max_seq_len, num_seqs = size(sequences)
    seq_lengths = ismissing(seq_lengths) ? CUDA.fill(max_seq_len, num_seqs) : seq_lengths

    @assert length(kmer_columns) >= num_seqs "The k-mer counts of $num_seqs sequences would not fit in $(length(kmer_columns)) columns."
    @assert length(kmer_columns) - offset >= num_seqs "The column offset is too high. The sequences don't fit."
    @assert length(seq_lengths) == num_seqs "The length of the seq_lengths vector ($(length(seq_lengths))) does not match the number of sequences ($num_seqs)."
    @assert maximum(seq_lengths) <= max_seq_len "The maximum sequence length provided ($(maximum(seq_lengths))) is higher than the size of the sequences matrix ($max_seq_len)."
    max_seq_len > k-1 || @warn "All sequences are shorter than k. No k-mers will be counted."

    mask = unsigned(4^k - 1)

    function count_kmers_column!(values, sequences, seq_lengths, num_seqs, k, mask)
        seq_idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x

        @inbounds if seq_idx <= num_seqs
            count_vector = view(values, :, seq_idx + offset)
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

    @cuda threads=threads blocks=blocks count_kmers_column!(values, sequences, seq_lengths, num_seqs, k, mask)

    kmer_columns
end

# maybe have a method that takes vector of sequences (as CuVectors) and converts to CuMatrix with sequences padded to max seq len and calls the method above

end
"""
    count_kmers!(kmer_count, sequence)

Mutates `kmer_count`.

If `reset` is `true`, the array will be zero-ed before counting.
"""
function count_kmers! end

#=  The k-mers in `kmers` must be represented as integers between 0 and length(kmer_count) - 1.
    This is not a very efficient method, since it takes an entire vector. It is mainly used for testing.
    Ideally the k-mers would be procedurally calculated in constant memory. =#
function count_kmers!(
    kmer_count::KmerCountVector{S, k},
    kmers::Vector{<:Integer};
    reset::Bool = true,
) where {S, k}
    @assert maximum(kmers) < S^k
    reset && zeros!(kmer_count)
    for kmer in kmers
        kmer_count[kmer + 1] += 1
    end
    kmer_count
end


"""
    count_kmers(sequence, k; zeros_func=zeros, reset=true)

Create a new S^k sized vector using `zeros_func` and count the k-mers in `kmers`.
"""
function count_kmers end

function count_kmers(
    ::Type{KmerCountVector{S, k, T}}, kmers::Vector{<:Integer};
    zeros_func::Function = zeros,
) where {S, k, T}
    kmer_count = KmerCountVector{S, k, T}(zeros_func)
    count_kmers!(kmer_count, kmers, reset=false)
    return kmer_count
end


"""
    count_kmers_gpu(sequences, k)

Count the k-mers in `sequences` using the GPU.
"""
function count_kmers_gpu end
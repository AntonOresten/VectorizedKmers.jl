"""
    count_kmers!(kmer_count_vector, sequence)

Mutates `kmer_count_vector`.

If `reset` is `true`, the array will be zero-ed before counting.
"""
function count_kmers! end

#=  The k-mers in `kmers` must be represented as integers between 0 and length(kmer_count_vector) - 1.
    This is not a very efficient method, since it takes an entire vector. It is mainly used for testing.
    Ideally the k-mers would be procedurally calculated in constant memory. =#
function count_kmers!(
    kmer_count_vector::KmerCountVector{S, k},
    kmers::Vector{<:Integer};
    reset::Bool = true,
) where {S, k}
    @assert maximum(kmers) < S^k
    reset && zeros!(kmer_count_vector)
    counts = kmer_count_vector.counts
    for kmer in kmers
        counts[kmer + 1] += 1
    end
    kmer_count_vector
end


"""
    count_kmers(sequence, k; zeros_func=zeros, reset=true)

Create a new S^k sized vector using `zeros_func` and count the k-mers in `kmers`.
"""
function count_kmers end

function count_kmers(
    kmers::Vector{<:Integer},
    S::Integer,
    k::Integer,
    T::Type{<:Real} = Int,
    zeros_func::Function = zeros,
)
    kmer_count_vector = KmerCountVector{S, k}(T, zeros_func)
    count_kmers!(kmer_count_vector, kmers, reset=false)
    kmer_count_vector
end
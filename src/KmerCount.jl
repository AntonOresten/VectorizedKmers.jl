"""
    AbstractKmerCount{A, K, T <: Real} <: AbstractVector{T}

Abstract type for k-mer counts. `A` is the alphabet size, `k` is the k-mer size,
and `T` is the element type of any array of counts that a subtype might have.
"""
abstract type AbstractKmerCount{A, K, T <: Real} <: AbstractVector{T} end

@inline Base.size(kmer_count::AbstractKmerCount) = size(kmer_count.counts)
@inline Base.length(kmer_count::AbstractKmerCount) = length(kmer_count.counts)
@inline Base.getindex(kmer_count::AbstractKmerCount, i::Integer) = kmer_count.counts[i]
@inline Base.setindex!(kmer_count::AbstractKmerCount, v::Real, i::Integer) = kmer_count.counts[i] = v

"""
    KmerCount{A, K, T} <: AbstractKmerCount{A, K, T}

A k-mer count vector. `A` is the alphabet size, `k` is the k-mer size,
and `T` is the element type of `counts`.
"""
struct KmerCount{A, K, T} <: AbstractKmerCount{A, K, T}
    counts::AbstractVector{T}

    function KmerCount{A, K, T}(counts::AbstractVector{T}) where {A, K, T}
        new{A, K, T}(counts)
    end

    function KmerCount{A, K, T}() where {A, K, T}
        KmerCount{A, K, T}(zeros(T, A^K))
    end
end

"""
    change_k(kmer_count, new_K)

Change the k-mer size of `kmer_count` to `new_K`. The new K value must be less than or equal to the old one.
If a 4-mer count gets changed to a 3-mer count, the count of `ACGT` will be added to both `ACG` and `CGT`
"""
function change_k(
    kmer_count::KmerCount{A, K, T},
    new_K::Integer,
) where {A, K, T}
    right_shifts = K - new_K
    @assert right_shifts > 0 "Attempt at changing K-mer count from $K to $new_K failed. New K must be less than or equal to old K."
    mask = A^new_K
    new_kmer_count = KmerCount{A, new_K, T}()
    for (kmer_plus_one, n) in enumerate(kmer_count)
        iszero(n) && continue
        kmer = kmer_plus_one - 1
        new_kmer_count[kmer % mask + 1] += n
        for _ in 1:right_shifts
            kmer ÷= A
            new_kmer_count[kmer % mask + 1] += n
        end
    end
    new_kmer_count
end

"""
    count_kmers!(kmer_count, kmers; reset=true)

Mutate the `counts` vector in `kmer_count` by adding the counts of each kmer in `kmers`.
The K-mers in `kmers` must be represented as integers between 0 and length(kmer_count) - 1.

If `reset` is `true`, the `counts` vector will be zero-ed before counting.

This is not a very efficient method, since it takes an entire vector. It is mainly used for testing.
Ideally the K-mers would be procedurally calculated in constant memory.
"""
function count_kmers!(kmer_count::KmerCount{A, K, T}, kmers::Vector{<:Integer}; reset::Bool=true) where {A, K, T}
    reset && fill!(kmer_count, zero(T))
    for kmer in kmers
        kmer_count[kmer + 1] += 1
    end
    kmer_count
end
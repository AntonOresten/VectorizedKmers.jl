"""
    AbstractKmerCount{A, K, T <: Real, V <: AbstractVector{T}}

Abstract type for K-mer counts.
`A` is the alphabet size,
`K` is the K-mer size,
and `T` is the element type of the underlying `counts` field,
which in turn has type `V{T}`.
"""
abstract type AbstractKmerCount{A, K, T <: Real, V <: AbstractVector{T}} <: AbstractVector{T} end

@inline Base.size(kmer_count::AbstractKmerCount) = size(kmer_count.counts)
@inline Base.length(kmer_count::AbstractKmerCount) = length(kmer_count.counts)
@inline Base.getindex(kmer_count::AbstractKmerCount, i::Integer) = kmer_count.counts[i]
@inline Base.setindex!(kmer_count::AbstractKmerCount, v::Real, i::Integer) = kmer_count.counts[i] = v
@inline Base.eltype(::AbstractKmerCount{A, K, T}) where {A, K, T} = T

"""
    KmerCount{A, K, T, V} <: AbstractKmerCount{A, K, T, V}

A concrete type for K-mer counts with vector type `Base.Vector{T}`.
"""
struct KmerCount{A, K, T, V} <: AbstractKmerCount{A, K, T, V}
    counts::V

    function KmerCount{A, K, T}(counts::V) where {A, K, T, V <: AbstractVector}
        new{A, K, T, V}(counts)
    end

    function KmerCount{A, K, T}() where {A, K, T}
        KmerCount{A, K, T}(zeros(T, A^K))
    end
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
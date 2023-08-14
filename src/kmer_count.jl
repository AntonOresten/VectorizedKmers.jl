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

A concrete type for K-mer counts with vector type `V` and element type `T`.
"""
struct KmerCount{A, K, T, V} <: AbstractKmerCount{A, K, T, V}
    counts::V

    function KmerCount{A, K}(counts::V) where {A, K, T <: Real, V <: AbstractVector{T}}
        @assert length(counts) == A^K
        new{A, K, T, V}(counts)
    end

    function KmerCount{A, K, T}(zeros_func::Function = zeros) where {A, K, T}
        KmerCount{A, K}(zeros_func(T, A^K))
    end
end

KmerCount{A, K}(zf::Function = zeros) where {A, K} = KmerCount{A, K, Int}(zf)

"""
    count_kmers!(kmer_count, kmers; reset=true)

Mutate the `counts` vector in `kmer_count` by counting K-mers in `kmers`.
The K-mers in `kmers` must be represented as integers between 0 and length(kmer_count) - 1.

If `reset` is `true`, the `counts` vector will be zero-ed before counting.
"""
function count_kmers!(
    kmer_count::KmerCount{A, K, T},
    kmers::Vector{<:Integer};
    reset::Bool = true,
) where {A, K, T}
    reset && fill!(kmer_count.counts, zero(T))
    for kmer in kmers
        kmer_count[kmer + 1] += 1
    end
    kmer_count
end
#=  This is not a very efficient method, since it takes an entire vector. It is mainly used for testing.
    Ideally the K-mers would be procedurally calculated in constant memory. =#

"""
    count_kmers(KmerCount{A, K, T}, kmers; zeros_func=zeros, reset=true)

Create a new A^K sized vector using `zeros_func` and count the K-mers in `kmers`.
The K-mers in `kmers` must be represented as integers between `0` and `length(kmer_count) - 1`.

If `reset` is `true`, the `counts` vector will be zero-ed before counting.
"""
function count_kmers(
    ::KmerCount{A, K, T}, kmers::Vector{<:Integer};
    zeros_func::Function = zeros,
) where {A, K, T}
    kmer_count = KmerCount{A, K, T}(zeros_func)
    count_kmers!(kmer_count, kmers, reset=false)
    kmer_count
end
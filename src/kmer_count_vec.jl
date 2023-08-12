"""
    AbstractKmerCountVector{A, K, T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{KmerCount{A, K, T}}

A container for K-mer counts, where K-mer counts are stored together as rows or columns in a matrix.
`A` is the alphabet size, `K` is the K-mer size, `T` is the element type of the counts,
and `M` is the type of the matrix in which the K-mer counts are stored.
"""
abstract type AbstractKmerCountVector{A, K, T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{KmerCount{A, K, T}} end

struct KmerCountColumns{A, K, T, M} <: AbstractKmerCountVector{A, K, T, M}
    counts::M

    function KmerCountColumns{A, K, T}(counts::M) where {A, K, T, M <: AbstractMatrix{T}}
        @assert size(counts, 1) == A^K
        new{A, K, T, M}(counts)
    end
end

@inline Base.size(kmer_count::KmerCountColumns) = (size(kmer_count.counts, 2),)
@inline Base.length(kmer_count::KmerCountColumns) = size(kmer_count.counts, 2)

@inline function Base.getindex(kmer_count::KmerCountColumns{A, K, T}, i::Integer) where {A, K, T}
    KmerCount{A, K, T}(view(kmer_count.counts, :, i))
end

struct KmerCountRows{A, K, T, M} <: AbstractKmerCountVector{A, K, T, M}
    counts::M

    function KmerCountRows{A, K, T}(counts::M) where {A, K, T, M <: AbstractMatrix{T}}
        @assert size(counts, 2) == A^K
        new{A, K, T, M}(counts)
    end
end

@inline Base.size(kmer_count::KmerCountRows) = (size(kmer_count.counts, 1),)
@inline Base.length(kmer_count::KmerCountRows) = size(kmer_count.counts, 1)
@inline Base.getindex(kmer_count::KmerCountRows{A, K, T}, i::Integer) where {A, K, T} = KmerCount{A, K, T}(view(kmer_count.counts, i, :))

function KmerCountColumns(kcc::KmerCountRows{A, K, T}) where {A, K, T}
    KmerCountColumns{A, K, T}(transpose(kcc.counts))
end

function KmerCountRows(kcc::KmerCountColumns{A, K, T}) where {A, K, T}
    KmerCountRows{A, K, T}(transpose(kcc.counts))
end

Base.transpose(kcc::KmerCountColumns) = KmerCountRows(kcc)
Base.transpose(kcr::KmerCountRows) = KmerCountColumns(kcr)

Base.adjoint(kcc::KmerCountColumns) = transpose(kcc)
Base.adjoint(kcr::KmerCountRows) = transpose(kcr)
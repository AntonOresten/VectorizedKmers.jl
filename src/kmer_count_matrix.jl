"""
    AbstractKmerCountMatrix{S, k, T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{KmerCountVector{S, k, T, V} where {V <: AbstractVector{T}}}

A container for k-mer counts, where k-mer counts are stored together as rows or columns in a matrix.
`A` is the alphabet size, `k` is the k-mer size, `T` is the element type of the counts,
and `M` is the type of the matrix in which the k-mer counts are stored.

Acts both like a vector of `KmerCountVector`s and like a matrix of k-mer counts.
"""
abstract type AbstractKmerCountMatrix{S, k, T, M} <: AbstractKmerCounts{2, S, k, T, M} end

@inline Base.size(kcm::AbstractKmerCountMatrix) = size(kcm.counts)
@inline Base.eachcol(kcm::AbstractKmerCountMatrix) = eachcol(kcm.counts)
@inline Base.eachrow(kcm::AbstractKmerCountMatrix) = eachrow(kcm.counts)
@inline Base.eachindex(kcm::AbstractKmerCountMatrix) = Base.OneTo(length(kcm))
@inline eachvec(kcm::AbstractKmerCountMatrix) = (kcm[i] for i in eachindex(kcm))


"""
    KmerCountColumns{S, k, T, M} <: AbstractKmerCountMatrix{S, k, T, M}

A container for k-mer counts, where k-mer counts are stored together as columns in a matrix.
This is more efficient than storing k-mer counts as rows in a matrix, since the elements in a column are contiguous in memory.
"""
struct KmerCountColumns{S, k, T, M} <: AbstractKmerCountMatrix{S, k, T, M}
    counts::M

    function KmerCountColumns{S, k}(counts::M) where {S, k, T <: Real, M <: AbstractMatrix{T}}
        @assert size(counts, 1) == S^k
        new{S, k, T, M}(counts)
    end
end

@inline Base.length(kcc::KmerCountColumns) = size(kcc, 2)

# overrides getindex(::AbstractKmerCounts, ::Any) from kmer_count.jl
@inline function Base.getindex(kcc::KmerCountColumns{S, k}, i::Integer) where {S, k}
    KmerCountVector{S, k}(view(kcc.counts, :, i))
end


"""
    KmerCountRows{S, k, T, M} <: AbstractKmerCountMatrix{S, k, T, M}

A container for k-mer counts, where k-mer counts are stored together as rows in a matrix.
This is not as efficient as storing k-mer counts as columns in a matrix, since the elements in a row are not contiguous in memory.
"""
struct KmerCountRows{S, k, T, M} <: AbstractKmerCountMatrix{S, k, T, M}
    counts::M

    function KmerCountRows{S, k}(counts::M) where {S, k, T <: Real, M <: AbstractMatrix{T}}
        @assert size(counts, 2) == S^k
        new{S, k, T, M}(counts)
    end
end

@inline Base.length(kcr::KmerCountRows) = size(kcr, 1)

# overrides getindex(::AbstractKmerCounts, ::Any) from kmer_count.jl
@inline function Base.getindex(kmer_count::KmerCountRows{S, k}, i::Integer) where {S, k}
    KmerCountVector{S, k}(view(kmer_count.counts, i, :))
end


function KmerCountColumns(kcc::KmerCountRows{S, k}) where {S, k}
    KmerCountColumns{S, k}(transpose(kcc.counts))
end

function KmerCountRows(kcr::KmerCountColumns{S, k}) where {S, k}
    KmerCountRows{S, k}(transpose(kcr.counts))
end

@inline Base.transpose(kcr::KmerCountRows) = KmerCountColumns(kcr)
@inline Base.transpose(kcc::KmerCountColumns) = KmerCountRows(kcc)

@inline Base.adjoint(kcm::AbstractKmerCountMatrix) = transpose(kcm)
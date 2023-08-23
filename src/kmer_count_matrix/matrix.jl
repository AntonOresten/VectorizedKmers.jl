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
@inline Base.getindex(kcm::AbstractKmerCountMatrix, r) = [kcm[i] for i in r]
@inline eachvec(kcm::AbstractKmerCountMatrix) = (kcm[i] for i in 1:length(kcm))

include("columns.jl")
include("rows.jl")
include("conversion.jl")
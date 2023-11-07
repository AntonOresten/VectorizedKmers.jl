"""
    AbstractKmerArray{N, S, k, T, A <: AbstractArray{T, N}} <: AbstractArray{T, N}

Abstract type representing k-mer count in N-dimensional arrays. This type serves as the foundation for various concrete k-mer count structures.

# Type Parameters
- `N`: Dimensionality of the array
- `S`: Alphabet size
- `k`: k-mer size
- `T`: Element type
- `A`: Array type, subtype of `AbstractArray{T, N}`
"""
abstract type AbstractKmerArray{N, S, k, T, A <: AbstractArray{T, N}} <: AbstractArray{T, N} end

const AbstractKmerVector = AbstractKmerArray{1}
const AbstractKmerMatrix = AbstractKmerArray{2}

@inline Base.size(ka::AbstractKmerArray) = size(ka.values)
@inline Base.length(ka::AbstractKmerArray) = length(ka.values)
@inline Base.getindex(ka::AbstractKmerArray, i) = ka.values[i]
@inline Base.getindex(ka::AbstractKmerArray, i, j) = ka.values[i, j]
@inline Base.setindex!(ka::AbstractKmerArray, v, i) = ka.values[i] = v
@inline Base.setindex!(ka::AbstractKmerArray, v, i, j) = ka.values[i, j] = v
@inline zeros!(ka::AbstractKmerArray) = fill!(ka.values, zero(eltype(ka)))

# also see the ::KmerVectors{D} method in src/vectors.jl
Base.hash(ka::AbstractKmerArray{N, S}, h::UInt) where {N, S} = hash(ka.values, h ⊻ S)

Base.:(==)(ka1::AbstractKmerArray, ka2::AbstractKmerArray) = hash(ka1) == hash(ka2)
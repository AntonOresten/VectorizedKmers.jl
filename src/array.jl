"""
    AbstractKmerArray{N, S, k, T <: Real, A <: AbstractArray{T, N}} <: AbstractArray{T, N}

Abstract type representing k-mer count in N-dimensional arrays. This type serves as the foundation for various concrete k-mer count structures.

# Type Parameters
- `N`: Dimensionality of the array.
- `S`: Alphabet size.
- `k`: k-mer size.
- `T`: Element type, subtype of `Real`.
- `A`: Array type, subtype of `AbstractArray{T, N}`.

## Rationale for Type Parameter Ordering
The order of the type parameters is purposefully chosen based on the logical importance of comparisons:

1. `N`: Dimensionality is paramount. Types with differing dimensions cannot be compared directly. They would somehow
   need to be transformed to the same dimension.
2. `S`: Alphabet size is a crucial differentiator. It's nonsensical to compare types with distinct alphabet sizes,
   even if k-mer sizes and other attributes are the same.
3. `k`: k-mer size is the next important type parameter. Types with varying k-mer sizes shouldn't be directly compared 
   even if they share other attributes.
4. `T`: Element type might imply certain transformations or states (e.g., normalization when using floats). Types with 
   different element types can have distinct implications and hence shouldn't be equated.
5. `A`: Lastly, the specific array type (like `Vector` or `SparseVector`) is the most flexible comparison point. Types 
   might be comparable if they share all other parameters, even if their array types differ.

This order ensures that the most significant distinctions between types are made first, allowing for more meaningful
and nuanced comparisons.
"""
abstract type AbstractKmerArray{N, S, k, T <: Real, A <: AbstractArray{T, N}} <: AbstractArray{T, N} end

const AbstractKmerVector = AbstractKmerArray{1}
const AbstractKmerMatrix = AbstractKmerArray{2}

@inline Base.size(ka::AbstractKmerArray) = size(ka.values)
@inline Base.length(ka::AbstractKmerArray) = length(ka.values)
@inline Base.getindex(ka::AbstractKmerArray, i) = ka.values[i]
@inline Base.getindex(ka::AbstractKmerArray, i, j) = ka.values[i, j]
@inline Base.setindex!(ka::AbstractKmerArray, v::Real, i) = ka.values[i] = v
@inline Base.setindex!(ka::AbstractKmerArray, v::Real, i, j) = ka.values[i, j] = v
@inline get_S(::AbstractKmerArray{N, S}) where {N, S} = S
@inline get_k(::AbstractKmerArray{N, S, k}) where {N, S, k} = k
@inline zeros!(ka::AbstractKmerArray) = fill!(ka.values, zero(eltype(ka)))

# also see the ::KmerVectors{D} method in src\vectors.jl
Base.hash(ka::AbstractKmerArray, h::UInt) = hash(get_S(ka), hash(ka.values, h))

function Base.:(==)(ka1::AbstractKmerArray, ka2::AbstractKmerArray)
   hash(ka1) == hash(ka2)
end
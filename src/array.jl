"""
    AbstractKmerCountArray{N, S, k, T <: Real, A <: AbstractArray{T, N}} <: AbstractArray{T, N}

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
2. `S`: Alphabet size is crucial. It's nonsensical to compare types with distinct alphabet sizes, even if k-mer sizes 
   and other attributes are the same.
3. `k`: k-mer size is the next important characteristic. Types with varying k-mer sizes shouldn't be directly compared 
   even if they share other attributes.
4. `T`: Element type might imply certain transformations or states (e.g., normalization when using floats). Types with 
   different element types can have distinct implications and hence shouldn't be equated.
5. `A`: Lastly, the specific array type (like `Vector` or `SparseVector`) is the most flexible comparison point. Types 
   might be comparable if they share all other parameters, even if their array types differ.

This order ensures that the most significant distinctions between types are made first, allowing for more meaningful
and nuanced comparisons.
"""
abstract type AbstractKmerCountArray{N, S, k, T <: Real, A <: AbstractArray{T, N}} <: AbstractArray{T, N} end

const AbstractKmerCountScalar = AbstractKmerCountArray{0}
const AbstractKmerCountVector = AbstractKmerCountArray{1}
const AbstractKmerCountMatrix = AbstractKmerCountArray{2}

@inline Base.size(kca::AbstractKmerCountArray) = size(kca.counts)
@inline Base.length(kca::AbstractKmerCountArray) = length(kca.counts)
@inline Base.getindex(kca::AbstractKmerCountArray, i) = kca.counts[i]
@inline Base.getindex(kca::AbstractKmerCountArray, i, j) = kca.counts[i, j]
@inline Base.setindex!(kca::AbstractKmerCountArray, v::Real, i) = kca.counts[i] = v
@inline Base.setindex!(kca::AbstractKmerCountArray, v::Real, i, j) = kca.counts[i, j] = v
@inline get_S(::AbstractKmerCountArray{N, S}) where {N, S} = S
@inline get_k(::AbstractKmerCountArray{N, S, k}) where {N, S, k} = k
@inline counts(kca::AbstractKmerCountArray) = kca.counts
@inline zeros!(kca::AbstractKmerCountArray) = fill!(kca.counts, zero(eltype(kca)))

Base.hash(kca::AbstractKmerCountArray, h::UInt) = hash(typeof(kca), hash(kca.counts, h))

function Base.:(==)(kca1::AbstractKmerCountArray, kca2::AbstractKmerCountArray)
    all((
        get_S(kca1) == get_S(kca2),
        kca1.counts == kca2.counts,
        typeof(kca1) == typeof(kca2) || typeof(kca1.counts) != typeof(kca2.counts), 
    ))
end
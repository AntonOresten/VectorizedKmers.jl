"""
    AbstractKmerCounts

type parameters are in order of importance: N, S, k, T, A
rationale:
even if S, k, T, A are the same, it doesn't make sense to directly compare types with different dimension N
even if k, T, A are the same, it doesn't make sense to directly compare types with different alphabet size S
even if T, A are the same, it doesn't make sense to directly compare types with different k-mer size k
even if array type can't be the same if T is different, it doesn't make sense to directly compare types with different
element type T because that might have some meaning like having been normalized if it's float or something
lastly: it does make sense to compare types with different array type A, like Vector and SparseVector,
if all other parameters are the same
"""
abstract type AbstractKmerCounts{N, S, k, T <: Real, A <: AbstractArray{T, N}} <: AbstractArray{T, N} end

@inline Base.getindex(kc::AbstractKmerCounts, i) = kc.counts[i]
@inline Base.getindex(kc::AbstractKmerCounts, i, j) = kc.counts[i, j]
@inline Base.setindex!(kc::AbstractKmerCounts, v::Real, i) = kc.counts[i] = v
@inline Base.setindex!(kc::AbstractKmerCounts, v::Real, i, j) = kc.counts[i, j] = v
@inline get_S(::AbstractKmerCounts{N, S}) where {N, S} = S
@inline get_k(::AbstractKmerCounts{N, S, k}) where {N, S, k} = k
@inline counts(kc::AbstractKmerCounts) = kc.counts
@inline zeros!(kc::AbstractKmerCounts) = fill!(kc.counts, zero(eltype(kc)))

Base.hash(kc::AbstractKmerCounts, h::UInt) = hash(typeof(kc), hash(kc.counts, h))

function Base.:(==)(kc1::AbstractKmerCounts, kc2::AbstractKmerCounts)
    all((
        get_S(kc1) == get_S(kc2),
        kc1.counts == kc2.counts,
        typeof(kc1) == typeof(kc2) || typeof(kc1.counts) != typeof(kc2.counts), 
    ))
end
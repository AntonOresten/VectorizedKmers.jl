"""
    KmerVectors

A container for k-mer vectors, where the values are stored together as columns or rows in a matrix.

!!! note
    It might be more efficient to store the values as columns in a matrix, since the elements in a column are contiguous in memory.

`D` is the dimension of the matrix in which the vectors are next to each other (1 for rows, 2 for columns).
`S` is the alphabet size,
`k` is the k-mer length,
and `T` is the element type of the underlying `values` field,
which in turn has type `V`.
"""
struct KmerVectors{D, S, k, T, M} <: AbstractKmerMatrix{S, k, T, M}
    values::M

    function KmerVectors{D, S, k}(values::M) where {D, S, k, T <: Real, M <: AbstractMatrix{T}}
        @assert D isa Integer && 1 <= D <= 2 "Type parameter `D` must be an integer between 1 and 2."
        @assert size(values, 3-D) == S^k
        new{D, S, k, T, M}(values)
    end

    function KmerVectors{D, S, k}(n::Integer; T::Type{<:Real}=Int, zeros::Function=zeros) where {D, S, k}
        dims = (D == 1) ? (n, S^k) : (S^k, n)
        KmerVectors{D, S, k}(zeros(T, dims))
    end
end

const KmerRows = KmerVectors{1}
const KmerColumns = KmerVectors{2}

# used to distinguish between KmerColumns{4, 2}(4) and KmerRows{4, 1}(16)
Base.hash(kvs::KmerVectors{D}, h::UInt) where D = hash(get_S(kvs), hash(kvs.values, h ⊻ D))

@inline Base.length(kvs::KmerVectors{D}) where D = size(kvs, D)

@inline function Base.getindex(krs::KmerRows{S, k}, i) where {S, k}
    KmerVector{S, k}(view(krs.values, i, :))
end

@inline function Base.getindex(kcs::KmerColumns{S, k}, i) where {S, k}
    KmerVector{S, k}(view(kcs.values, :, i))
end

@inline _eachvec_iter(kvs::KmerVectors{D}) where D = (D == 1) ? eachrow(kvs.values) : eachcol(kvs.values)
@inline eachvec(kvs::KmerVectors{D, S, k}) where {D, S, k} = (KmerVector{S, k}(v) for v in _eachvec_iter(kvs))
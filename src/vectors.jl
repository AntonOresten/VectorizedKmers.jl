"""
    KmerCountVectors

A container for k-mer counts, where k-mer counts are stored together as columns or rows in a matrix.
It is more efficient to store k-mer counts as columns in a matrix, since the elements in a column are contiguous in memory.

`D` is the dimension of the matrix in which the vectors are next to each other (1 for rows, 2 for columns).
`S` is the alphabet size,
`k` is the k-mer length,
and `T` is the element type of the underlying `counts` field,
which in turn has type `V`.
"""
struct KmerCountVectors{D, S, k, T, M} <: AbstractKmerCountMatrix{S, k, T, M}
    counts::M

    function KmerCountVectors{D, S, k}(counts::M) where {D, S, k, T <: Real, M <: AbstractMatrix{T}}
        @assert D isa Integer && 1 <= D <= 2 "Type parameter `D` must be an integer between 1 and 2."
        @assert size(counts, 3-D) == S^k
        new{D, S, k, T, M}(counts)
    end

    function KmerCountVectors{D, S, k}(n::Integer, T::Type{<:Real}=Int, zeros_func::Function=zeros) where {D, S, k}
        dims = (D == 1) ? (n, S^k) : (S^k, n)
        KmerCountVectors{D, S, k}(zeros_func(T, dims))
    end
end

const KmerCountRows = KmerCountVectors{1}
const KmerCountColumns = KmerCountVectors{2}

@inline Base.length(kcvs::KmerCountVectors{D}) where D = size(kcvs, D)

@inline function Base.getindex(kcr::KmerCountRows{S, k}, i) where {S, k}
    KmerCountVector{S, k}(view(kcr.counts, i, :))
end

@inline function Base.getindex(kcc::KmerCountColumns{S, k}, i) where {S, k}
    KmerCountVector{S, k}(view(kcc.counts, :, i))
end

@inline _eachvec_iter(kcvs::KmerCountVectors{D}) where D = (D == 1) ? eachrow(kcvs.counts) : eachcol(kcvs.counts)
@inline eachvec(kcvs::KmerCountVectors{D, S, k}) where {D, S, k} = (KmerCountVector{S, k}(v) for v in _eachvec_iter(kcvs))
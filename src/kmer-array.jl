import StaticArraysCore: StaticArray
import OffsetArrays: OffsetArray

ktuple(N::T, K::Int) where T = NTuple{K, T}(N for _ in 1:K)
offset_axes(N::Int, K::Int, offset::Int=-1) = ktuple(1+offset:N+offset, K)
hypercubify(A::AbstractArray, N::Int, K::Int) = reshape(A, ktuple(N, K))

"""
    KmerArray{N, K, T <: Real, A <: AbstractArray{T, K}} <: StaticArray{NTuple{K, N}, T, K}

- `N` is the alphabet size
- `K` is the K-mer size
- `T` is the element type
- `A` is the array type
"""
struct KmerArray{N, K, T <: Real, A <: AbstractArray{T, K}} <: StaticArray{NTuple{K, N}, T, K}
    offset_values::OffsetArray{T, K, A}

    function KmerArray{N, K, T, A}(offset_values::OffsetArray{T, K, A}) where {N, K, T <: Real, A <: AbstractArray{T, K}}
        all(isequal(N), ktuple(N, K)) || throw(ArgumentError("all dimensions of `offset_values` must be equal to `N`"))
        size(offset_values) == ktuple(N, K) || throw(ArgumentError("size(offset_values) must be $(ktuple(N, K))"))
        axes(offset_values) == offset_axes(N, K) || throw(ArgumentError("axes(offset_values) must be $(offset_axes(N, K))"))
        return new{N, K, T, A}(offset_values)
    end
end

KmerArray{N, K, T, A}(values::A) where {N, K, T <: Real, A <: AbstractArray{T, K}} = KmerArray{N, K, T, A}(OffsetArray(values, offset_axes(N, K)))

KmerArray{N, K}(values::A) where {N, K, T <: Real, A <: AbstractArray{T, K}} = KmerArray{N, K, T, A}(values)
KmerArray{N, K}(values::AbstractArray) where {N, K} = KmerArray{N, K}(hypercubify(values, N, K))

KmerArray(values::AbstractArray) = KmerArray{size(values, 1), ndims(values)}(values)
KmerArray(N::Int, K::Int, T::Type{<:Real}=Int, zeros::Function=zeros) = KmerArray{N, K}(zeros(T, ktuple(N, K)))

@inline Base.size(::KmerArray{N, K}) where {N, K} = ktuple(N, K)
@inline Base.length(::KmerArray{N, K}) where {N, K} = N^K
@inline Base.axes(::KmerArray{N, K}) where {N, K} = offset_axes(N, K)

@inline Base.values(ka::KmerArray) = ka.offset_values.parent

@inline Base.to_index(::KmerArray{N, K}, kmer::AbstractVector{<:Integer}) where {N, K} = LinearIndices(ktuple(N, K))[reverse(kmer .+ 1)...] - 1

@inline Base.getindex(ka::KmerArray, inds::Vararg{Integer, N}) where N = ka.offset_values[inds...]
@inline Base.getindex(ka::KmerArray, i::Int) = ka.offset_values.parent[i + 1] # need to access parent since linear indexing gets offset if K == 1
@inline Base.getindex(ka::KmerArray, i::Integer) = ka[i % Int]
@inline Base.getindex(ka::KmerArray, i::AbstractVector{<:Integer}) = ka[Base.to_index(ka, i)]

@inline Base.setindex!(ka::KmerArray, v, inds::Vararg{Integer, N}) where N = (ka.offset_values[inds...] = v)
@inline Base.setindex!(ka::KmerArray, v, i::Int) = (ka.offset_values.parent[i + 1] = v)
@inline Base.setindex!(ka::KmerArray, v, i::Integer) = (ka[i % Int] = v)
@inline Base.setindex!(ka::KmerArray, v, i::AbstractVector{<:Integer}) = (ka[Base.to_index(ka, i)] = v)

Base.show(io::IO, ka::KmerArray) = print(io, "$(typeof(ka)) with size $(size(ka))")
Base.show(io::IO, ::MIME"text/plain", ka::KmerArray) = show(io, ka)

zeros!(ka::KmerArray) = fill!(ka.offset_values, zero(eltype(ka)))
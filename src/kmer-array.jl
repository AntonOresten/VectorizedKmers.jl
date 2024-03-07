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
struct KmerArray{N, K, T <: Real, A <: AbstractArray{T, K}} <: AbstractArray{T, K}
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

@inline Base.values(ka::KmerArray) = ka.offset_values.parent

@inline Base.size(::KmerArray{N, K}) where {N, K} = ktuple(N, K)
@inline Base.length(::KmerArray{N, K}) where {N, K} = N^K
@inline Base.axes(::KmerArray{N, K}) where {N, K} = offset_axes(N, K)

@inline axis_index(::KmerArray, m::Integer) = m

@inline Base.getindex(ka::KmerArray, kmer::Integer) = ka.offset_values.parent[kmer + 1] # need to access parent since linear indexing gets offset if K == 1
@inline Base.getindex(ka::KmerArray{N, K}, axis_indices::CartesianIndex{K}) where {N, K} = ka.offset_values[axis_indices]
@inline Base.getindex(ka::KmerArray{N, K}, axis_indices::Vararg{Int, K}) where {N, K} = ka.offset_values[axis_indices...]
@inline Base.getindex(ka::KmerArray, kmer) = ka[(axis_index(ka, m) for m in Iterators.reverse(kmer))...]

@inline Base.setindex!(ka::KmerArray, v, kmer::Integer) = (ka.offset_values.parent[kmer + 1] = v)
@inline Base.setindex!(ka::KmerArray{N, K}, v, axis_indices::CartesianIndex{K}) where {N, K} = (ka.offset_values[axis_indices] = v)
@inline Base.setindex!(ka::KmerArray{N, K}, v, axis_indices::Vararg{Int, K}) where {N, K} = (ka.offset_values[axis_indices...] = v)
@inline Base.setindex!(ka::KmerArray, v, kmer) = (ka[(axis_index(ka, m) for m in Iterators.reverse(kmer))...] = v)

Base.similar(ka::KmerArray, ::Type{T}=eltype(ka), dims::Dims=size(ka)) where T = KmerArray(similar(ka.offset_values, T, dims))

Base.show(io::IO, ka::KmerArray) = print(io, "$(typeof(ka)) with size $(size(ka))")
Base.show(io::IO, ::MIME"text/plain", ka::KmerArray) = show(io, ka)

zeros!(ka::KmerArray) = fill!(ka.offset_values, zero(eltype(ka)))

default_alphabet_size(::Type{T}) where T = error("$(T) does not have a defined alphabet size. Please define `default_alphabet_size(::Type{<:$(T)})` or insert the alphabet size as a second argument in the `count_kmers` function call.")
default_alphabet_size(::T) where T = default_alphabet_size(T)
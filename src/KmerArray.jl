import StaticArraysCore: StaticArray

ktuple(N::T, K::Int) where T = NTuple{K, T}(N for _ in 1:K)

"""
    KmerArray{N, K, T <: Real, A <: AbstractArray{T, K}} <: StaticArray{NTuple{K, N}, T, K}

- `N` is the alphabet size
- `K` is the K-mer size
- `T` is the element type
- `A` is the array type
"""
struct KmerArray{N, K, T <: Real, A <: AbstractArray{T, K}} <: StaticArray{NTuple{K, N}, T, K}
    values::A

    function KmerArray{N, K, T, A}(values::A) where {N, K, T <: Real, A <: AbstractArray{T, K}}
        size(values) == ktuple(N, K) || throw(ArgumentError("size(values) must be $(ktuple(N, K))"))
        axes(values) == ktuple(1:N, K) || throw(ArgumentError("axes(values) must not be offset"))
        return new{N, K, T, A}(values)
    end
end

KmerArray{N, K}(values::A) where {N, K, T <: Real, A <: AbstractArray{T, K}} = KmerArray{N, K, T, A}(values)
KmerArray{N, K}(values::AbstractArray) where {N, K} = KmerArray{N, K}(reshape(values, ktuple(N, K)))
KmerArray{N}(values::AbstractArray) where N = KmerArray{N, Int(log(N, length(values)))}(values)

KmerArray(values::AbstractArray) = KmerArray{size(values, 1), ndims(values)}(values)
KmerArray(N::Int, K::Int, T::Type{<:Real}=Int, zeros::Function=zeros) = KmerArray{N, K}(zeros(T, ktuple(N, K)))

@inline Base.values(ka::KmerArray) = ka.values

@inline Base.size(::KmerArray{N, K}) where {N, K} = ktuple(N, K)
@inline Base.length(::KmerArray{N, K}) where {N, K} = N^K
@inline Base.axes(::KmerArray{N, K}) where {N, K} = ktuple(0:N-1, K)

@inline axis_index(::KmerArray, m::Integer) = m
@inline deconstruct(ka::KmerArray{N}, kmer) where N = (axis_index(ka, m) for m in Iterators.reverse(kmer))

@inline Base.getindex(ka::KmerArray, kmer::Integer) = ka.values[kmer + 1]
@inline Base.getindex(ka::KmerArray{N, K}, axis_indices::Vararg{Integer, K}) where {N, K} = ka.values[(axis_indices .+ 1)...]
@inline Base.getindex(ka::KmerArray{N, K}, axis_indices::CartesianIndex{K}) where {N, K} = ka[Tuple(axis_indices)...]
@inline Base.getindex(ka::KmerArray, kmer) = ka[deconstruct(ka, kmer)...]

@inline Base.setindex!(ka::KmerArray, v, kmer::Integer) = (ka.values[kmer + 1] = v)
@inline Base.setindex!(ka::KmerArray{N, K}, v, axis_indices::Vararg{Integer, K}) where {N, K} = (ka.values[(axis_indices .+ 1)...] = v)
@inline Base.setindex!(ka::KmerArray{N, K}, v, axis_indices::CartesianIndex{K}) where {N, K} = (ka[Tuple(axis_indices)...] = v)
@inline Base.setindex!(ka::KmerArray, v, kmer) = (ka[deconstruct(ka, kmer)...] = v)

# this changes the behavior of copy(ka) to return a KmerArray
#Base.similar(ka::KmerArray, ::Type{T}=eltype(ka), dims::Dims=size(ka)) where T = KmerArray(similar(ka.values, T, dims))

Base.show(io::IO, ka::KmerArray) = print(io, "$(typeof(ka)) with size $(size(ka))")
Base.show(io::IO, ::MIME"text/plain", ka::KmerArray) = show(io, ka)

zeros!(ka::KmerArray) = fill!(ka.values, zero(eltype(ka)))

default_alphabet_size(::Type{T}) where T = error("$(T) does not have a defined alphabet size. Please define `default_alphabet_size(::Type{<:$(T)})` or insert the alphabet size as a second argument in the `count_kmers` function call.")
default_alphabet_size(::T) where T = default_alphabet_size(T)
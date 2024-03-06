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
    values::OffsetArray{T, K, A}

    function KmerArray{N, K, T, A}(values::OffsetArray{T, K, A}) where {N, K, T <: Real, A <: AbstractArray{T, K}}
        all(isequal(N), ktuple(N, K)) || throw(ArgumentError("all dimensions of `values` must be equal to `N`"))
        size(values) == ktuple(N, K) || throw(ArgumentError("size(values) must be $(ktuple(N, K))"))
        return new{N, K, T, A}(values)
    end
end

KmerArray{N, K, T, A}(values::A) where {N, K, T <: Real, A <: AbstractArray{T, K}} = KmerArray{N, K, T, A}(OffsetArray(values, offset_axes(N, K)))

KmerArray{N, K}(values::A) where {N, K, T <: Real, A <: AbstractArray{T, K}} = KmerArray{N, K, T, A}(values)
KmerArray{N, K}(values::AbstractArray) where {N, K} = KmerArray{N, K}(hypercubify(values, N, K))

KmerArray(values::AbstractArray) = KmerArray{size(values, 1), ndims(values)}(values)
KmerArray(N::Int, K::Int, T::Type{<:Real}=Int, zeros::Function=zeros) = KmerArray{N, K}(zeros(T, ktuple(N, K)))

Base.size(ka::KmerArray) = size(ka.values)
Base.length(ka::KmerArray) = length(ka.values)

Base.axes(::KmerArray{N, K}) where {N, K} = offset_axes(N, K)

Base.getindex(ka::KmerArray, i::Vararg{<:Integer, N}) where N = ka.values[i...]
Base.getindex(ka::KmerArray, i::Integer) = ka.values.parent[i+1]

Base.setindex!(ka::KmerArray, v, i::Vararg{<:Integer, N}) where N = (ka.values[i...] = v)
Base.setindex!(ka::KmerArray, v, i::Integer) = (ka.values.parent[i+1] = v)

Base.show(io::IO, ka::KmerArray) = print(io, "$(typeof(ka)) with size $(size(ka.values))")
Base.show(io::IO, ::MIME"text/plain", ka::KmerArray) = show(io, ka)

zeros!(ka::KmerArray) = fill!(ka.values, zero(eltype(ka)))
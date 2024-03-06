import StaticArraysCore: StaticArray

sizetuple(N::Int, K::Int) = NTuple{K, Int}(N for _ in 1:K)
sizetuple(::Type{NTuple{K, N}}) where {K, N} = sizetuple(N, K)
hypercubify(A::AbstractArray, N::Int, K::Int) = reshape(A, sizetuple(N, K))

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
        all(isequal(N), sizetuple(N, K)) || throw(ArgumentError("all dimensions of `values` must be equal to `N`"))
        size(values) == sizetuple(N, K) || throw(ArgumentError("size(values) must be $(sizetuple(N, K))"))
        return new{N, K, T, A}(values)
    end
end

KmerArray{N, K}(values::A) where {N, K, T <: Real, A <: AbstractArray{T, K}} = KmerArray{N, K, T, A}(values)
KmerArray{N, K}(values::AbstractArray) where {N, K} = KmerArray{N, K}(hypercubify(values, N, K))

KmerArray(values::AbstractArray) = KmerArray{size(values, 1), ndims(values)}(values)
KmerArray(N::Int, K::Int, T::Type{<:Real}=Int, zeros::Function=zeros) = KmerArray{N, K}(zeros(T, sizetuple(N, K)))

Base.size(ka::KmerArray) = size(ka.values)
Base.length(ka::KmerArray) = length(ka.values)

Base.getindex(ka::KmerArray, i::Vararg{Int, N}) where N = ka.values[i...]
Base.getindex(ka::KmerArray, i) = ka.values[i]

Base.setindex!(ka::KmerArray, v, i::Vararg{Int, N}) where N = (ka.values[i...] = v)
Base.setindex!(ka::KmerArray, v, i) = (ka.values[i] = v)

Base.show(io::IO, ka::KmerArray) = print(io, "$(typeof(ka)) with size $(size(ka.values))")
Base.show(io::IO, ::MIME"text/plain", ka::KmerArray) = show(io, ka)

zeros!(ka::KmerArray) = fill!(ka.values, zero(eltype(ka)))
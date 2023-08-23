"""
    KmerCountColumns{S, k, T, M} <: AbstractKmerCountMatrix{S, k, T, M}

A container for k-mer counts, where k-mer counts are stored together as columns in a matrix.
This is more efficient than storing k-mer counts as rows in a matrix, since the elements in a column are contiguous in memory.
"""
struct KmerCountColumns{S, k, T, M} <: AbstractKmerCountMatrix{S, k, T, M}
    counts::M

    function KmerCountColumns{S, k}(counts::M) where {S, k, T <: Real, M <: AbstractMatrix{T}}
        @assert size(counts, 1) == S^k
        new{S, k, T, M}(counts)
    end

    function KmerCountColumns{S, k}(n::Integer, T::Type{<:Real}=Int, zeros_func::Function=zeros) where {S, k}
        KmerCountColumns{S, k}(zeros_func(T, S^k, n))
    end
end

@inline Base.length(kcc::KmerCountColumns) = size(kcc, 2)

# overrides getindex(::AbstractKmerCounts, ::Any) from kmer_count.jl
@inline function Base.getindex(kcc::KmerCountColumns{S, k}, i::Integer) where {S, k}
    KmerCountVector{S, k}(view(kcc.counts, :, i))
end
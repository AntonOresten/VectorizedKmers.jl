"""
    KmerCountRows{S, k, T, M} <: AbstractKmerCountMatrix{S, k, T, M}

A container for k-mer counts, where k-mer counts are stored together as rows in a matrix.
This is not as efficient as storing k-mer counts as columns in a matrix, since the elements in a row are not contiguous in memory.
"""
struct KmerCountRows{S, k, T, M} <: AbstractKmerCountMatrix{S, k, T, M}
    counts::M

    function KmerCountRows{S, k}(counts::M) where {S, k, T <: Real, M <: AbstractMatrix{T}}
        @assert size(counts, 2) == S^k
        new{S, k, T, M}(counts)
    end

    function KmerCountRows{S, k}(n::Integer, T::Type{<:Real}=Int, zeros_func::Function=zeros) where {S, k}
        KmerCountRows{S, k}(zeros_func(T, n, S^k))
    end
end

@inline Base.length(kcr::KmerCountRows) = size(kcr, 1)

# overrides getindex(::AbstractKmerCounts, ::Any) from kmer_count.jl
@inline function Base.getindex(kcr::KmerCountRows{S, k}, i::Integer) where {S, k}
    KmerCountVector{S, k}(view(kcr.counts, i, :))
end
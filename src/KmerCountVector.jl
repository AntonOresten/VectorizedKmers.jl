abstract type AbstractKmerCountVector{A, k, T, M <: AbstractMatrix{T}} <: AbstractVector{KmerCount{A, k, T}} end

struct KmerCountColumns{A, k, T, M} <: AbstractKmerCountVector{A, k, T, M}
    data::M

    function KmerCountColumns{A, k, T}(data::M) where {A, k, T, M <: AbstractMatrix{T}}
        @assert size(data, 1) == A^k
        new{A, k, T, M}(data)
    end
end

@inline Base.size(kmer_count::KmerCountColumns) = (size(kmer_count.data, 2),)
@inline Base.length(kmer_count::KmerCountColumns) = size(kmer_count.data, 2)
@inline Base.getindex(kmer_count::KmerCountColumns{A, k, T}, i::Integer) where {A, k, T} = KmerCount{A, k, T}(view(kmer_count.data, :, i))

struct KmerCountRows{A, k, T, M} <: AbstractKmerCountVector{A, k, T, M}
    data::M

    function KmerCountRows{A, k, T}(data::M) where {A, k, T, M <: AbstractMatrix{T}}
        @assert size(data, 2) == A^k
        new{A, k, T, M}(data)
    end
end

@inline Base.size(kmer_count::KmerCountRows) = (size(kmer_count.data, 1),)
@inline Base.length(kmer_count::KmerCountRows) = size(kmer_count.data, 1)
@inline Base.getindex(kmer_count::KmerCountRows{A, k, T}, i::Integer) where {A, k, T} = KmerCount{A, k, T}(view(kmer_count.data, i, :))

abstract type AbstractKmerCount{A, k, T <: Real} <: AbstractVector{T} end

@inline Base.size(kmer_count::AbstractKmerCount) = size(kmer_count.data)
@inline Base.length(kmer_count::AbstractKmerCount) = length(kmer_count.data)
@inline Base.getindex(kmer_count::AbstractKmerCount, i::Integer) = kmer_count.data[i]
@inline Base.setindex!(kmer_count::AbstractKmerCount, v::Real, i::Integer) = kmer_count.data[i] = v

struct KmerCount{A, k, T, ArrayType <: AbstractArray} <: AbstractKmerCount{A, k, T}
    data::AbstractVector{T}

    function KmerCount{A, k, T}(data::AbstractVector{T}) where {A, k, T}
        new{A, k, T, typeof(data)}(data)
    end

    function KmerCount{A, k, T}() where {A, k, T}
        KmerCount{A, k, T}(zeros(T, A^k))
    end
end

function count_kmers!(kmer_count::KmerCount{A, k, T}, kmers::Vector{<:Integer}; reset::Bool=true) where {A, k, T}
    reset && fill!(kmer_count, zero(T))
    for kmer in kmers
        kmer_count[kmer + 1] += 1
    end
    kmer_count
end
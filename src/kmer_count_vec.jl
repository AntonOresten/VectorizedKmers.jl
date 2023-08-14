"""
    AbstractKmerCountVector{A, K, T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{KmerCount{A, K, T, V} where {V <: AbstractVector{T}}}

Yup... that's indeed an abomination of a type.
A container for K-mer counts, where K-mer counts are stored together as rows or columns in a matrix.
`A` is the alphabet size, `K` is the K-mer size, `T` is the element type of the counts,
and `M` is the type of the matrix in which the K-mer counts are stored.
"""
abstract type AbstractKmerCountVector{A, K, T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{KmerCount{A, K, T, V} where {V <: AbstractVector{T}}}  end


"""
    KmerCountColumns{A, K, T, M} <: AbstractKmerCountVector{A, K, T, M}

A container for K-mer counts, where K-mer counts are stored together as columns in a matrix.
This is more efficient than storing K-mer counts as rows in a matrix, since the elements in a column are contiguous in memory.
"""
struct KmerCountColumns{A, K, T, M} <: AbstractKmerCountVector{A, K, T, M}
    counts::M

    function KmerCountColumns{A, K}(counts::M) where {A, K, T <: Real, M <: AbstractMatrix{T}}
        @assert size(counts, 1) == A^K
        new{A, K, T, M}(counts)
    end
end

@inline Base.size(kmer_count::KmerCountColumns) = (size(kmer_count.counts, 2),)
@inline Base.length(kmer_count::KmerCountColumns) = size(kmer_count.counts, 2)

@inline function Base.getindex(kmer_count::KmerCountColumns{A, K}, i::Integer) where {A, K}
    KmerCount{A, K}(view(kmer_count.counts, :, i))
end


"""
    KmerCountRows{A, K, T, M} <: AbstractKmerCountVector{A, K, T, M}

A container for K-mer counts, where K-mer counts are stored together as rows in a matrix.
This is not as efficient as storing K-mer counts as columns in a matrix, since the elements in a row are not contiguous in memory.
"""
struct KmerCountRows{A, K, T, M} <: AbstractKmerCountVector{A, K, T, M}
    counts::M

    function KmerCountRows{A, K}(counts::M) where {A, K, T <: Real, M <: AbstractMatrix{T}}
        @assert size(counts, 2) == A^K
        new{A, K, T, M}(counts)
    end
end

@inline Base.size(kmer_count::KmerCountRows) = (size(kmer_count.counts, 1),)
@inline Base.length(kmer_count::KmerCountRows) = size(kmer_count.counts, 1)

@inline function Base.getindex(kmer_count::KmerCountRows{A, K}, i::Integer) where {A, K}
    KmerCount{A, K}(view(kmer_count.counts, i, :))
end


function KmerCountColumns(kcc::KmerCountRows{A, K}) where {A, K}
    KmerCountColumns{A, K}(transpose(kcc.counts))
end

function KmerCountRows(kcc::KmerCountColumns{A, K}) where {A, K}
    KmerCountRows{A, K}(transpose(kcc.counts))
end

@inline Base.transpose(kcc::KmerCountColumns) = KmerCountRows(kcc)
@inline Base.transpose(kcr::KmerCountRows) = KmerCountColumns(kcr)

@inline Base.adjoint(kcv::AbstractKmerCountVector) = transpose(kcv)
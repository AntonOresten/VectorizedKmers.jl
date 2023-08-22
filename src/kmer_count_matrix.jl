"""
    AbstractKmerCountMatrix{A, k, T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{KmerCountVector{A, k, T, V} where {V <: AbstractVector{T}}}

Yup... that's indeed an abomination of a type.
A container for k-mer counts, where k-mer counts are stored together as rows or columns in a matrix.
`A` is the alphabet size, `k` is the k-mer size, `T` is the element type of the counts,
and `M` is the type of the matrix in which the k-mer counts are stored.
"""
abstract type AbstractKmerCountMatrix{A, k, T <: Real, M <: AbstractMatrix{T}} <: AbstractVector{KmerCountVector{A, k, T, V} where {V <: AbstractVector{T}}}  end

@inline reset!(kcm::AbstractKmerCountMatrix) = fill!(kcm.counts, 0)

"""
    KmerCountColumns{A, k, T, M} <: AbstractKmerCountMatrix{A, k, T, M}

A container for k-mer counts, where k-mer counts are stored together as columns in a matrix.
This is more efficient than storing k-mer counts as rows in a matrix, since the elements in a column are contiguous in memory.
"""
struct KmerCountColumns{A, k, T, M} <: AbstractKmerCountMatrix{A, k, T, M}
    counts::M

    function KmerCountColumns{A, k}(counts::M) where {A, k, T <: Real, M <: AbstractMatrix{T}}
        @assert size(counts, 1) == A^k
        new{A, k, T, M}(counts)
    end
end

@inline Base.size(kmer_count::KmerCountColumns) = (size(kmer_count.counts, 2),)
@inline Base.length(kmer_count::KmerCountColumns) = size(kmer_count.counts, 2)

@inline function Base.getindex(kmer_count::KmerCountColumns{A, k}, i::Integer) where {A, k}
    KmerCountVector{A, k}(view(kmer_count.counts, :, i))
end


"""
    KmerCountRows{A, k, T, M} <: AbstractKmerCountMatrix{A, k, T, M}

A container for k-mer counts, where k-mer counts are stored together as rows in a matrix.
This is not as efficient as storing k-mer counts as columns in a matrix, since the elements in a row are not contiguous in memory.
"""
struct KmerCountRows{A, k, T, M} <: AbstractKmerCountMatrix{A, k, T, M}
    counts::M

    function KmerCountRows{A, k}(counts::M) where {A, k, T <: Real, M <: AbstractMatrix{T}}
        @assert size(counts, 2) == A^k
        new{A, k, T, M}(counts)
    end
end

@inline Base.size(kmer_count::KmerCountRows) = (size(kmer_count.counts, 1),)
@inline Base.length(kmer_count::KmerCountRows) = size(kmer_count.counts, 1)

@inline function Base.getindex(kmer_count::KmerCountRows{A, k}, i::Integer) where {A, k}
    KmerCountVector{A, k}(view(kmer_count.counts, i, :))
end


function KmerCountColumns(kcc::KmerCountRows{A, k}) where {A, k}
    KmerCountColumns{A, k}(transpose(kcc.counts))
end

function KmerCountRows(kcc::KmerCountColumns{A, k}) where {A, k}
    KmerCountRows{A, k}(transpose(kcc.counts))
end

@inline Base.transpose(kcc::KmerCountColumns) = KmerCountRows(kcc)
@inline Base.transpose(kcr::KmerCountRows) = KmerCountColumns(kcr)

@inline Base.adjoint(kcv::AbstractKmerCountMatrix) = transpose(kcv)
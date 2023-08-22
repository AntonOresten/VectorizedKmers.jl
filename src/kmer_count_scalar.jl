# S and k gives context to the `kmer` field
struct KmerCountScalar{S, k, T, A} <: AbstractKmerCounts{0, S, k, T, A}
    kmer::Integer
    counts::A

    function KmerCountScalar{S, k}(kmer::Integer, counts::A) where {S, k, T, A <: AbstractArray{T, 0}}
        new{S, k, T, A}(kmer, counts)
    end

    function KmerCountScalar{S, k, T, A}(kmer::Integer, count::T) where {S, k, T, A}
        new{S, k, T, A}(kmer, fill(count))
    end
end

@inline Base.size(::KmerCountScalar) = ()
@inline Base.length(::KmerCountScalar) = 1
@inline Base.getindex(kcs::KmerCountScalar) = kcs.counts[1]

@inline Base.repr(kcs::KmerCountScalar) = kcs[]
@inline Base.show(io::IO, kcs::KmerCountScalar) = print(io, repr(kcs))
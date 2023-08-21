"""
    KmerCountVector{A, k, T, V}

`A` is the alphabet size,
`k` is the k-mer size,
and `T` is the element type of the underlying `counts` field,
which in turn has type `V`.
"""
struct KmerCountVector{A, k, T <: Real, V <: AbstractVector{T}} <: AbstractVector{T}
    counts::V

    function KmerCountVector{A, k, T, V}(counts::V) where {A, k, T, V}
        @assert length(counts) == A^k
        new{A, k, T, V}(counts)
    end

    function KmerCountVector{A, k}(counts::V) where {A, k, T <: Real, V <: AbstractVector{T}}
        KmerCountVector{A, k, T, V}(counts)
    end

    function KmerCountVector{A, k, T}(zeros_func::Function = zeros) where {A, k, T}
        KmerCountVector{A, k}(zeros_func(T, A^k))
    end

    function KmerCountVector{A, k}(zeros_func::Function = zeros) where {A, k}
        KmerCountVector{A, k, Int}(zeros_func)
    end
end

@inline Base.size(kmer_count::KmerCountVector) = size(kmer_count.counts)
@inline Base.length(kmer_count::KmerCountVector) = length(kmer_count.counts)
@inline Base.getindex(kmer_count::KmerCountVector, i::Integer) = kmer_count.counts[i]
@inline Base.setindex!(kmer_count::KmerCountVector, v::Real, i::Integer) = kmer_count.counts[i] = v
@inline Base.eltype(::KmerCountVector{A, k, T}) where {A, k, T} = T
@inline get_A(::KmerCountVector{A}) where A = A
@inline get_k(::KmerCountVector{A, k}) where {A, k} = k
@inline reset!(kmer_count::KmerCountVector) = fill!(kmer_count.counts, 0)

function Base.summary(kc::KmerCountVector)
    string(typeof(kc))
end

function Base.show(io::IO, kc::KmerCountVector)
    print(io, summary(kc)*"($(kc.counts))")
end

"""
    count_kmers!(kmer_count, kmers; reset=true)

Mutate the `counts` vector in `kmer_count` by counting k-mers in `kmers`.
The k-mers in `kmers` must be represented as integers between 0 and length(kmer_count) - 1.

If `reset` is `true`, the `counts` vector will be zero-ed before counting.
"""
function count_kmers!(
    kmer_count::KmerCountVector{A, k},
    kmers::Vector{<:Integer};
    reset::Bool = true,
) where {A, k}
    @assert maximum(kmers) < A^k
    reset && reset!(kmer_count)
    for kmer in kmers
        kmer_count[kmer + 1] += 1
    end
    kmer_count
end
#=  This is not a very efficient method, since it takes an entire vector. It is mainly used for testing.
    Ideally the k-mers would be procedurally calculated in constant memory. =#

"""
    count_kmers(KmerCountVector{A, k, T}, kmers; zeros_func=zeros)

Create a new A^k sized vector using `zeros_func` and count the k-mers in `kmers`.
The k-mers in `kmers` must be represented as integers between `0` and `length(kmer_count) - 1`.

If `reset` is `true`, the `counts` vector will be zero-ed before counting.
"""
function count_kmers(
    ::Type{KmerCountVector{A, k, T}}, kmers::Vector{<:Integer};
    zeros_func::Function = zeros,
) where {A, k, T}
    kmer_count = KmerCountVector{A, k, T}(zeros_func)
    count_kmers!(kmer_count, kmers, reset=false)
    return kmer_count
end
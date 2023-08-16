"""
    KmerCount{A, K, T, V}

`A` is the alphabet size,
`K` is the K-mer size,
and `T` is the element type of the underlying `counts` field,
which in turn has type `V`.
"""
struct KmerCount{A, K, T <: Real, V <: AbstractVector{T}} <: AbstractVector{T}
    counts::V

    function KmerCount{A, K, T, V}(counts::V) where {A, K, T, V}
        @assert length(counts) == A^K
        new{A, K, T, V}(counts)
    end

    function KmerCount{A, K}(counts::V) where {A, K, T <: Real, V <: AbstractVector{T}}
        KmerCount{A, K, T, V}(counts)
    end

    function KmerCount{A, K, T}(zeros_func::Function = zeros) where {A, K, T}
        KmerCount{A, K}(zeros_func(T, A^K))
    end

    function KmerCount{A, K}(zeros_func::Function = zeros) where {A, K}
        KmerCount{A, K, Int}(zeros_func)
    end
end

@inline Base.size(kmer_count::KmerCount) = size(kmer_count.counts)
@inline Base.length(kmer_count::KmerCount) = length(kmer_count.counts)
@inline Base.getindex(kmer_count::KmerCount, i::Integer) = kmer_count.counts[i]
@inline Base.setindex!(kmer_count::KmerCount, v::Real, i::Integer) = kmer_count.counts[i] = v
@inline Base.eltype(::KmerCount{A, K, T}) where {A, K, T} = T
@inline get_A(::KmerCount{A}) where A = A
@inline get_K(::KmerCount{A, K}) where {A, K} = K
@inline reset!(kmer_count::KmerCount) = fill!(kmer_count.counts, 0)

function Base.summary(kc::KmerCount)
    string(typeof(kc))
end

function Base.show(io::IO, kc::KmerCount)
    print(io, summary(kc)*"($(kc.counts))")
end

"""
    count_kmers!(kmer_count, kmers; reset=true)

Mutate the `counts` vector in `kmer_count` by counting K-mers in `kmers`.
The K-mers in `kmers` must be represented as integers between 0 and length(kmer_count) - 1.

If `reset` is `true`, the `counts` vector will be zero-ed before counting.
"""
function count_kmers!(
    kmer_count::KmerCount,
    kmers::Vector{<:Integer};
    reset::Bool = true,
)
    reset && reset!(kmer_count)
    for kmer in kmers
        kmer_count[kmer + 1] += 1
    end
    kmer_count
end
#=  This is not a very efficient method, since it takes an entire vector. It is mainly used for testing.
    Ideally the K-mers would be procedurally calculated in constant memory. =#

"""
    count_kmers(KmerCount{A, K, T}, kmers; zeros_func=zeros)

Create a new A^K sized vector using `zeros_func` and count the K-mers in `kmers`.
The K-mers in `kmers` must be represented as integers between `0` and `length(kmer_count) - 1`.

If `reset` is `true`, the `counts` vector will be zero-ed before counting.
"""
function count_kmers(
    ::Type{KmerCount{A, K, T}}, kmers::Vector{<:Integer};
    zeros_func::Function = zeros,
) where {A, K, T}
    kmer_count = KmerCount{A, K, T}(zeros_func)
    count_kmers!(kmer_count, kmers, reset=false)
    return kmer_count
end
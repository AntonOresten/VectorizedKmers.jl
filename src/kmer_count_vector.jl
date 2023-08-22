"""
    KmerCountVector{S, k, T, V}

`A` is the alphabet size,
`k` is the k-mer size,
and `T` is the element type of the underlying `counts` field,
which in turn has type `V`.
"""
struct KmerCountVector{S, k, T, V} <: AbstractKmerCounts{1, S, k, T, V}
    counts::V

    function KmerCountVector{S, k, T, V}(counts::V) where {S, k, T, V}
        @assert length(counts) == S^k
        new{S, k, T, V}(counts)
    end

    function KmerCountVector{S, k}(counts::V) where {S, k, T <: Real, V <: AbstractVector{T}}
        KmerCountVector{S, k, T, V}(counts)
    end

    function KmerCountVector{S, k, T}(zeros_func::Function = zeros) where {S, k, T}
        KmerCountVector{S, k}(zeros_func(T, S^k))
    end

    function KmerCountVector{S, k}(zeros_func::Function = zeros) where {S, k}
        KmerCountVector{S, k, Int}(zeros_func)
    end
end

@inline Base.size(kmer_count::KmerCountVector) = size(kmer_count.counts)
@inline Base.length(kmer_count::KmerCountVector) = length(kmer_count.counts)
@inline Base.getindex(kmer_count::KmerCountVector, i::Integer) = kmer_count.counts[i]

function Base.summary(kc::KmerCountVector)
    string(typeof(kc))
end

function Base.show(io::IO, kc::KmerCountVector)
    max_elements_to_show = 10
    len_counts = length(kc.counts)
    
    print(io, summary(kc) * "([")

    if len_counts <= max_elements_to_show
        print(io, join(kc.counts, ", "))
    else
        print(io, join(kc.counts[1:max_elements_to_show], ", "), " … (", len_counts - max_elements_to_show, " more)")
    end

    print(io, "])")
end

# The k-mers in `kmers` must be represented as integers between 0 and length(kmer_count) - 1.
#=  This is not a very efficient method, since it takes an entire vector. It is mainly used for testing.
    Ideally the k-mers would be procedurally calculated in constant memory. =#
function count_kmers!(
    kmer_count::KmerCountVector{S, k},
    kmers::Vector{<:Integer};
    reset::Bool = true,
) where {S, k}
    @assert maximum(kmers) < S^k
    reset && zeros!(kmer_count)
    for kmer in kmers
        kmer_count[kmer + 1] += 1
    end
    kmer_count
end

function count_kmers(
    ::Type{KmerCountVector{S, k, T}}, kmers::Vector{<:Integer};
    zeros_func::Function = zeros,
) where {S, k, T}
    kmer_count = KmerCountVector{S, k, T}(zeros_func)
    count_kmers!(kmer_count, kmers, reset=false)
    return kmer_count
end
"""
    KmerCountVector{S, k, T, V} <: AbstractKmerCountVector{S, k, T, V}

`A` is the alphabet size,
`k` is the k-mer size,
and `T` is the element type of the underlying `counts` field,
which in turn has type `V`.
"""
struct KmerCountVector{S, k, T, V} <: AbstractKmerCountVector{S, k, T, V}
    counts::V

    function KmerCountVector{S, k}(counts::V) where {S, k, T <: Real, V <: AbstractVector{T}}
        @assert length(counts) == S^k
        new{S, k, T, V}(counts)
    end

    function KmerCountVector{S, k}(T::Type{<:Real}=Int, zeros::Function=zeros) where {S, k}
        KmerCountVector{S, k}(zeros(T, S^k))
    end
end

function Base.summary(kcv::KmerCountVector)
    string(typeof(kcv))
end

function Base.show(io::IO, kcv::KmerCountVector)
    max_elements_to_show = 10
    len_counts = length(kcv.counts)
    
    print(io, summary(kcv) * "([")

    if len_counts <= max_elements_to_show
        print(io, join(kcv.counts, ", "))
    else
        print(io, join(view(kcv.counts, 1:max_elements_to_show), ", "), " … (", len_counts - max_elements_to_show, " more)")
    end

    print(io, "])")
end

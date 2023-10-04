"""
    KmerVector{S, k, T, V} <: AbstractKmerVector{S, k, T, V}

`A` is the alphabet size,
`k` is the k-mer size,
and `T` is the element type of the underlying `values` field,
which in turn has type `V`.
"""
struct KmerVector{S, k, T, V} <: AbstractKmerVector{S, k, T, V}
    values::V

    function KmerVector{S, k}(values::V) where {S, k, T <: Real, V <: AbstractVector{T}}
        @assert length(values) == S^k
        return new{S, k, T, V}(values)
    end

    function KmerVector{S, k}(; T::Type{<:Real}=Int, zeros::Function=zeros) where {S, k}
        return KmerVector{S, k}(zeros(T, S^k))
    end
end

Base.summary(kv::KmerVector) = string(typeof(kv))

function Base.show(io::IO, kv::KmerVector)
    max_elements_to_show = 10
    len_values = length(kv.values)
    
    print(io, summary(kv) * "([")

    if len_values <= max_elements_to_show
        print(io, join(kv.values, ", "))
    else
        print(io, join(view(kv.values, 1:max_elements_to_show), ", "), " … (", len_values - max_elements_to_show, " more)")
    end

    print(io, "])")

    return nothing
end

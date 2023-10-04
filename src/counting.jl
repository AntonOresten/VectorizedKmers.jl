"""
    count_kmers!(kmer_vector, sequence; reset=true)

Mutates `kmer_vector`.

If `reset` is `true`, the array will be zero-ed before counting.
"""
function count_kmers! end

# This method is used for testing.
# sequence is wrapped by RefValue cause otherwise it gets mistaken for a vector of sequences.
@inline function count_kmers!(
    kv::KmerVector{S, k}, sequence::Base.RefValue{Vector{T}};
    reset::Bool = true,
) where {S, k, T <: Integer}
    sequence = sequence[]
    @assert 0 <= maximum(sequence) < S
    reset && zeros!(kv)
    length(sequence) < k && return kv
    mask = S^k
    kmer = zero(T)
    for m in @view sequence[1:k-1]
        kmer = kmer * S + m
    end
    for m in @view sequence[k:end]
        kmer = kmer * S % mask + m
        kv.values[kmer + 1] += 1
    end
    return kv
end

@inline function count_kmers!(
    kvs::KmerVectors{D, S, k}, sequences::Vector{SequenceType};
    offset::Integer = 0, reset::Bool = true,
) where {D, S, k, SequenceType}
    kv_gen = Iterators.drop(eachvec(kvs), offset)
    for (kv, sequence) in zip(kv_gen, sequences)
        count_kmers!(kv, sequence, reset=reset)
    end
    return kvs
end


"""
    count_kmers(sequence, S, k; T=Int, zeros=zeros)

Creates a new S^k-sized vector using the supplied `zeros` function and counts the k-mers in `sequence`.
Since S is the alphabet, and the elements in sequence are integers, the maximum value in `sequence` must be less than `S`
"""
function count_kmers end

@inline function count_kmers(
    sequence::SequenceType, S::Integer, k::Integer;
    T::Type{<:Real} = Int, zeros::Function = zeros,
) where SequenceType
    kv = KmerVector{S, k}(T=T, zeros=zeros)
    count_kmers!(kv, sequence, reset=false)
    return kv
end

@inline function count_kmers(
    sequences::Vector{SequenceType}, S::Integer, k::Integer;
    T::Type{<:Real} = Int, zeros::Function = zeros, D = 2,
) where SequenceType
    kvs = KmerVectors{D, S, k}(length(sequences), T=T, zeros=zeros)
    count_kmers!(kvs, sequences, reset=false)
    return kvs
end


function alphabet_size(T)
    error("$(T) does not have a defined alphabet size. Please define `alphabet_size(::Type{<:$(T)})` or insert the alphabet size as a second argument in the `count_kmers` function call.")
end

@inline function count_kmers(
    sequence::SequenceType, k::Integer;
    T::Type{<:Real} = Int, zeros::Function = zeros,
) where SequenceType
    S = alphabet_size(SequenceType)
    return count_kmers(sequence, S, k, T=T, zeros=zeros)
end

@inline function count_kmers(
    sequences::Vector{SequenceType}, k::Integer;
    T::Type{<:Real} = Int, zeros::Function = zeros,
) where SequenceType
    S = alphabet_size(SequenceType)
    return count_kmers(sequences, S, k, T=T, zeros=zeros)
end

# for backward compatibility
@inline count_kmers(seq, k::Integer, T::Type{<:Real}, zeros::Function=zeros) = count_kmers(seq, k, T=T, zeros=zeros)
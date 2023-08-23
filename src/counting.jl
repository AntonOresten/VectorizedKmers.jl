"""
    count_kmers!(kmer_count_vector, sequence)

Mutates `kmer_count_vector`.

If `reset` is `true`, the array will be zero-ed before counting.
"""
function count_kmers! end

# This method is used for testing.
# sequence is wrapped by RefValue cause otherwise it gets mistaken for a vector of sequences.
function count_kmers!(
    kcv::KmerCountVector{S, k}, sequence::Base.RefValue{Vector{T}};
    reset::Bool = true,
) where {S, k, T <: Integer}
    sequence = sequence[]
    @assert 0 <= maximum(sequence) < S
    reset && zeros!(kcv)
    length(sequence) < k && return kcv
    counts = kcv.counts
    mask = unsigned(S^k)

    kmer = zero(UInt)
    for m in view(sequence, 1:k-1)
        kmer = kmer * S + m
    end

    for m in view(sequence, k:length(sequence))
        kmer = kmer * S % mask + m
        counts[kmer + 1] += 1
    end

    kcv
end

function count_kmers!(
    kcc::KmerCountVectors{D, S, k}, sequences::Vector{SequenceType};
    offset::Integer = 0, reset::Bool = true,
) where {D, S, k, SequenceType}
    kcv_gen = Iterators.drop(eachvec(kcc), offset)
    for (kcv, sequence) in zip(kcv_gen, sequences)
        count_kmers!(kcv, sequence, reset=reset)
    end
    kcc
end


"""
    count_kmers(sequence::Vector{<:Integer}, S, k; zeros=zeros, reset=true)

This is the default k-mer counting method.
It create a new S^k-sized vector using `zeros` and counts the k-mers in `sequence`.
Since S is the alphabet, and the elements in sequence are integers, the maximum value in `sequence` must be less than `S`
"""
function count_kmers end

function count_kmers(
    sequence::SequenceType, S::Integer, k::Integer;
    T::Type{<:Real} = Int, zeros::Function = zeros,
) where SequenceType
    kcv = KmerCountVector{S, k}(T=T, zeros=zeros)
    count_kmers!(kcv, sequence, reset=false)
    kcv
end

function count_kmers(
    sequences::Vector{SequenceType}, S::Integer, k::Integer;
    T::Type{<:Real} = Int, zeros::Function = zeros, D = 2,
) where SequenceType
    kcvs = KmerCountVectors{D, S, k}(length(sequences), T=T, zeros=zeros)
    count_kmers!(kcvs, sequences, reset=false)
    kcvs
end

alphabet_size(T::DataType) = error("$(T) does not have a defined alphabet size. Please define `alphabet_size(::Type{<:$(T)})` or insert the alphabet size as a second argument in the `count_kmers` function call.")

function count_kmers(
    sequence::SequenceType, k::Integer;
    T::Type{<:Real} = Int, zeros::Function = zeros,
) where SequenceType
    S = alphabet_size(SequenceType)
    count_kmers(sequence, S, k, T=T, zeros=zeros)
end

function count_kmers(
    sequences::Vector{SequenceType}, k::Integer;
    T::Type{<:Real} = Int, zeros::Function = zeros,
) where SequenceType
    S = alphabet_size(SequenceType)
    count_kmers(sequences, S, k, T=T, zeros=zeros)
end
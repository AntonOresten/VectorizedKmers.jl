```@meta
CurrentModule = VectorizedKmers
DocTestSetup = quote
    using VectorizedKmers
end
```

# Vectors of k-mer counts

The `KmerVector` type has four type parameters, but you only really need to care about the first two: `S`, the alphabet size, and `k`, the k-mer length. So, to count the 6-mers of a DNA sequence, you would use `KmerVector{4, 6}`. For each of these k-mer counts, memory for a vector of size `S^k` is allocated, unless a vector type like `SparseVector` is used. This brings us to the two other type parameters: `T`, which is the vector element type, and `V`, which is the type of the actual vector.

Let's see it in action! Here we import `BioSequences` to unlock a method of `count_kmers` that works on the `LongDNA` type. In this example, we count the 1-mers of the sequence `GATTACA`. The result is a `KmerVector{4, 1, Int64, Vector{Int64}}`, which is a vector of 4 `Int64` elements.

```jldoctest
julia> using BioSequences

julia> kv = count_kmers(dna"GATTACA", 1)
4-element KmerVector{4, 1, Int64, Vector{Int64}}:
 3
 1
 1
 2
```

We can also change the element type to save memory, at the cost of a smaller maximum count value.

```jldoctest
julia> count_kmers(dna"GATTACA", 1, T=UInt16)
4-element KmerVector{4, 1, UInt16, Vector{UInt16}}:
 0x0003
 0x0001
 0x0001
 0x0002
```

!!! note
    Be careful when using element types with fewer bits, such as `UInt16`. Even though you might not expect any one k-mer to occur more than 65,535 times, some vector operations such as `LinearAlgebra.dot` and `Distances.sqeuclidean` will still use `UInt16` when summing up terms, which might lead to integer overflow.

The `count_kmers!` function will modify a `KmerVector` instance in-place. This is useful when you want to count k-mers in a sequence without allocating a new vector.

Let's create a `KmerVector` instance with `SparseVector` as its vector type. We can do this by passing `spzeros` to the `KmerVector{S, k, T}` constructor (it's just `zeros` by default):

```jldoctest
julia> using SparseArrays

julia> kv2 = KmerVector{4, 1}(T=Int64, zeros=spzeros)
4-element KmerVector{4, 1, Int64, SparseVector{Int64, Int64}}:
 0
 0
 0
 0

julia> count_kmers!(kv2, dna"GATTACA")
4-element KmerVector{4, 1, Int64, SparseVector{Int64, Int64}}:
 3
 1
 1
 2

julia> kv == kv2
```
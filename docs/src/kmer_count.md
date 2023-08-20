```@meta
CurrentModule = VectorizedKmers
DocTestSetup = quote
    using VectorizedKmers
end
```

# Counting k-mers with the `KmerCountVector` type

The `KmerCountVector` type has four type parameters, but you only really need to care about the first two: `A`, the alphabet size, and `k`, the k-mer length. So, to count the 6-mers of a DNA sequence, you would use `KmerCountVector{4, 6}`. For each of these k-mer counts, memory for a vector of size `A^k` is allocated, unless a vector type like `SparseVector` is used. This brings us to the two other type parameters: `T`, which is the vector element type, and `V`, which is the type of the actual vector.

Let's see it in action! Here we import `BioSequences` to unlock a method of `count_kmers` that works on the `LongDNA` type. In this example, we count the 1-mers of the sequence `GATTACA`. The result is a `KmerCountVector{4, 1, Int64, Vector{Int64}}`, which is a vector of 4 `Int64` elements.

```jldoctest
julia> using BioSequences

julia> count_kmers(dna"GATTACA", 1)
4-element KmerCountVector{4, 1, Int64, Vector{Int64}}:
 3
 1
 1
 2
```

We can also set the element type to `UInt16` to save memory, at the cost of a smaller maximum count value.

```jldoctest
julia> count_kmers(dna"GATTACA", 1, UInt16)
4-element KmerCountVector{4, 1, UInt16, Vector{UInt16}}:
 0x0003
 0x0001
 0x0001
 0x0002
```

!!! note
    Be careful when using element types with fewer bits, such as `UInt16`. Even though you might not expect any one k-mer to occur more than 65,535 times, some vector operations such as `LinearAlgebra.dot` and `Distances.sqeuclidean` will still use `UInt16` when summing up terms, which might lead to integer overflow.

The default `count_kmers` method is a bit different. Under the hood, it uses `count_kmers!` to modify a KmerCountVector instance in-place. This is useful when you want to count k-mers in a sequence without allocating a new vector. `count_kmers` doesn't modify a vector though, so it needs to create a new instance of a given type:

```jldoctest
julia> kc1 = count_kmers(KmerCountVector{4, 1, Int64}, [2, 0, 3, 3, 0, 1, 0])
4-element KmerCountVector{4, 1, Int64, Vector{Int64}}:
 3
 1
 1
 2
```

This default method is supposed to be as generic as possible, which is why it takes the k-mers in the form of integers already, but that's not very efficient, since that whole array would have to be allocated. Ideally, k-mers would be procedurally generated in constant memory, as is the case for the `KmerCountVector(::LongDNA, ::Integer)` method.

To reiterate, the default `count_kmers` method takes a type as its first argument so that it can create a new instance of that type, which includes creating a new vector of zeros. The `count_kmers!` method on the other hand, takes an instance of `KmerCountVector` as its first argument, and modifies it in-place, which is more flexible.

Let's create a `KmerCountVector` instance with `SparseVector` as its vector type. We can do this by passing `spzeros` to the `KmerCountVector{A, k, T}` constructor (it's just `zeros` by default):

```jldoctest
julia> using SparseArrays

julia> kc2 = KmerCountVector{4, 1, Int64}(spzeros)
4-element KmerCountVector{4, 1, Int64, SparseVector{Int64, Int64}}:
 0
 0
 0
 0

julia> count_kmers!(kc2, [2, 0, 3, 3, 0, 1, 0])
4-element KmerCountVector{4, 1, Int64, SparseVector{Int64, Int64}}:
 3
 1
 1
 2

julia> kc1 == kc2
true
```
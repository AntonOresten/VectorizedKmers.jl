```@meta
CurrentModule = VectorizedKmers
DocTestSetup = quote
    using VectorizedKmers
end
```

# Vectors of k-mer counts

How can we efficiently store multiple k-mer counts of sequences? We *could* use a normal vector: `Base.Vector{<:KmerCount}`, but remember that `KmerCount` can wrap any `AbstractVector`, including rows/columns of matrices, which means that we could store the k-mer counts of multiple sequences next to each other in a matrix (all $k$-mer counts will have a size of $A^k$). This is exactly  what the `KmerCountVector` type is for. It has two subtypes: `KmerCountColumns` and `KmerCountRows`, which wrap `AbstractMatrix` types, and store the k-mer counts as columns or rows of the matrix, respectively.

Let's create a `KmerCountColumns` meant to store the 1-mer counts of three DNA sequences. The alphabet size for DNA is 4, so each KmerCount will have a size of $A^k=4^1=4$. We'll initialize it with a matrix of zeros:

```jldoctest
julia> k = 1;

julia> n = 3;

julia> kcc = KmerCountColumns{4, k}(zeros(Int, 4^k, n))
3-element KmerCountColumns{4, 1, Int64, Matrix{Int64}}:
 KmerCount{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])
 KmerCount{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])
 KmerCount{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])
```

Oof! That does not look pretty... Let's break down what's happening here:
- The matrix consists of three different 1-mer counts, stored in columns.
- `KmerCountColumns` is an `AbstractVector`, so it gets displayed in the same way that a vector normally would, with one element per row.
- Each element of the `KmerCountColumns` is a `KmerCount` wrapped around a view of a column of the underlying matrix, hence the `SubArray` type.

We can easily access each KmerCount by index:

```jldoctest
julia> kcc[1]
4-element KmerCount{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}:
 0
 0
 0
 0
```

We can see the actual underlying matrix by looking at the `counts` field:

```jldoctest
julia> kcc.counts
4×3 Matrix{Int64}:
 0  0  0
 0  0  0
 0  0  0
 0  0  0
```

That looks much better!

### Let's count some k-mers

Now that we understand how `KmerCountVector` works, we can count k-mers of sequences using the `count_kmers!` function, which takes a `KmerCount` and a sequence as input.

Here, our `KmerCountColumns` instance, which we call `kcc`, is configured to count the 1-mers of DNA. Let's count the 1-mers of the sequence `GATTACA` in the first column of `kcc`:

```jldoctest
julia> using BioSequences

julia> count_kmers!(kcc[1], dna"GATTACA");

julia> kcc.counts
4×3 Matrix{Int64}:
 3  0  0
 1  0  0
 1  0  0
 2  0  0
```

Splendid! We can see that the 1-mers `A`, `C`, `G` and `T` occur 3, 1, 1 and 2 times, respectively, in the sequence `GATTACA`.

Although easy to visualize, 1-mer-counts of short DNA sequences are not very interesting, and are also not very useful...
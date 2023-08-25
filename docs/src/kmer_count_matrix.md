```@meta
CurrentModule = VectorizedKmers
DocTestSetup = quote
    using VectorizedKmers
end
```

# Matrices of k-mer counts

How can we efficiently store multiple k-mer count vectors of sequences? We *could* use a regular vector: `Base.Vector{<:KmerVector}`, but remember that `KmerVector` can wrap any `AbstractVector`, including rows/columns of matrices, which means that we can store the k-mer counts of multiple sequences next to each other in a matrix (all $k$-mer counts will have a size of $S^k$). This can be done using the `KmerColumns` and `KmerRows` types, which wrap `AbstractMatrix` types, and store the k-mer counts as columns and rows of the matrix, respectively.

Let's create an instance of `KmerColumns`, and configure it for storing the 1-mer counts of three DNA sequences. The alphabet size for DNA is 4, so each KmerVector will have a size of $S^k=4^1=4$. We'll initialize it with a matrix of zeros:

```jldoctest
julia> k = 1;

julia> n = 3;

julia> kc = KmerColumns{4, k}(n)
4×3 KmerColumns{4, 1, Int64, Matrix{Int64}}:
 0  0  0
 0  0  0
 0  0  0
 0  0  0
```

We can create a generator of `KmerVector`s using `eachvec`, and collect it to get a vector of `KmerVector`s.

```jldoctest
julia> collect(eachvec(kc))
3-element Vector{KmerVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}}:
 KmerVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])
 KmerVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])
 KmerVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])
```

Those types sure don't look pretty. Let's break down what's happening here:
- The matrix can be thought of as three different 1-mer count vectors of length 4, stored in columns.
- Each vector in `KmerColumns` is a `KmerVector` wrapped around a view of a column of the underlying matrix, hence the `SubArray` type.
- When collected, the generator becomes a vector of `KmerVector`s.

We can also access each column by index, which gives us a `KmerVector` wrapped around a view of the corresponding column:

```jldoctest
julia> kc[1]
4-element KmerVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}:
 0
 0
 0
 0
```

### Let's count some k-mers

Now that we understand how storing k-mer counts in matrices works, we can count k-mers of sequences using a method of `count_kmers!` that takes a `KmerVector` and a sequence as input.

Here, our `KmerColumns` instance, which we call `kc`, is configured to count the 1-mers of DNA. Let's count the 1-mers of the sequence `GATTACA` in the first column of `kc`:

```jldoctest
julia> using BioSequences

julia> count_kmers!(kc, [dna"GATTACA"])
4×3 KmerColumns{4, 1, Int64, Matrix{Int64}}:
 3  0  0
 1  0  0
 1  0  0
 2  0  0

julia> count_kmers!(kc, [dna"ACGT", dna"ATAT"], offset=1)
4×3 KmerColumns{4, 1, Int64, Matrix{Int64}}:
 3  1  2
 1  1  0
 1  1  0
 2  1  2
```

Splendid!

Although easy to visualize, 1-mer-counts of short DNA sequences are not very interesting, and are also not very useful...
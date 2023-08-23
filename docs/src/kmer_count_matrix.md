```@meta
CurrentModule = VectorizedKmers
DocTestSetup = quote
    using VectorizedKmers
end
```

# Matrices of k-mer counts

How can we efficiently store multiple k-mer count vectors of sequences? We *could* use a regular vector: `Base.Vector{<:KmerCountVector}`, but remember that `KmerCountVector` can wrap any `AbstractVector`, including rows/columns of matrices, which means that we can store the k-mer counts of multiple sequences next to each other in a matrix (all $k$-mer counts will have a size of $S^k$). This can be done using the `KmerCountColumns` and `KmerCountRows` types, which wrap `AbstractMatrix` types, and store the k-mer counts as columns and rows of the matrix, respectively.

Let's create an instance of `KmerCountColumns`, and configure it for storing the 1-mer counts of three DNA sequences. The alphabet size for DNA is 4, so each KmerCountVector will have a size of $S^k=4^1=4$. We'll initialize it with a matrix of zeros:

```jldoctest
julia> k = 1;

julia> n = 3;

julia> kcc = KmerCountColumns{4, k}(n)
4×3 KmerCountColumns{4, 1, Int64, Matrix{Int64}}:
 0  0  0
 0  0  0
 0  0  0
 0  0  0
```

We can create a generator of `KmerCountVector`s using `eachvec`, and collect it to get a vector of `KmerCountVector`s.

```jldoctest
julia> collect(eachvec(kcc))
3-element Vector{KmerCountVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}}:
 KmerCountVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])
 KmerCountVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])
 KmerCountVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}([0, 0, 0, 0])
```

Those types sure don't look pretty. Let's break down what's happening here:
- The matrix can be thought of as three different 1-mer count vectors of length 4, stored in columns.
- Each vector in `KmerCountColumns` is a `KmerCountVector` wrapped around a view of a column of the underlying matrix, hence the `SubArray` type.
- When collected, the generator becomes a vector of `KmerCountVector`s.

We can also access each column by index, which gives us a `KmerCountVector` wrapped around a view of the corresponding column:

```jldoctest
julia> kcc[1]
4-element KmerCountVector{4, 1, Int64, SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}:
 0
 0
 0
 0
```

### Let's count some k-mers

Now that we understand how storing k-mer counts in matrices works, we can count k-mers of sequences using a method of `count_kmers!` that takes a `KmerCountVector` and a sequence as input.

Here, our `KmerCountColumns` instance, which we call `kcc`, is configured to count the 1-mers of DNA. Let's count the 1-mers of the sequence `GATTACA` in the first column of `kcc`:

```jldoctest
julia> using BioSequences

julia> count_kmers!(kcc, [dna"GATTACA"])
4×3 KmerCountColumns{4, 1, Int64, Matrix{Int64}}:
 3  0  0
 1  0  0
 1  0  0
 2  0  0

julia> count_kmers!(kcc, [dna"ACGT", dna"ATAT"], offset=1)
4×3 KmerCountColumns{4, 1, Int64, Matrix{Int64}}:
 3  1  2
 1  1  0
 1  1  0
 2  1  2
```

Splendid!

Although easy to visualize, 1-mer-counts of short DNA sequences are not very interesting, and are also not very useful...
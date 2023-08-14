```@meta
CurrentModule = VectorizedKmers
DocTestSetup = quote
    using VectorizedKmers
end
```

# The `KmerCount` type

The `KmerCount` type has four type parameters, but you only really need to care about the first two: `A`, the alphabet size, and `K`, the K-mer length. So, to count the 6-mers of a DNA sequence, you would use `KmerCount{4, 6}`. For each of these K-mer counts, memory for a vector of size `A^K` is allocated, unless a vector type like `SparseVector` is used. This brings us to the two other type parameters: `T`, which is the vector element type, and `V`, which is the type of the actual vector.

Let's see it in action:

```jldoctest
julia> kc = KmerCount{4, 2}(); # creates a Vector{Int} of zeros with length 4^2

julia> using BioSequences # a weak dependency that lets us count kmers of LongDNA{4} sequences

julia> count_kmers!(kc, dna"ACGT"); # AC is 0001, CG is 0110, GT is 1011

julia> @show kc; # indices offset by 1. thanks julia...
kc = [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0]

julia> count_kmers!(kc, dna"ACGT"); # count again

julia> @show kc; # the vector gets reset to zeros before counting again
kc = [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0]

julia> count_kmers!(kc, dna"ACGT", reset=false); # avoid reset with reset=false

julia> @show kc; # look! we counted 2-mers of ACGT twice
kc = [0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0]
```
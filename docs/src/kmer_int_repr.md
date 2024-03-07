# Integer representation of k-mers

This package relies on representing k-mers as integers for indexing. This page goes over how that works exactly.

## DNA sequences

For DNA, each non-ambiguous nucleotide is assigned a number between 0 and 3:

| Nucleotide | Base-4 | Base-2 |
|------------|--------|--------|
| A          | 0      | 00     |
| C          | 1      | 01     |
| G          | 2      | 10     |
| T          | 3      | 11     |

!!! note
    Any ordering works, but this is the one used by [BioSequences.jl](https://github.com/BioJulia/BioSequences.jl). It also has some nice properties, like being in alphabetical order, and that XOR-ing a base with 3 gives you its complement.

We could theoretically convert any DNA sequence to an integer, but 64-bit unsigned integers limit us to 32-mers.

Consider the DNA sequence `GATTACA`. If we convert it to an integer using the table above, we get $2033010_4 = 10001111000100_2 = 9156_{10}$, so the integer value of `GATTACA` is 9156. Since Julia uses 1-based indexing, we would add 1 to this value to get the index for the value in a vector associated with `GATTACA`.

Writing a function for this might look like the following:

```jldoctest
julia> const DNA_ENCODING_VECTOR = zeros(UInt, 127);

julia> for (i, char) in enumerate("ACGT")
           DNA_ENCODING_VECTOR[char % Int8] = i - 1
       end;

julia> function kmer_to_int(kmer::String)
           kmer_int = zero(UInt)
           for char in kmer
               kmer_int = (kmer_int << 2) | DNA_ENCODING_VECTOR[char % Int8]
           end
           return kmer_int
       end;

julia> f(k;i=0)=[i=4i|(c%Int-1-(c=='C'))&3 for c=k][end]; # or if you're into code golf:
```

!!! note
    Consider using types like `BioSequences.LongDNA{4}` instead of `String`, as strings are not optimized for performance.

This function would not be very efficient in a practical setting (even though we're using super cool bit manipulation), since we're convert each k-mer individually, instead of having some kind of sliding window. Moreover, the function takes the k-mer in the form of a `String`, which is not ideal. The function should work as intended, though. Let's test it:

```jldoctest
julia> kmer_to_int("GATTACA")
9156

julia> f("GATTACA")
9156
```

This package does not use vectors for this, but rather K-dimensional arrays, with each dimension having a size of 4, or whatever else the alphabet size might be. Linear indexing of these arrays is still possible, and they can be lazily reshaped to vectors.

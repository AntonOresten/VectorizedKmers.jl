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

Consider the DNA sequence `GATTACA`. If we convert it to an integer using the table above, we get $2033010_4 = 10001111000100_2 = 9156_{10}$, so the integer value of `GATTACA` is 9156. Since Julia uses 1-based indexing, we would add 1 to this value to get the index of `GATTACA` in the vector.

Writing a function for this might look like the following:

```jldoctest
const DNA_ENCODING_VECTOR = zeros(UInt, 127)

for (i, char) in enumerate("ACGT")
    DNA_ENCODING_VECTOR[char % Int8] = i - 1
end

function kmer_to_int(kmer::String)
    kmer_int = zero(UInt)
    for char in kmer
        kmer_int = (kmer_int << 2) | DNA_ENCODING_VECTOR[char % Int8]
    end
    return kmer_int
end

# or if you're into code golf:
f(k;i=0)=[i=4i|(c%Int-1-(c=='C'))&3 for c=k][end]
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

## Amino acid sequences

Amino acid sequences are a little more difficult to deal with since we have multiple reading frames, and since the alphabet is much larger, the vectors grow in size much quicker as $k$ increases. However, we can still represent them as integers, like we did with DNA in the form of Strings, just without the fancy-looking bit-manipulation.

BioSequences.jl has 28 amino acids in its AminoAcidAlphabet, so we can represent each amino acid as an integer between 0 and 27.

```jldoctest
julia> using BioSequences

julia> length(AminoAcidAlphabet())
28
```

The AminoAcidAlphabet consists of amino acids with 8 bits each, so we can reinterpret them as 8-bit integers.

```jldoctest
julia> reinterpret.(Int8, [AA_A, AA_M, AA_I, AA_N, AA_O])
5-element Vector{Int8}:
  0
 12
  9
  2
 20
```

Let's say we want to convert the amino acid sequence `AMINO` to an integer. As seen above, the amino acids in the sequence have values of `0`, `12`, `9`, `2`, and `20` respectively. Thus, the integer value of the k-mer should be:

$$
0 \cdot 28^4 + 12 \cdot 28^3 + 9 \cdot 28^2 + 2 \cdot 28^1 + 20 \cdot 28^0 = 270556
$$

We can write a function for this:

```jldoctest
function kmer_to_int(kmer::LongAA)
    kmer_int = zero(UInt)
    for aa in kmer
        kmer_int = kmer_int * 28 + reinterpret(UInt8, aa)
    end
    kmer_int
end
```

To test it, we can use the `aa"..."` string macro to create a `LongAA` instance:

```jldoctest
julia> kmer_to_int(aa"AMINO")
270556
```
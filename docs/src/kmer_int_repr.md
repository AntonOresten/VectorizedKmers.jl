# Integer representation of k-mers

This package relies on representing K-mers as integers for indexing, and understanding how it works is recommended (unless you're only using higher-level API stuff).

## DNA sequences

For DNA, each non-ambiguous nucleotide is assigned a number between 0 and 3:

| Nucleotide | Base-4 | Base-2 |
|------------|--------|--------|
| A          | 0      | 00     |
| C          | 1      | 01     |
| G          | 2      | 10     |
| T          | 3      | 11     |

!!! note
    Any ordering works, but this one is the one used by [BioSequences.jl](https://github.com/BioJulia/BioSequences.jl), and it also has some nice properties, like being in alphabetical order, and that XOR-ing a base with 3 gives you its complement.

We could theoretically convert any DNA sequence to an integer, but 64-bit integers limit us to 32-mers.

Consider the DNA sequence `GATTACA`. If we convert it to an integer using the table above, we get $2033010_4 = 10001111000100_2 = 9156_{10}$, so the integer value of `GATTACA` is 9156. Since Julia uses 1-based indexing, we would add 1 to this value to get the index of `GATTACA` in the vector.

If we want to write a function for this, we may do something like the following:

```jldoctest
const DNA_ENCODING_VECTOR = zeros(Int, 127)

for (i, char) in enumerate("ACGT")
    # `c % Int8` is the ASCII value of the character `c`.
    DNA_ENCODING_VECTOR[char % Int8] = i - 1
end

function kmer_to_int(kmer::String)
    kmer_int = 0
    for char in kmer
        kmer_int = (kmer_int << 2) + DNA_ENCODING_VECTOR[charh % Int8]
    end
    kmer_int
end

# or if you like playing code golf:
f(k;i=0)=[i=(c%Int-1-(c=='C'))&3|i<<2 for c=k][end]
```

!!! note
    Strings are bad! Don't use strings! Please! They're bad! Very bad! Chars have variable length when part of Strings in Julia, so indexing and taking lengths and stuff is kinda slow! Don't use strings! Again: very bad! Use something like LongDNA{4}, please!

This function would not be very efficient in a practical setting (even though we're using super cool bit manipulation stuff), and converting each k-mer individually, instead of having some kind of sliding window. Moreover, the function takes the k-mer in the form of a `String`, which is not ideal. The function should work as intended, though. Let's test it:

```jldoctest
julia> kmer_to_int("GATTACA")
9156
```

## Amino acid sequences

Amino acid sequences are a little more difficult to deal with since there are a lot more of them, and the vectors would grow in size even quicker. However, we can still represent them as integers, but we can't use bit manipulation anymore.

BioSequences.jl has 28 amino acids in its AminoAcidAlphabet, so we can represent each amino acid as an integer between 0 and 27.
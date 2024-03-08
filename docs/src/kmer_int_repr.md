# Integer representation of K-mers

This package relies on representing K-mers as integers for indexing.

For DNA, each non-ambiguous nucleotide is assigned a number between 0 and 3:

| Nucleotide | Base-4 | Base-2 |
|------------|--------|--------|
| A          | 0      | 00     |
| C          | 1      | 01     |
| G          | 2      | 10     |
| T          | 3      | 11     |

Any ordering works, but this is the one used by [BioSequences.jl](https://github.com/BioJulia/BioSequences.jl). It also has some nice properties, like being in alphabetical order, and that XOR-ing a base with 3 gives you its complement.

We could theoretically convert any DNA sequence to an integer, but 64-bit unsigned integers limit us to 32-mers.

Consider the DNA sequence `GATTACA`. If we convert it to an integer using the table above, we get $2033010_4 = 10001111000100_2 = 9156_{10}$, so the integer value of `GATTACA` is 9156. Since Julia uses 1-based indexing, we would add 1 to this value to get the index for the value in a vector associated with `GATTACA`.

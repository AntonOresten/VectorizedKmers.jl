# <img width="25%" src="./docs/src/assets/logo.png" align="right" /> VectorizedKmers

[![Latest Release](https://img.shields.io/github/release/AntonOresten/VectorizedKmers.jl.svg)](https://github.com/AntonOresten/VectorizedKmers.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://AntonOresten.github.io/VectorizedKmers.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://AntonOresten.github.io/VectorizedKmers.jl/dev/)
[![Status](https://github.com/AntonOresten/VectorizedKmers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/AntonOresten/VectorizedKmers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/AntonOresten/VectorizedKmers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/AntonOresten/VectorizedKmers.jl)

VectorizedKmers.jl is a Julia package primarily designed for fast $K$-mer counting of biological sequences. The core idea is that $K$-mers with an alphabet size of $N$ are essentially integers in base $N$, and can be used as indices in a vector of size $N^K$ to count the corresponding $K$-mers.

This data structure can be used to quickly approximate distances between sequences. Notably, the squared Euclidean distance was used to approximate edit distance in [this paper](https://doi.org/10.1093/nar/gkz657). The dot product has also proven to be a useful metric for comparing correlation between sequences.

## Examples

```julia
julia> using VectorizedKmers, BioSequences

julia> kmer_array = count_kmers(dna"AACCGGTT", 2)
KmerArray{4, 2, Int64, Matrix{Int64}} with size (4, 4)

julia> kmer_array |> values
4×4 Matrix{Int64}:
 1  0  0  0
 1  1  0  0
 0  1  1  0
 0  0  1  1

julia> kmer_array[dna"AC"]
1

julia> kmer_array[dna"CA"]
0
```

## Limitations

The main downside of counting $K$-mers this way is that the arrays grow exponentially with respect to $K$. The 31-mer array of a DNA sequence would have a length of $4^{31} = 4,611,686,018,427,387,904$, which is equivalent to four exbibytes of memory, if the values are stored with 8-bit integers — which is just not feasible, really. Not only does allocating a lot of memory take up a lot of memory, but it can also take a substantial amount of time! This method of counting $K$-mers therefore works best for lower $K$-values.

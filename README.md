# <img width="25%" src="./docs/src/assets/logo.png" align="right" /> VectorizedKmers

[![Latest Release](https://img.shields.io/github/release/anton083/VectorizedKmers.jl.svg)](https://github.com/anton083/VectorizedKmers.jl/releases/latest)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://anton083.github.io/VectorizedKmers.jl/stable/)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://anton083.github.io/VectorizedKmers.jl/dev/)
[![Status](https://github.com/anton083/VectorizedKmers.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/anton083/VectorizedKmers.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/anton083/VectorizedKmers.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/anton083/VectorizedKmers.jl)

VectorizedKmers.jl is a Julia package primarily designed for fast $k$-mer counting of biological sequences. The core idea is that $k$-mers with an alphabet size of $S$ are essentially integers in base $S$, and can be used as indices in a vector of size $S^k$ to count the corresponding $k$-mers.

The `KmerVector` type is a wrapper for `AbstractVector`, which means that these vector k-mer counts are not limited to Julia's `Base.Vector` type; other kinds of vectors can be used as well, such as `CUDA.CuVector`, `SparseArrays.SparseVector` or even matrix views.

To efficiently group k-mer counts together, `KmerRows` and `KmerColumns` stores them in a matrix as rows or columns. These types can wrap any `AbstractMatrix`, such as `Matrix` or `CuMatrix`, and accessing its elements by index returns a `KmerVector` wrapped around a view of a row or column of the original matrix.

This data structure can be used to quickly approximate distances between sequences. Most notably, the squared Euclidean distance was used to approximate edit distance in [this paper](https://doi.org/10.1093/nar/gkz657). The dot product has also proven to be a useful metric for comparing correlation between sequences. Performing matrix multiplication on a `KmerRows` and a `KmerColumns` can therefore be used to rapidly compute a correlation score between all pairs of sequences, especially on a GPU.

## Limitations
The main downside of counting $k$-mers this way is that, unless [SparseArrays.jl](https://github.com/JuliaSparse/SparseArrays.jl) is used, the vectors will grow immensely large as $k$ increases. The 31-mer count vector of a DNA sequence would have a length of $4^{31} = 4,611,686,018,427,387,904$, which is equivalent to four exbibytes of memory, if the counts are stored as 8-bit integers — which is just not feasible, really. Not only does allocating a lot of memory take up a lot of memory, but it can also take a substantial amount of time. This method of counting $k$-mers therefore works best for lower $k$-values (unless [SparseArrays.jl](https://github.com/JuliaSparse/SparseArrays.jl) is used).

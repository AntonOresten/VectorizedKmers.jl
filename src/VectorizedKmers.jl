module VectorizedKmers

export
    # AbstractKmerCountArray
    AbstractKmerCountArray,
    get_S,
    get_k,
    counts,
    zeros!,

    # Vector
    AbstractKmerCountVector,
    KmerCountVector,

    # Matrix
    AbstractKmerCountMatrix,
    KmerCountVectors,
    KmerCountColumns,
    KmerCountRows,
    eachvec,

    # counting
    count_kmers!,
    count_kmers

include("array.jl")
include("vector.jl")
include("vectors.jl")
include("conversion.jl")
include("counting.jl")

end
module VectorizedKmers

export
    # single.jl
    AbstractKmerCount,
    KmerCount,

    # vector.jl
    AbstractKmerCountVector,
    KmerCountColumns,
    KmerCountRows,

    count_kmers!

include("KmerCount.jl")
include("KmerCountVector.jl")

end
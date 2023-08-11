module VectorizedKmers

    export AbstractKmerCount
    export KmerCount

    export AbstractKmerCountVector
    export KmerCountColumns, KmerCountRows

    export count_kmers!

include("KmerCount.jl")
include("KmerCountVector.jl")

end
module VectorizedKmers

    export KmerCountVector
    export get_A, get_k, reset!

    export AbstractKmerCountMatrix
    export KmerCountColumns, KmerCountRows

    export count_kmers!, count_kmers

include("kmer_count_vector.jl")
include("kmer_count_matrix.jl")

end
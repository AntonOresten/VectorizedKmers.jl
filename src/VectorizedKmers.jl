module VectorizedKmers

    export KmerCount
    export get_A, get_k, reset!

    export AbstractKmerCountVector
    export KmerCountColumns, KmerCountRows

    export count_kmers!, count_kmers

include("kmer_count.jl")
include("kmer_count_vec.jl")

end
module VectorizedKmers

    export AbstractKmerCount
    export KmerCount

    export AbstractKmerCountVector
    export KmerCountColumns, KmerCountRows

    export count_kmers!, count_kmers

include("kmer_count.jl")
include("kmer_count_vec.jl")

end
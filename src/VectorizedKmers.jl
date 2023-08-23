module VectorizedKmers

    export count_kmers!, count_kmers, count_kmers_gpu

    export AbstractKmerCounts
    export get_S, get_k, counts, zeros!

    export KmerCountVector

    export AbstractKmerCountMatrix
    export KmerCountColumns, KmerCountRows
    export eachvec

include("kmer_count.jl")
include("kmer_count_vector.jl")
include("kmer_count_matrix/matrix.jl")
include("count_kmers.jl")

end
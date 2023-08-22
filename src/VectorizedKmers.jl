module VectorizedKmers

    export count_kmers!, count_kmers, count_kmers_gpu

    export AbstractKmerCounts
    export get_S, get_k, counts, zeros!

    export KmerCountVector

    export AbstractKmerCountMatrix
    export KmerCountColumns, KmerCountRows

include("kmer_count.jl")
include("kmer_count_scalar.jl")
include("kmer_count_vector.jl")
include("kmer_count_matrix.jl")
include("count_kmers.jl")

end
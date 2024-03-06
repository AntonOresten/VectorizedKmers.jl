module VectorizedKmers

export KmerArray
export count_kmers!, count_kmers

include("kmer-array.jl")
include("counting.jl")

end
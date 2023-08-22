using VectorizedKmers
using Test

@testset "VectorizedKmers.jl" begin

    include("kmer_count.jl")
    include("kmer_count_vector.jl")
    include("kmer_count_matrix.jl")
    include("count_kmers.jl")

    @testset "extensions" begin
        include("ext/BioSeqCUDAExt.jl")
        include("ext/BioSequencesExt.jl")
        include("ext/CUDAExt.jl")
    end

end
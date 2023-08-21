using VectorizedKmers
using Test

@testset "VectorizedKmers.jl" begin

    include("kmer_count.jl")
    include("kmer_count_vec.jl")

    @testset "extensions" begin
        include("ext/BioSeqCUDAExt.jl")
        include("ext/BioSequencesExt.jl")
        include("ext/CUDAExt.jl")
    end

end
using VectorizedKmers
using Test

@testset "VectorizedKmers.jl" begin

    include("kmer-array.jl")
    include("counting.jl")

    @testset "extensions" begin
        include("ext/BioSequencesExt.jl")
    end

end
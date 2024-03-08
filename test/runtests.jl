using VectorizedKmers
using Test

@testset "VectorizedKmers.jl" begin

    include("KmerArray.jl")
    include("count.jl")

    @testset "extensions" begin
        include("ext/BioSequencesExt.jl")
    end

end
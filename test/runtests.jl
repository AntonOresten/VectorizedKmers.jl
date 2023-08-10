using VectorizedKmers
using Test

@testset "VectorizedKmers.jl" begin

    include("KmerCount.jl")
    include("KmerCountVector.jl")

    @testset "ext" begin
        include("ext/BioSequencesExt.jl")
        include("ext/CUDAExt.jl")
    end

end
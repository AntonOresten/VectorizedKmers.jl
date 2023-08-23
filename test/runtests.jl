using VectorizedKmers
using Test

@testset "VectorizedKmers.jl" begin

    include("array.jl")
    include("vector.jl")
    include("vectors.jl")
    include("counting.jl")

    @testset "extensions" begin
        include("ext/BioSeqCUDAExt.jl")
        include("ext/BioSequencesExt.jl")
        include("ext/CUDAExt.jl")
    end

end
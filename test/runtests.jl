using VectorizedKmers
using Test

@testset "VectorizedKmers.jl" begin

    include("kmer-array.jl")
    include("counting.jl")

    include("ext/ext.jl")

end
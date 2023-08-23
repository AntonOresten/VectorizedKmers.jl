using VectorizedKmers
using Test

@testset "VectorizedKmers.jl" begin

    include("array.jl")
    include("vector.jl")
    include("vectors.jl")
    include("conversion.jl")
    include("counting.jl")

    include("ext/ext.jl")

end
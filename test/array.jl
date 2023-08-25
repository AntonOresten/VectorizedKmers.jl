@testset "array.jl" begin

    @test hash(KmerColumns{4, 1}(16)) != hash(KmerRows{4, 2}(4))

end
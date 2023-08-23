@testset "array.jl" begin

    @test hash(KmerCountColumns{4, 1}(16)) != hash(KmerCountRows{4, 2}(4))

end
@testset "conversion.jl" begin
    
    @testset "KmerVector mapslices" begin
        kcs = KmerColumns{4, 1}([1 3; 1 1; 1 1; 1 2])
        @test KmerVector(sum, kcs) == KmerVector{4, 1}([4, 2, 2, 3])
    end

end
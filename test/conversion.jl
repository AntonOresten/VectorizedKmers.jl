@testset "conversion.jl" begin
    
    @testset "KmerCountVector mapslices" begin
        kcc = KmerCountColumns{4, 1}([1 3; 1 1; 1 1; 1 2])
        @test KmerCountVector(sum, kcc) == KmerCountVector{4, 1}([4, 2, 2, 3])
    end

end
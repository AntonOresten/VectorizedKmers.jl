@testset "counting.jl" begin

    @testset "count_kmers!" begin
        kc = KmerCountVector{4, 1}()
        @test count_kmers!(kc, [2, 0, 3, 3, 0, 1, 0]) == kc
        @test kc.counts == [3, 1, 1, 2]
    end

    @testset "count_kmers" begin
        @test count_kmers([2, 0, 3, 3, 0, 1, 0], 4, 1) == KmerCountVector{4, 1}([3, 1, 1, 2])
    end

end
@testset "counting.jl" begin

    @testset "count_kmers!" begin
        kv = KmerArray(4, 1)
        @test count_kmers!(kv, [2, 0, 3, 3, 0, 1, 0]) == kv
        @test kv == KmerArray([3, 1, 1, 2])
    end

    @testset "count_kmers" begin
        @test count_kmers([2, 0, 3, 3, 0, 1, 0], 4, 1) == KmerArray([3, 1, 1, 2])
    end

    @test_throws ErrorException VectorizedKmers.alphabet_size(String)

end
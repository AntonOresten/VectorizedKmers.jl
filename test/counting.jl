@testset "counting.jl" begin

    @testset "count_kmers!" begin
        kv = KmerVector{4, 1}()
        @test count_kmers!(kv, Ref([2, 0, 3, 3, 0, 1, 0])) == kv
        @test kv.values == [3, 1, 1, 2]
    end

    @testset "count_kmers" begin
        @test count_kmers(Ref([2, 0, 3, 3, 0, 1, 0]), 4, 1) == KmerVector{4, 1}([3, 1, 1, 2])
        @test count_kmers([Ref([2, 0, 3, 3, 0, 1, 0])], 4, 1) == KmerColumns{4, 1}(reshape([3, 1, 1, 2], :, 1))
        @test sum(count_kmers(Ref([2, 0, 3, 3, 0, 1, 0]), 4, 2).values) == 6
    end

    @test_throws ErrorException VectorizedKmers.alphabet_size(DataType)

end
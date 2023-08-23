@testset "counting.jl" begin

    @testset "count_kmers!" begin
        kc = KmerCountVector{4, 1}()
        @test count_kmers!(kc, Ref([2, 0, 3, 3, 0, 1, 0])) == kc
        @test kc.counts == [3, 1, 1, 2]
    end

    @testset "count_kmers" begin
        @test count_kmers(Ref([2, 0, 3, 3, 0, 1, 0]), 4, 1) == KmerCountVector{4, 1}([3, 1, 1, 2])
        @test count_kmers([Ref([2, 0, 3, 3, 0, 1, 0])], 4, 1) == KmerCountColumns{4, 1}(reshape([3, 1, 1, 2], :, 1))
        @test sum(count_kmers(Ref([2, 0, 3, 3, 0, 1, 0]), 4, 2).counts) == 6
    end

    @test_throws ErrorException VectorizedKmers.alphabet_size(DataType)

end
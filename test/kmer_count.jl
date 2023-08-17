@testset "kmer_count.jl" begin

    @testset "KmerCount" begin
        kc = KmerCount{4, 1}([2, 1, 1, 2])
        @test size(kc) == (4,)
        @test length(kc) == 4
        @test kc[1] == 2
        @test (kc[1] = 3) == 3
        @test kc[1] == 3
        @test eltype(kc) == Int
        @test length(KmerCount{4, 2}()) == 16

        @test get_A(kc) == 4
        @test get_k(kc) == 1

        @test summary(kc) == "KmerCount{4, 1, Int64, Vector{Int64}}"
        @test repr(kc) == "KmerCount{4, 1, Int64, Vector{Int64}}([3, 1, 1, 2])"
    end

    @testset "count_kmers!" begin
        kc = KmerCount{4, 1}()
        @test count_kmers!(kc, [2, 0, 3, 3, 0, 1, 0]) == kc
        @test kc.counts == [3, 1, 1, 2]
    end

    @testset "count_kmers" begin
        @test count_kmers(KmerCount{4, 1, Int}, [2, 0, 3, 3, 0, 1, 0]) == KmerCount{4, 1}([3, 1, 1, 2])
    end

end
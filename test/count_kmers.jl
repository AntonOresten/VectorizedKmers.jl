@testset "count_kmers.jl" begin

    @testset "count_kmers!" begin
        kc = KmerCountVector{4, 1}()
        @test count_kmers!(kc, [2, 0, 3, 3, 0, 1, 0]) == kc
        @test kc.counts == [3, 1, 1, 2]
    end

    @testset "count_kmers" begin
        @test count_kmers(KmerCountVector{4, 1, Int}, [2, 0, 3, 3, 0, 1, 0]) == KmerCountVector{4, 1}([3, 1, 1, 2])
    end

end
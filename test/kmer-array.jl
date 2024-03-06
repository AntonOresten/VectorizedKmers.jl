@testset "kmer-array.jl" begin

    @testset "KmerArray" begin
        kv = KmerArray{4, 1}([2, 1, 1, 2])
        @test size(kv) == (4,)
        @test length(kv) == 4
        @test kv[1] == 2
        kv[1] = 3
        @test kv[1] == 3
        @test repr(kv) == "KmerArray{4, 1, Int64, Vector{Int64}} with size (4,)"
    end

end
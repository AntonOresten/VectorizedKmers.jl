@testset "vector.jl" begin

    @testset "KmerVector" begin
        kv = KmerVector{4, 1}([2, 1, 1, 2])
        @test size(kv) == (4,)
        @test length(kv) == 4
        @test kv[1] == 2
        @test (kv[1] = 3) == 3
        @test kv[1] == 3
        @test eltype(kv) == Int
        @test length(KmerVector{4, 2}()) == 16

        @test get_S(kv) == 4
        @test get_k(kv) == 1

        @test summary(kv) == "KmerVector{4, 1, Int64, Vector{Int64}}"
        @test repr(kv) == "KmerVector{4, 1, Int64, Vector{Int64}}([3, 1, 1, 2])"
        @test repr(KmerVector{4, 2}()) == "KmerVector{4, 2, Int64, Vector{Int64}}([0, 0, 0, 0, 0, 0, 0, 0, 0, 0 … (6 more)])"
    end

end
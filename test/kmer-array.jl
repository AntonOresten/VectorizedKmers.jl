@testset "kmer-array.jl" begin

    @testset "KmerArray" begin
        ka = KmerArray{4, 1}([2, 1, 1, 2])
        @test ka == KmerArray([2, 1, 1, 2])
        @test size(ka) == (4,)
        @test length(ka) == 4
        @test ka[1] == 2
        ka[1] = 3
        @test ka[1] == 3
        @test repr(ka) == "KmerArray{4, 1, Int64, Vector{Int64}} with size (4,)"

        ka = KmerArray(4, 2)
        ka[1] = 3
        @test ka[1] == 3
        ka[1, 2] = 4
        @test ka[1, 2] == 4

        @test repr(MIME("text/plain"), ka) == repr(ka)
    end

end
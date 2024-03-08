@testset "KmerArray.jl" begin

    @testset "KmerArray" begin
        ka = KmerArray{4, 1}([2, 1, 1, 2])
        @test ka == KmerArray([2, 1, 1, 2])
        @test size(ka) == (4,)
        @test length(ka) == 4
        @test ka[0] == 2
        ka[0] = 3
        @test ka[0] == 3
        @test repr(ka) == "KmerArray{4, 1, Int64, Vector{Int64}} with size (4,)"

        ka = KmerArray(4, 2)
        @test values(ka) == zeros(Int, 4, 4)

        ka[0] = 3
        @test ka[0] == 3

        ka[1, 0] = 4
        @test ka[1, 0] == 4

        ka[CartesianIndex(2, 0)] = 5
        @test ka[CartesianIndex(2, 0)] == 5
        
        ka[[0, 3]] = 6
        @test ka[[0, 3]] == 6

        @test repr(MIME("text/plain"), ka) == repr(ka)
    end

end
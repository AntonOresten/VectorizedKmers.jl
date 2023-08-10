@testset "vector.jl" begin

    @testset "KmerCountVector" begin

        @testset "KmerCountColumns" begin
            @test KmerCountColumns{4, 1, Int}(zeros(Int, (4, 3))) |> length == 3
        end

        @testset "KmerCountRows" begin
            @test KmerCountRows{4, 1, Int}(zeros(Int, (3, 4))) |> length == 3
        end

    end

end
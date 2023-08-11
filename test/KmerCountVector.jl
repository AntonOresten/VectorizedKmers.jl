@testset "vector.jl" begin

    @testset "KmerCountVector" begin

        @testset "KmerCountColumns" begin
            @test KmerCountColumns{4, 2, Int}(zeros(Int, (16, 3))) |> length == 3
        end

        @testset "KmerCountRows" begin
            @test KmerCountRows{4, 2, Int}(zeros(Int, (3, 16))) |> length == 3
        end

    end

end
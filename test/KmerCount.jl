@testset "single.jl" begin

    @testset "KmerCount" begin
        @test KmerCount{4, 1, Int}() |> length == 4
    end

end
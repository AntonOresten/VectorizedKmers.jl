@testset "single.jl" begin

    @testset "KmerCount" begin
        @test KmerCount{4, 2, Int}() |> length == 16
    end

end
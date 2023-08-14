@testset "kmer_count_vec.jl" begin

    @testset "KmerCountVector" begin

        @testset "KmerCountColumns" begin
            kcc = KmerCountColumns{4, 2}(zeros(Int, (16, 3)))
            @test length(kcc) == 3
        end

        @testset "KmerCountRows" begin
            kcr = KmerCountRows{4, 2}(zeros(Int, (3, 16)))
            @test length(kcr) == 3
        end

    end

end
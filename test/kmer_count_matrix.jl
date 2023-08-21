@testset "kmer_count_vec.jl" begin

    @testset "AbstractKmerCountMatrix" begin

        @testset "KmerCountColumns" begin
            kcc = KmerCountColumns{4, 2}(zeros(Int, (16, 3)))
            @test size(kcc) == (3,)
            @test length(kcc) == 3
            @test kcc[1] == KmerCountVector{4, 2}()

            @test transpose(kcc) == kcc'

            @test all(iszero, reset!(kcc))
        end

        @testset "KmerCountRows" begin
            kcr = KmerCountRows{4, 2}(zeros(Int, (3, 16)))
            @test size(kcr) == (3,)
            @test length(kcr) == 3
            @test kcr[1] == KmerCountVector{4, 2}()

            @test transpose(kcr) == kcr'
        end

    end

end
@testset "vectors.jl" begin

    @testset "AbstractKmerCountMatrix" begin

        @testset "KmerCountColumns" begin
            kcc = KmerCountColumns{4, 2}(zeros(Int, (16, 3)))
            @test counts(kcc) == kcc.counts
            @test size(kcc) == (16, 3)
            @test kcc[1, 1] == 0
            @test length(kcc) == 3
            @test kcc[1] == KmerCountVector{4, 2}()
            @test transpose(kcc) == kcc'
            @test all(iszero, zeros!(kcc))
            kcc[1, 1] = 1
            @test kcc[1, 1] == 1
            @test eachcol(kcc)[1] == kcc[1].counts
        end

        @testset "KmerCountRows" begin
            kcr = KmerCountRows{4, 2}(zeros(Int, (3, 16)))
            @test counts(kcr) == kcr.counts
            @test size(kcr) == (3, 16)
            @test kcr[1, 1] == 0
            @test length(kcr) == 3
            @test kcr[1] == KmerCountVector{4, 2}()
            @test transpose(kcr) == kcr'
            @test all(iszero, zeros!(kcr))
            kcr[1, 1] = 1
            @test kcr[1, 1] == 1
            @test eachrow(kcr)[1] == kcr[1].counts
        end

    end

end
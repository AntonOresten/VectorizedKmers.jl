@testset "vectors.jl" begin

    @testset "AbstractKmerMatrix" begin

        @testset "KmerColumns" begin
            kcs = KmerColumns{4, 2}(zeros(Int, (16, 3)))
            @test size(kcs) == (16, 3)
            @test kcs[1, 1] == 0
            @test length(kcs) == 3
            @test kcs[1] == KmerVector{4, 2}()
            @test transpose(kcs) == kcs'
            @test all(iszero, zeros!(kcs))
            kcs[1, 1] = 1
            @test kcs[1, 1] == 1
            @test eachcol(kcs)[1] == kcs[1].values
            @test first(VectorizedKmers.eachvec(kcs)) == kcs[1]
        end

        @testset "KmerRows" begin
            krs = KmerRows{4, 2}(zeros(Int, (3, 16)))
            @test size(krs) == (3, 16)
            @test krs[1, 1] == 0
            @test length(krs) == 3
            @test krs[1] == KmerVector{4, 2}()
            @test transpose(krs) == krs'
            @test all(iszero, zeros!(krs))
            krs[1, 1] = 1
            @test krs[1, 1] == 1
            @test eachrow(krs)[1] == krs[1].values
            @test first(VectorizedKmers.eachvec(krs)) == krs[1]
        end

    end

end
using BioSequences

@testset "BioSequencesExt.jl" begin

    @testset "count_kmers!" begin

        @test count_kmers(dna"ACGT", 1).counts == [1, 1, 1, 1]
        @test count_kmers(dna"ACGT", 1) == count_kmers(LongDNA{2}(dna"ACGT"), 1)

        kc = KmerCountVector{4, 1}()
        @test count_kmers!(kc, dna"ACGT").counts == [1, 1, 1, 1]

        kc = KmerCountVector{4, 2}()
        @test count_kmers!(kc, dna"ACGT").counts == [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0]

    end
    
end
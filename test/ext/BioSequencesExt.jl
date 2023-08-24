using BioSequences

@testset "BioSequencesExt.jl" begin

    @testset "count_kmers" begin

        @testset "Single sequences" begin
            @test count_kmers(dna"ACGT", 1).counts == [1, 1, 1, 1]
            @test count_kmers(rna"ACGU", 1).counts == [1, 1, 1, 1]
            @test count_kmers(view(dna"ACGT", 1:3), 1).counts == [1, 1, 1, 0]
            @test count_kmers(view(rna"ACGU", 2:4), 1).counts == [0, 1, 1, 1]
            @test count_kmers(dna"ACGT", 1) == count_kmers(dna"ACGT", 1, UInt)
            @test count_kmers(dna"ACGT", 1) == count_kmers(LongDNA{2}(dna"ACGT"), 1)

            @test count_kmers(view(LongDNA{2}(dna"A"^32*dna"T"), 33:33), 1) == [0, 0, 0, 1]
            @test count_kmers(view(dna"A"^16*dna"T", 17:17), 1) == [0, 0, 0, 1]

            kc = KmerCountVector{4, 1}()
            @test count_kmers!(kc, dna"ACGT").counts == [1, 1, 1, 1]

            kc = KmerCountVector{4, 2}()
            @test count_kmers!(kc, dna"ACGT").counts == [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0]
        end

        @testset "Multiple sequences" begin
            seqs = [dna"ACGT", dna"GATTACA"]
            k = 1
            kcc = KmerCountColumns{4, 1}(length(seqs))
            @test count_kmers!(kcc, seqs).counts == [1 3; 1 1; 1 1; 1 2]
            @test kcc.counts == count_kmers(seqs, k).counts
        end

    end
    
end
using BioSequences

@testset "BioSequencesExt.jl" begin

    @testset "count_kmers" begin

        @testset "Single sequences" begin
            @test count_kmers(dna"ACGT", 1).values == [1, 1, 1, 1]
            @test count_kmers(rna"ACGU", 1).values == [1, 1, 1, 1]
            @test count_kmers(view(dna"ACGT", 1:3), 1).values == [1, 1, 1, 0]
            @test count_kmers(view(rna"ACGU", 2:4), 1).values == [0, 1, 1, 1]
            @test count_kmers(dna"ACGT", 1) == count_kmers(dna"ACGT", 1, UInt)
            @test count_kmers(dna"ACGT", 1) == count_kmers(LongDNA{2}(dna"ACGT"), 1)

            @test count_kmers(view(LongDNA{2}(dna"A"^32*dna"T"), 33:33), 1) == [0, 0, 0, 1]
            @test count_kmers(view(dna"A"^16*dna"T", 17:17), 1) == [0, 0, 0, 1]

            kv = KmerVector{4, 1}()
            @test count_kmers!(kv, dna"ACGT").values == [1, 1, 1, 1]

            kv = KmerVector{4, 2}()
            @test count_kmers!(kv, dna"ACGT").values == [0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0]

            seq = randdnaseq(100)
            @test all([count_kmers(seq[i:j], 5) == count_kmers(view(seq, i:j), 5) for i in 1:100, j in 1:100])

            @test count_kmers!(KmerVector{4, 5}(falses(4^5)), seq).values == count_kmers(seq, 5).values .% 2
            @test count_kmers!(KmerVector{4, 5}(falses(4^5)), view(seq, 50:100)).values == count_kmers(seq[50:100], 5).values .% 2
        end

        @testset "Multiple sequences" begin
            seqs = [dna"ACGT", dna"GATTACA"]
            k = 1
            kcs = KmerColumns{4, 1}(length(seqs))
            @test count_kmers!(kcs, seqs).values == [1 3; 1 1; 1 1; 1 2]
            @test kcs.values == count_kmers(seqs, k).values
        end

    end
    
end
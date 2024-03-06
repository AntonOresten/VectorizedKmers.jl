using BioSequences

@testset "BioSequencesExt.jl" begin

    @testset "count_kmers" begin

        @test values(count_kmers(dna"ACGT", 1)) == [1, 1, 1, 1]
        @test count_kmers(dna"ACGT", 1) == KmerArray([1, 1, 1, 1])
        @test count_kmers(rna"ACGU", 1) == KmerArray([1, 1, 1, 1])
        @test count_kmers(view(dna"ACGT", 1:3), 1) == KmerArray([1, 1, 1, 0])
        @test count_kmers(view(rna"ACGU", 2:4), 1) == KmerArray([0, 1, 1, 1])
        @test count_kmers(dna"ACGT", 1) == count_kmers(dna"ACGT", 1, UInt)
        @test count_kmers(dna"ACGT", 1) == count_kmers(LongDNA{2}(dna"ACGT"), 1)

        @test count_kmers(view(LongDNA{2}(dna"A"^32*dna"T"), 33:33), 1) == KmerArray([0, 0, 0, 1])
        @test count_kmers(view(dna"A"^16*dna"T", 17:17), 1) == KmerArray([0, 0, 0, 1])

        kv = KmerArray(4, 1)
        @test count_kmers!(kv, dna"ACGT") == KmerArray([1, 1, 1, 1])
        @test kv[LongDNA{2}(dna"A")] == 1
        @test kv[dna"A"] == 1

        kv = KmerArray(4, 2)
        @test count_kmers!(kv, dna"ACGT") == KmerArray([0; 1; 0; 0;; 0; 0; 1; 0;; 0; 0; 0; 1;; 0; 0; 0; 0])
        @test kv[LongDNA{2}(dna"AC")] == 1
        @test kv[dna"AC"] == 1

        seq = randdnaseq(50)
        @test all([count_kmers(seq[i:j], 3) == count_kmers(view(seq, i:j), 3) for i in 1:50, j in 1:50])

        @test values(count_kmers!(KmerArray{4, 5}(falses(4^5)), seq)) == values(count_kmers(seq, 5)) .% 2
        @test values(count_kmers!(KmerArray{4, 5}(falses(4^5)), view(seq, 1:50))) == values(count_kmers(seq[1:50], 5)) .% 2

        aa_seq = randaaseq(48)*aa"AY"
        @test all([count_kmers(aa_seq[i:j], 2) == count_kmers(view(aa_seq, i:j), 2) for i in 1:50, j in 1:50])
        @test count_kmers(aa_seq, 2)[aa"AY"] >= 1

    end
    
end
using BioSequences, CUDA

if CUDA.functional()

@testset "BioSeqCUDAExt.jl" begin

    # these tests only pass in the REPL for some reason
    @testset "count_kmers! LongDNA{4}" begin
        kcc = KmerCountColumns{4, 2}(CUDA.zeros(Int, 4^2, 3))

        count_kmers!(kcc, [dna"GATTACA", dna"ACGTACGT"])
        @test collect(kcc[1]) == [0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1]
        @test collect(kcc[2]) == [0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 0, 0, 0]

        count_kmers!(kcc, [dna"ACGTACGT"], offset=1, reset=false)
        @test collect(kcc[2]) == [0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 0, 0, 0] .* 2

        @test collect(kcc[3]) == zeros(Int, 16)

        count_kmers!(kcc, LongDNA{4}[], reset=true)
        @test all(iszero, kcc.counts)
    end

end

end
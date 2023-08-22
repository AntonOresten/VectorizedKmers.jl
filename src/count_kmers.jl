"""
    count_kmers!(kmer_count, sequence)

Mutates `kmer_count`.

If `reset` is `true`, the array will be zero-ed before counting.
"""
function count_kmers! end

"""
    count_kmers(sequence, k; zeros_func=zeros, reset=true)

Create a new S^k sized vector using `zeros_func` and count the k-mers in `kmers`.
"""
function count_kmers end

"""
    count_kmers_gpu(sequences, k)

Count the k-mers in `sequences` using the GPU.
"""
function count_kmers_gpu end
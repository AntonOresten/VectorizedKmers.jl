@inline Base.transpose(krs::KmerVectors{D, S, k}) where {D, S, k} = KmerVectors{3-D, S, k}(transpose(krs.values))
@inline Base.adjoint(kcm::KmerVectors) = transpose(kcm)

function KmerVector(f::Function, kvs::KmerVectors{D, S, k}) where {D, S, k}
    KmerVector{S, k}(reshape(mapslices(f, kvs.values, dims=D), :))
end
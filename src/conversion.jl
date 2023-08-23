@inline Base.transpose(kcr::KmerCountVectors{D, S, k}) where {D, S, k} = KmerCountVectors{3-D, S, k}(transpose(kcr.counts))
@inline Base.adjoint(kcm::KmerCountVectors) = transpose(kcm)

function KmerCountVector(f::Function, kcvs::KmerCountVectors{D, S, k}) where {D, S, k}
    KmerCountVector{S, k}(reshape(mapslices(f, kcvs.counts, dims=D), :))
end
@inline Base.transpose(kcr::KmerCountRows{S, k}) where {S, k} = KmerCountColumns{S, k}(transpose(kcr.counts))
@inline Base.transpose(kcc::KmerCountColumns{S, k}) where {S, k} = KmerCountRows{S, k}(transpose(kcc.counts))
@inline Base.adjoint(kcm::AbstractKmerCountMatrix) = transpose(kcm)

function KmerCountVector(f::Function, kcc::KmerCountColumns{S, k}) where {S, k}
    KmerCountVector{S, k}(reshape(mapslices(f, kcc.counts, dims=2), :))
end

function KmerCountVector(f::Function, kcr::KmerCountRows{S, k}) where {S, k}
    KmerCountVector{S, k}(reshape(mapslices(f, kcr.counts, dims=1), :))
end
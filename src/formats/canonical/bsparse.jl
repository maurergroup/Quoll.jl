
struct RealBSparseMetadata{E} <: AbstractOperatorMetadata
    sparsity::RealBlockSparsity
    basis::BasisSetMetadata{E}
    z1z2_ij2interval::Dictionary{NTuple{2, ChemicalSpecies}, Dictionary{NTuple{2, Int}, UnitRange{Int}}}
    # TODO: possibly construct this by only taking BasisSetMetadata, Atoms, Sparsity,
    # but create an outer constructor which could compute z1z2_ij2interval from those internally.
    # TODO: possibly create a constructor which could create a matrix with zero values
    # TODO: add additional field with spin metadata? This is related to considerations of
    # the ability to add arbitrary element metadata to `keydata`, how could this be extendable?
    # Would any additional metadata in `keydata` elements would have to be added to RealBSparseMetadata as a field?
    # This can be achieved by creating a function which would take relevant information such as basisf and spin from
    # RealBSparseMetadata but would also allow to add custom metadata (no need for an extra struct)
    # TODO: It would also be nice to be able to change metadata on `keydata` once RealBSparseOperator has been
    # defined already (however this wouldn't be possible without creating a new instance because RealBSparseOperator
    # is parametric with respect to KT, but that's possible without copies)
end

# TODO: for now assume each (i, j) pair maps to a 2D slice,
# but eventually this could break, e.g. in the case of gradients
# (assuming we treat gradients as a type of operatorkind instead
# of defining a separate struct).
# However, this wouldn't be too hard to fix because operatorkind
# immediately implies a particular dimension (which also suggests
# dispatching on dimension separately would not be required)
# I will usually work with the case where metadata in keydata are named tuples
# i.e. KT == Tuple{Vector{@NamedTuple{basisf::BasisMetadata{E}, spin::SpinMetadata}}, ...}
struct RealBSparseOperator{O<:AbstractOperatorKind, T<:AbstractFloat, AT, KT} <: AbstractBSparseOperator
    kind::O
    data::Dictionary{NTuple{2, ChemicalSpecies}, Array{T, 3}}
    keydata::Dictionary{NTuple{2, Int}, KeyedArray{T, 3, AT, KT}}
    metadata::RealBSparseMetadata
    # TODO: possibly rename data -> rawdata, keydata -> data
    # if we figure out that working with keydata directly results
    # in negligible overhead (even if we do probably shouldn't rename,
    # it would become inconsistent with DeepH and FHIaimsCSC formats)
end

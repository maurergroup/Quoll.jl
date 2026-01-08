### FACTORIES ###

function build_operator_data(
    ::Type{<:DenseRecipData}, metadata::M, value::T, initialised::Bool
) where {M<:AbstractMetadata,T<:Number}
    basisset = op_basisset(metadata)
    atom2nbasis = get_atom2nbasis(basisset)
    nb = sum(atom2nbasis)

    data = Array{T,2}(undef, nb, nb)
    initialised && fill!(data, value)

    return wrap_data(M, data)
end

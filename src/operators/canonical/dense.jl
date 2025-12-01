struct RealDenseMetadata{
    A<:AbstractSystem,
    S<:RealDenseSparsity,
    B<:BasisSetMetadata,
    P<:Union{SpinsMetadata, Nothing}
} <: AbstractDenseMetadata{A, S, B, P}
    atoms::A
    sparsity::S
    basisset::B
    spins::P
end

struct RecipDenseMetadata{
    A<:AbstractSystem,
    S<:RecipDenseSparsity,
    B<:BasisSetMetadata,
    P<:Union{SpinsMetadata, Nothing},
    K<:SVector{3}
} <: AbstractDenseMetadata{A, S, B, P}
    atoms::A
    sparsity::S
    basisset::B
    spins::P
    kpoint::K
end

struct DenseOperator{
    O<:OperatorKind,
    T<:Number,
    D<:AbstractArray{T},
    M<:AbstractDenseMetadata,
    KD<:AtomPairAnyDict
} <: AbstractCanonicalOperator{O, T, D, M}
    kind::O
    data::D
    metadata::M
    keydata::KD
end

function DenseOperator(kind::OperatorKind, metadata::AbstractDenseMetadata; uninit = false, value = 0.0, type = nothing)
    data = initialise_data(metadata, uninit = uninit, value = value)

    return DenseOperator(kind, data, metadata)
end

function DenseOperator(kind::OperatorKind, data::AbstractArray, metadata::AbstractDenseMetadata)
    keydata = construct_keydata(data, metadata)

    return DenseOperator(kind, data, metadata, keydata)
end

function initialise_data(metadata::RecipDenseMetadata; uninit = false, value::T₁ = 0.0, type::Type{T₂} = Nothing) where {T₁, T₂}
    basisset = get_basisset(metadata)

    uninit && @argcheck !(T₂ <: Nothing)

    if !(T₂ <: Nothing)
        T = T₂
        converted_value = convert(T₂, value)
    else
        T = T₁
        converted_value = value
    end

    atom2nbasis = get_atom2nbasis(basisset)
    nbasis = sum(atom2nbasis)
    data = Array{T, 2}(undef, nbasis, nbasis)
    uninit || fill!(data, converted_value)

    return data
end

function construct_keydata(data, metadata::RecipDenseMetadata)
    atoms = get_atoms(metadata)
    basisset = get_basisset(metadata)
    
    atom2basis = get_atom2basis(basisset)
    unique_species = unique(species(atoms, :))
    species2atom = get_species2atom(atoms)

    keydata_values = []
    ij_list = NTuple{2, Int}[]

    for z1 in unique_species, z2 in unique_species

        μ, ν = construct_axis_metadata.((z1, z2), Ref(metadata))

        for iat in species2atom[z1], jat in species2atom[z2]
            push!(ij_list, (iat, jat))
            push!(
                keydata_values,
                KeyedArray(
                    view(data, atom2basis[iat], atom2basis[jat]),
                    μ = μ,
                    ν = ν,
                )
            )
        end
    end
    keydata = Dictionary{NTuple{2, Int}, typeof(first(keydata_values))}(ij_list, keydata_values)

    return keydata
end
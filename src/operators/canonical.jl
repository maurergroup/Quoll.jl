struct CanonicalSource <: AbstractSource end

### SPHERICAL HARMONICS CONVENTION ###

#! format: off
"""
    default_shconv(source) -> SHConvention

Return the default spherical harmonics convention for the given source. Memoized.
"""
@memoize function default_shconv(source::CanonicalSource)
    return SHConvention(
        [
                        [1],
                     [1, 2, 3],
                  [1, 2, 3, 4, 5], 
               [1, 2, 3, 4, 5, 6, 7],
            [1, 2, 3, 4, 5, 6, 7, 8, 9],
        ],
        [
                        [1], 
                     [1, 1, 1],
                  [1, 1, 1, 1, 1],
               [1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1],
        ],
    )
end
#! format: on

### DATA ALIASES ###

const CanonicalBlockRealData{T} = DataContainer{
    T,3,<:SpeciesPairDictionaryArray{T,3},<:CanonicalSource,<:BlockRealSparsity
}
const CanonicalBlockRealKeyData{T} = DataContainer{
    T,3,<:AtomPairDictionary{<:KeyedArray{T,3}},<:CanonicalSource,<:BlockRealSparsity
}
const CanonicalDenseRecipKeyData{T} = DataContainer{
    T,2,<:AtomPairDictionary{<:KeyedArray{T,2}},<:CanonicalSource,<:DenseRecipSparsity
}

### METADATA ALIASES ###

# CanonicalBlockRealMetadata

const CanonicalBlockNoSpinRealMetadata{
O<:OperatorKind,
X<:CanonicalSource,
S<:BlockRealSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem
} = RealMetadata{O,X,S,B,Y,A}

const CanonicalBlockSpinRealMetadata{
O<:OperatorKind,
X<:CanonicalSource,
S<:BlockRealSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem,
P<:SpinsMetadata
} = SpinRealMetadata{O,X,S,B,Y,A,P}

const CanonicalBlockRealMetadata{
O<:OperatorKind,
X<:CanonicalSource,
S<:BlockRealSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem,
P<:SpinsMetadata
} = Union{
    <:CanonicalBlockNoSpinRealMetadata{O,X,S,B,Y,A},
    <:CanonicalBlockSpinRealMetadata{O,X,S,B,Y,A,P},
}
op_data_type(::Type{<:CanonicalBlockRealMetadata}) = CanonicalBlockRealData
op_keydata_type(::Type{<:CanonicalBlockRealMetadata}) = CanonicalBlockRealKeyData
op_source_type(::Type{<:CanonicalBlockRealMetadata}) = CanonicalSource
op_sparsity_type(::Type{<:CanonicalBlockRealMetadata}) = BlockRealSparsity

# CanonicalDenseRecipMetadata

const CanonicalDenseNoSpinRecipMetadata{
O<:OperatorKind,
X<:CanonicalSource,
S<:DenseRecipSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem,
K<:SVector{3}
} = RecipMetadata{O,X,S,B,Y,A,K}

const CanonicalDenseSpinRecipMetadata{
O<:OperatorKind,
X<:CanonicalSource,
S<:DenseRecipSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem,
P<:SpinsMetadata,
K<:SVector{3}
} = SpinRecipMetadata{O,X,S,B,Y,A,P,K}

const CanonicalDenseRecipMetadata{
O<:OperatorKind,
X<:CanonicalSource,
S<:DenseRecipSparsity,
B<:BasisSetMetadata,
Y<:SHConvention,
A<:AbstractSystem,
P<:SpinsMetadata,
K<:SVector{3}
} = Union{
    <:CanonicalDenseNoSpinRecipMetadata{O,X,S,B,Y,A,K},
    <:CanonicalDenseSpinRecipMetadata{O,X,S,B,Y,A,P,K},
}
op_data_type(::Type{<:CanonicalDenseRecipMetadata}) = DenseRecipData
op_keydata_type(::Type{<:CanonicalDenseRecipMetadata}) = CanonicalDenseRecipKeyData
op_source_type(::Type{<:CanonicalDenseRecipMetadata}) = CanonicalSource
op_sparsity_type(::Type{<:CanonicalDenseRecipMetadata}) = DenseRecipSparsity

### FACTORIES ###

function build_data(
    ::Type{<:CanonicalBlockRealData}, metadata::M, value::T, initialised::Bool
) where {M<:AbstractMetadata,T<:Number}
    basisset = op_basisset(metadata)
    species2nbasis = get_species2nbasis(basisset)
    z1z2_ij2interval = get_z1z2_ij2interval(metadata)

    data = Dictionary{NTuple{2,ChemicalSpecies},Array{T,3}}()
    for ((z1, z2), ij2interval) in pairs(z1z2_ij2interval)
        array = Array{T,3}(
            undef, species2nbasis[z1], species2nbasis[z2],
            # End interval value of the last (i, j) pair
            # (relies on the fact that ordered dictionary type is used)
            ij2interval[end][end],
        )
        initialised && fill!(array, value)
        insert!(data, (z1, z2), array)
    end

    return wrap_data(M, data)
end

function build_keydata(
    ::Type{<:CanonicalBlockRealKeyData}, metadata::M, data::CanonicalBlockRealData
) where {M<:AbstractMetadata}
    data_body = unwrap_data(data)
    sparsity = op_sparsity(metadata)
    z1z2_ij2interval = get_z1z2_ij2interval(metadata)

    keydata_values = []
    ij_contiguous = NTuple{2,Int}[]

    for (z1, z2) in keys(z1z2_ij2interval)
        μ, ν = construct_axis_metadata.(Ref(metadata), (z1, z2))

        # TODO: AxisKeys.findindex treats each element of StaticArray instance
        # as separate values to search for. Using tuples instead is a workaround
        # but it's not ideal
        for ij in keys(z1z2_ij2interval[(z1, z2)])
            push!(ij_contiguous, ij)
            push!(
                keydata_values,
                KeyedArray(
                    view(data_body[(z1, z2)], :, :, z1z2_ij2interval[(z1, z2)][ij]);
                    μ = μ,
                    ν = ν,
                    R = Tuple.(sparsity.ij2images[ij]),
                ),
            )
        end
    end

    KD_V = typeof(first(keydata_values))
    keydata = Dictionary{NTuple{2,Int},KD_V}(ij_contiguous, keydata_values)

    return wrap_data(M, keydata)
end

function build_keydata(
    ::Type{<:CanonicalDenseRecipKeyData}, metadata::M, data::DenseRecipData
) where {M<:AbstractMetadata}
    data_body = unwrap_data(data)
    atoms = op_atoms(metadata)
    basisset = op_basisset(metadata)

    atom2basis = get_atom2basis(basisset)
    unique_species = unique(species(atoms, :))
    species2atom = get_species2atom(atoms)

    keydata_values = []
    ij_list = NTuple{2,Int}[]

    for z1 in unique_species, z2 in unique_species
        μ, ν = construct_axis_metadata.(Ref(metadata), (z1, z2))

        for iat in species2atom[z1], jat in species2atom[z2]
            push!(ij_list, (iat, jat))
            push!(
                keydata_values,
                KeyedArray(
                    view(data_body, atom2basis[iat], atom2basis[jat]);
                    μ = μ,
                    ν = ν,
                ),
            )
        end
    end

    KD_V = typeof(first(keydata_values))
    keydata = Dictionary{NTuple{2,Int},KD_V}(ij_list, keydata_values)

    return wrap_data(M, keydata)
end

### SYNCHRONISATION ###

# The same method could be used for many other data types based on dictionaries
function synchronise_data!(data::CanonicalBlockRealData, comm::MPI.Comm)
    body = unwrap_data(data)
    for key in keys(body)
        MPI.Allreduce!(body[key], +, comm)
    end
end

### MISC METHODS ###

function construct_axis_metadata(
    metadata::M, z1::ChemicalSpecies
) where {M<:AbstractMetadata}
    return construct_axis_metadata(trait(SpinTrait, M), metadata, z1)
end

function construct_axis_metadata(::NoSpin, metadata::AbstractMetadata, z1::ChemicalSpecies)
    basisset = op_basisset(metadata)
    return StructArray(;
        orbital=basis_species(basisset, z1)
    )
end

function construct_axis_metadata(::HasSpin, metadata::AbstractMetadata, z1::ChemicalSpecies)
    basisset = op_basisset(metadata)
    spins = op_spins(metadata)
    return StructArray(;
        orbital=basis_species(basisset, z1),
        spin=spins_species(spins, z1),
    )
end

"""
    get_z1z2_ij2interval(atoms_or_metadata_or_operator)
    get_z1z2_ij2interval(atoms, sparsity)

Compute a species-pair-keyed dictionary mapping each `(z₁, z₂)` pair to an ordered
dictionary of `(i, j) => UnitRange` intervals. The intervals index contiguously into the
third dimension of the canonical block-real data array for that species pair. Used internally
to lay out and access per-image data slices.
"""
function get_z1z2_ij2interval(operator::AbstractOperator)
    return get_z1z2_ij2interval(op_atoms(operator), op_sparsity(operator))
end

function get_z1z2_ij2interval(metadata::AbstractMetadata)
    return get_z1z2_ij2interval(op_atoms(metadata), op_sparsity(metadata))
end

function get_z1z2_ij2interval(
    atoms::AbstractSystem, sparsity::BlockRealSparsity
)
    unique_species = unique(species(atoms, :))

    z1z2_ij2interval = Dict{
        NTuple{2,ChemicalSpecies},Dictionary{NTuple{2,Int},UnitRange{Int}}
    }()
    for z1 in unique_species, z2 in unique_species
        ij_filtered = filter(keys(sparsity.ij2images)) do ij
            i, j = ij
            species(atoms, i) == z1 && species(atoms, j) == z2
        end
        isempty(ij_filtered) && continue

        n_images = [length(sparsity.ij2images[ij]) for ij in ij_filtered]
        interval_end = cumsum(n_images)
        interval_start = interval_end .- n_images .+ 1

        # Using ordered dictionary here to ensure contiguous looping
        # in subsequent uses of ij2interval
        ij2interval = Dictionary(
            Tuple.(ij_filtered),
            (
                interval_start[i_interv]:interval_end[i_interv] for
                i_interv in axes(n_images, 1)
            ),
        )

        z1z2_ij2interval[(z1, z2)] = ij2interval
    end

    return Base.ImmutableDict(pairs(z1z2_ij2interval)...)
end

function get_z1z2_ij2offset(operator::AbstractOperator)
    return get_z1z2_ij2offset(op_atoms(operator), op_sparsity(operator))
end

function get_z1z2_ij2offset(metadata::AbstractMetadata)
    return get_z1z2_ij2offset(op_atoms(metadata), op_sparsity(metadata))
end

function get_z1z2_ij2offset(atoms::AbstractSystem, sparsity::BlockRealSparsity)
    # Converting from Base.ImmutableDict to Dictionary
    z1z2_ij2interval = dictionary(get_z1z2_ij2interval(atoms, sparsity))
    return get_z1z2_ij2offset(z1z2_ij2interval)
end

function get_z1z2_ij2offset(z1z2_ij2interval::AbstractDictionary)
    map(z1z2_ij2interval) do ij2interval
        map(ij2interval) do interval
            first(interval) - 1
        end
    end
end

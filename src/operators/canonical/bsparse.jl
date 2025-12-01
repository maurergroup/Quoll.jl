using StructArrays

# TODO: add additional field with spin metadata? This is related to considerations of
# the ability to add arbitrary element metadata to `keydata`, how could this be extendable?
# Would any additional metadata in `keydata` elements would have to be added to BSparseMetadata as a field?
# This can be achieved by creating a function which would take relevant information such as basisf and spin from
# RealBSparseMetadata but would also allow to add custom metadata (no need for an extra struct)
# It would also be nice to be able to change metadata on `keydata` once RealBSparseOperator has been
# defined already (however this wouldn't be possible without creating a new instance because BSparseOperator
# is parametric with respect to KT, but that's possible without copies)

struct RealBSparseMetadata{
    A<:AbstractSystem,
    S<:RealBlockSparsity,
    B<:BasisSetMetadata,
    P<:Union{SpinsMetadata, Nothing}
} <: AbstractBSparseMetadata{A, S, B, P}
    atoms::A
    sparsity::S
    basisset::B
    spins::P
end

struct RecipBSparseMetadata{
    A<:AbstractSystem,
    S<:RecipBlockSparsity,
    B<:BasisSetMetadata,
    P<:Union{SpinsMetadata, Nothing},
    K<:SVector{3}
} <: AbstractBSparseMetadata{A, S, B, P}
    atoms::A
    sparsity::S
    basisset::B
    spins::P
    kpoint::K
end

# TODO: Consider changing <:SpeciesPairAnyDict{<:AbstractArray{T}}
struct BSparseOperator{
    O<:OperatorKind,
    T<:Number,
    D<:SpeciesPairAnyDict{<:AbstractArray{T}},
    M<:AbstractBSparseMetadata,
    KD<:AtomPairAnyDict
} <: AbstractCanonicalOperator{O, T, D, M}
    kind::O
    data::D
    metadata::M
    keydata::KD
end

get_keydata(operator::BSparseOperator) = operator.keydata

# Constructor to initialize the operator with some data
function BSparseOperator(kind::OperatorKind, metadata::AbstractBSparseMetadata; uninit = false, value = 0.0, type = nothing)

    z1z2_ij2interval = get_z1z2_ij2interval(metadata)

    data = initialise_data(metadata, z1z2_ij2interval, uninit = uninit, value = value, type = type)
    keydata = construct_keydata(data, metadata, z1z2_ij2interval)

    return BSparseOperator(kind, data, metadata, keydata)
end

function BSparseOperator(kind::OperatorKind, data::SpeciesPairAnyDict, metadata::AbstractBSparseMetadata)

    z1z2_ij2interval = get_z1z2_ij2interval(metadata)
    keydata = construct_keydata(data, metadata, z1z2_ij2interval)

    return BSparseOperator(kind, data, metadata, keydata)
end

function initialise_data(metadata::RealBSparseMetadata, z1z2_ij2interval;
    uninit = false, value::T₁ = 0.0, type::Type{T₂} = Nothing) where {T₁, T₂}
    
    uninit && @argcheck !(T₂ <: Nothing)

    if !(T₂ <: Nothing)
        T = T₂
        converted_value = convert(T₂, value)
    else
        T = T₁
        converted_value = value
    end

    basisset = get_basisset(metadata)

    data = Dictionary{NTuple{2, ChemicalSpecies}, Array{T, 3}}()
    species2nbasis = get_species2nbasis(basisset)

    for ((z1, z2), ij2interval) in pairs(z1z2_ij2interval)
        array = Array{T, 3}(
            undef, species2nbasis[z1], species2nbasis[z2],
            # End interval value of the last (i, j) pair
            # (relies on the fact that ordered dictionary type is used)
            ij2interval[end][end] 
        )
        uninit || fill!(array, converted_value)
        insert!(data, (z1, z2), array)
    end

    return data
end

function construct_keydata(data, metadata::RealBSparseMetadata, z1z2_ij2interval)
    sparsity = get_sparsity(metadata)

    keydata_values = []
    ij_contiguous = NTuple{2, Int}[]

    for (z1, z2) in keys(z1z2_ij2interval)
        
        μ, ν = construct_axis_metadata.((z1, z2), Ref(metadata))

        # TODO: AxisKeys.findindex treats each element of StaticArray instance
        # as separate values to search for. Using tuples instead is a workaround
        # but it's not ideal
        for ij in keys(z1z2_ij2interval[(z1, z2)])
            push!(ij_contiguous, ij)
            push!(
                keydata_values,
                KeyedArray(
                    view(data[(z1, z2)], :, :, z1z2_ij2interval[(z1, z2)][ij]),
                    μ = μ,
                    ν = ν,
                    R = Tuple.(sparsity.ij2images[ij])
                )
            )
        end
    end
    keydata = Dictionary{NTuple{2, Int}, typeof(first(keydata_values))}(ij_contiguous, keydata_values)

    return keydata
end

function synchronise_data!(operator::BSparseOperator, comm::MPI.Comm)
    data = get_data(operator)
    for key in keys(data)
        MPI.Allreduce!(data[key], +, comm)
    end
end

get_z1z2_ij2interval(operator::BSparseOperator) = get_z1z2_ij2interval(get_atoms(operator), get_sparsity(operator))
get_z1z2_ij2interval(metadata::AbstractBSparseMetadata) = get_z1z2_ij2interval(get_atoms(metadata), get_sparsity(metadata))

function get_z1z2_ij2interval(atoms::AbstractSystem, sparsity::RealBlockSparsity)
    unique_species = unique(species(atoms, :))

    z1z2_ij2interval = Dict{NTuple{2, ChemicalSpecies}, Dictionary{NTuple{2, Int}, UnitRange{Int}}}()
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
            (interval_start[i_interv]:interval_end[i_interv] for i_interv in axes(n_images, 1))
        )

        z1z2_ij2interval[(z1, z2)] = ij2interval
    end

    return Base.ImmutableDict(pairs(z1z2_ij2interval)...)
end

get_z1z2_ij2offset(operator::BSparseOperator) = get_z1z2_ij2offset(get_atoms(operator), get_sparsity(operator))
get_z1z2_ij2offset(metadata::AbstractBSparseMetadata) = get_z1z2_ij2offset(get_atoms(metadata), get_sparsity(metadata))

function get_z1z2_ij2offset(atoms::AbstractSystem, sparsity::AbstractBlockSparsity)
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

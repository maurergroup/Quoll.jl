using StructArrays

# TODO: add additional field with spin metadata? This is related to considerations of
# the ability to add arbitrary element metadata to `keydata`, how could this be extendable?
# Would any additional metadata in `keydata` elements would have to be added to BSparseMetadata as a field?
# This can be achieved by creating a function which would take relevant information such as basisf and spin from
# RealBSparseMetadata but would also allow to add custom metadata (no need for an extra struct)
# It would also be nice to be able to change metadata on `keydata` once RealBSparseOperator has been
# defined already (however this wouldn't be possible without creating a new instance because BSparseOperator
# is parametric with respect to KT, but that's possible without copies)

struct RealBSparseMetadata{A<:AbstractSystem, S<:RealBlockSparsity, B<:BasisSetMetadata, P<:Union{SpinsMetadata, Nothing}} <: AbstractBSparseMetadata{A, S, B, P}
    atoms::A
    sparsity::S
    basisset::B
    spins::P
end

struct RecipBSparseMetadata{A<:AbstractSystem, S<:RecipBlockSparsity, B<:BasisSetMetadata, P<:Union{SpinsMetadata, Nothing}} <: AbstractBSparseMetadata{A, S, B, P}
    atoms::A
    sparsity::S
    basisset::B
    spins::P
    kpoint::SVector{3, Float64}
end

# I will usually work with the case where metadata in keydata are named tuples
# i.e. KT == Tuple{Vector{@NamedTuple{basisf::BasisMetadata{E}, spin::SpinMetadata}}, ...}
# TODO: Consider changing <:SpeciesPairAnyDict{<:AbstractArray{T}}
struct BSparseOperator{O<:OperatorKind, T<:Number, D<:SpeciesPairAnyDict{<:AbstractArray{T}}, M<:AbstractBSparseMetadata, KD<:AtomPairAnyDict} <: AbstractBSparseOperator{O, T, D, M, KD}
    kind::O
    data::D
    metadata::M
    keydata::KD
end

get_keydata(operator::BSparseOperator) = operator.keydata

# Constructor to initialize the operator with zero values
# TODO: possibly split into smaller functions? Also this function could accept
# additional metadata to go into keydata
function BSparseOperator(kind::OperatorKind, metadata::AbstractBSparseMetadata; float::Type{T} = Float64) where T

    z1z2_ij2interval = get_z1z2_ij2interval(metadata)
    basisset = get_basisset(metadata)
    spins = get_spins(metadata)
    sparsity = get_sparsity(metadata)

    # Construct data
    data = Dictionary{NTuple{2, ChemicalSpecies}, Array{T, 3}}()
    species2nbasis = get_species2nbasis(basisset)

    for ((z1, z2), ij2interval) in pairs(z1z2_ij2interval)
        insert!(
            data, (z1, z2),
            zeros(
                T, species2nbasis[z1], species2nbasis[z2],
                # End interval value of the last (i, j) pair
                # (relies on the fact that ordered dictionary type is used)
                ij2interval[end][end] 
            )
        )
    end

    # Construct keydata
    keydata_values = []
    ij_contiguous = NTuple{2, Int}[]

    for (z1, z2) in keys(z1z2_ij2interval)
        
        if metadata.spins !== nothing
            μ = StructArray(
                orbital = basis_species(basisset, z1),
                spin = spins_species(spins, z1),
            )
            ν = StructArray(
                orbital = basis_species(basisset, z2),
                spin = spins_species(spins, z2),
            )
        else
            μ = StructArray(
                orbital = basis_species(basisset, z1),
            )
            ν = StructArray(
                orbital = basis_species(basisset, z2),
            )
        end

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

    return BSparseOperator(kind, data, metadata, keydata)
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

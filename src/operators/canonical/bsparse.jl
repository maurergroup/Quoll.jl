using StructArrays

struct RealBSparseMetadata{A<:AbstractSystem, E} <: AbstractOperatorMetadata
    atoms::A
    sparsity::RealBlockSparsity
    basisset::BasisSetMetadata{E}
    spins::Union{SpinsMetadata, Nothing}
    z1z2_ij2interval::Dictionary{NTuple{2, ChemicalSpecies}, Dictionary{NTuple{2, Int}, UnitRange{Int}}}
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
# Right now I assume N = 3 in all methods, but I could modify this
# type by making it parametric on N and then modifying the methods
# which I could modify based on OperatorKind instead of N
# because a particular OperatorKind immediately implies particular N
# 
# I will usually work with the case where metadata in keydata are named tuples
# i.e. KT == Tuple{Vector{@NamedTuple{basisf::BasisMetadata{E}, spin::SpinMetadata}}, ...}
struct RealBSparseOperator{O<:OperatorKind, T<:AbstractFloat, A<:AbstractSystem, E, AT, KT} <: AbstractBSparseOperator
    kind::O
    data::Dictionary{NTuple{2, ChemicalSpecies}, Array{T, 3}}
    keydata::Dictionary{NTuple{2, Int}, KeyedArray{T, 3, AT, KT}}
    metadata::RealBSparseMetadata{A, E}
end

get_keydata(operator::RealBSparseOperator) = operator.keydata
get_z1z2_ij2interval(operator::RealBSparseOperator) = operator.metadata.z1z2_ij2interval
get_z1z2_ij2interval(metadata::RealBSparseMetadata) = metadata.z1z2_ij2interval

# Constructor to initialize the operator with zero values
# TODO: possibly split into smaller functions?
function RealBSparseOperator(kind::OperatorKind, metadata::RealBSparseMetadata; T = Float64)

    z1z2_ij2interval = get_z1z2_ij2interval(metadata)
    basisset = get_basisset(metadata)
    spins = get_spins(metadata)
    sparsity = get_sparsity(metadata)

    # Construct data
    data = Dictionary{NTuple{2, ChemicalSpecies}, Array{T, 3}}()
    species2nbasis = get_species2nbasis(basisset)

    for ((z1, z2), interval) in pairs(z1z2_ij2interval)
        insert!(
            data, (z1, z2),
            zeros(
                T, species2nbasis[z1], species2nbasis[z2],
                # End interval value of the last (i, j) pair
                interval[end][end] 
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

    return RealBSparseOperator(kind, data, keydata, metadata)
end

function get_z1z2_ij2interval(atoms::AbstractSystem, sparsity::RealBlockSparsity)
    unique_species = unique(species(atoms, :))

    z1z2_ij2interval = Dictionary{NTuple{2, ChemicalSpecies}, Dictionary{NTuple{2, Int}, UnitRange{Int}}}()
    for z1 in unique_species, z2 in unique_species
        ij_filtered = filter(
            ij -> (species(atoms, ij[1]), species(atoms, ij[2])) == (z1, z2),
            keys(sparsity.ij2images)
        )
        isempty(ij_filtered) && continue

        n_images = [length(sparsity.ij2images[ij]) for ij in ij_filtered]
        interval_end = cumsum(n_images)
        interval_start = interval_end .- n_images .+ 1

        # Using ordered dictionary here to ensure contiguous looping
        # in subsequent uses of ij2interval
        ij2interval = Dictionary(
            ij_filtered,
            (interval_start[i_interv]:interval_end[i_interv] for i_interv in axes(n_images, 1))
        )
        insert!(z1z2_ij2interval, (z1, z2), ij2interval)
    end

    return z1z2_ij2interval
end

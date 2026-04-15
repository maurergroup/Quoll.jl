"""
    BasisMetadata{E}

Metadata for a single basis function: chemical species `z`, principal quantum number `n`,
angular momentum `l`, magnetic quantum number `m` (always in wiki SH convention), and
optional `extras` (e.g. degeneracy labels). Two `BasisMetadata` are equal if all fields
match; they [`share_subblock`](@ref) if everything except `m` matches.
"""
struct BasisMetadata{E}
    z::ChemicalSpecies
    n::Int
    l::Int
    m::Int
    extras::E
end

function Base.hash(b::BasisMetadata, h::UInt)
    return hash(
        getfield(b, :z),
        hash(
            getfield(b, :n),
            hash(
                getfield(b, :l),
                hash(
                    getfield(b, :m),
                    hash(
                        getfield(b, :extras),
                        hash(:BasisMetadata, h),
                    ),
                ),
            ),
        ),
    )
end

function Base.:(==)(b1::BasisMetadata, b2::BasisMetadata)
    return (
        getfield(b1, :z) == getfield(b2, :z) &&
        getfield(b1, :n) == getfield(b2, :n) &&
        getfield(b1, :l) == getfield(b2, :l) &&
        getfield(b1, :m) == getfield(b2, :m) &&
        getfield(b1, :extras) == getfield(b2, :extras)
    )
end

function Base.isequal(b1::BasisMetadata, b2::BasisMetadata)
    return (
        isequal(getfield(b1, :z), getfield(b2, :z)) &&
        isequal(getfield(b1, :n), getfield(b2, :n)) &&
        isequal(getfield(b1, :l), getfield(b2, :l)) &&
        isequal(getfield(b1, :m), getfield(b2, :m)) &&
        isequal(getfield(b1, :extras), getfield(b2, :extras))
    )
end

"""
    share_subblock(b1::BasisMetadata, b2::BasisMetadata) -> Bool

Return `true` if `b1` and `b2` belong to the same angular momentum subblock (same species,
`n`, `l`, and extras — only `m` may differ).
"""
function share_subblock(b1::BasisMetadata, b2::BasisMetadata)
    return (
        isequal(getfield(b1, :z), getfield(b2, :z)) &&
        isequal(getfield(b1, :n), getfield(b2, :n)) &&
        isequal(getfield(b1, :l), getfield(b2, :l)) &&
        isequal(getfield(b1, :extras), getfield(b2, :extras))
    )
end

function Base.show(io::IO, basisf::BasisMetadata)
    return print(io, "$(basisf.z)$(basisf.n)$(basisf.l)$(basisf.m)")
end

"""
    BasisSetMetadata{B,A}

A basis set. Contains a `basis` dictionary mapping each `ChemicalSpecies` to its vector
of [`BasisMetadata`](@ref) entries, and an `atom2species` vector mapping atom indices
to species.
"""
struct BasisSetMetadata{B<:SpeciesAnyDict,A<:AbstractVector{ChemicalSpecies}}
    basis::B
    atom2species::A
end

function basis_atom(basisset::BasisSetMetadata, iat::Int)
    return basisset.basis[basisset.atom2species[iat]]
end

function basis_species(basisset::BasisSetMetadata, species::ChemicalSpecies)
    return basisset.basis[species]
end

get_unique_species(basisset::BasisSetMetadata) = unique(basisset.atom2species)

"""
    get_angular_momenta(basis) -> Vector{Int}

Return the angular momentum `l` of each subblock's representative basis function.
"""
function get_angular_momenta(basis::AbstractVector{T}) where {T<:BasisMetadata}
    return getfield.(get_basis_per_subblock(basis), :l)
end

"""
    get_subblock_ranges(basisset::BasisSetMetadata)
    get_subblock_ranges(basis_z::AbstractVector{<:BasisMetadata})

Return contiguous index ranges for each angular momentum subblock. The `BasisSetMetadata`
form returns a `Dictionary{ChemicalSpecies, Vector{UnitRange{Int}}}`. It is assumed that
all basis functions belonging to a sublock are next to each other.
"""
function get_subblock_ranges(basisset::BasisSetMetadata)
    ranges_dict = Dictionary{ChemicalSpecies,Vector{UnitRange{Int}}}()
    for z in keys(basisset.basis)
        insert!(ranges_dict, z, get_subblock_ranges(basis_species(basisset, z)))
    end
    return ranges_dict
end

function get_subblock_ranges(basis_z::AbstractVector{<:BasisMetadata})
    # Get a list which shows whether the next basis function belongs to the same subblock.
    # If false, then the next basis function belongs to a different subblock,
    # which makes it straightforward to get ranges of subblocks. An additional zero
    # is pushed to the end of the list to make the subsequent loop easier
    next_in_subblock_list = map(
        i -> share_subblock(basis_z[i], basis_z[i + 1]),
        Base.Iterators.take(eachindex(basis_z), length(basis_z) - 1),
    )
    push!(next_in_subblock_list, false)

    ranges = UnitRange{Int}[]
    interval_start = firstindex(basis_z)

    for (i, next_in_subblock) in zip(eachindex(basis_z), next_in_subblock_list)
        if !next_in_subblock
            interval_end = i
            push!(ranges, interval_start:interval_end)

            # Update the start of the interval
            interval_start = interval_end + 1
        end
    end

    return ranges
end

"""
    get_indices_in_subblock(basisset::BasisSetMetadata)
    get_indices_in_subblock(basis_z::AbstractVector{<:BasisMetadata})

Return the 1-based position of each basis function within its subblock. Unlike `m`, this
reflects the actual ordering after SH convention reordering.
"""
function get_indices_in_subblock(basisset::BasisSetMetadata)
    indices_dict = Dictionary{ChemicalSpecies,Vector{Int}}()
    for z in keys(basisset.basis)
        insert!(indices_dict, z, get_indices_in_subblock(basis_species(basisset, z)))
    end
    return indices_dict
end

function get_indices_in_subblock(basis_z::AbstractVector{<:BasisMetadata})
    ranges = get_subblock_ranges(basis_z)
    norm_ranges = [first(ranges)...]

    for i in Base.Iterators.drop(eachindex(ranges), 1)
        cumsum_prev = last(ranges[i - 1])
        range = (first(ranges[i]) - cumsum_prev):(last(ranges[i]) - cumsum_prev)
        push!(norm_ranges, range...)
    end
    return norm_ranges
end

"""
    convert_basisset_shconv(basisset, shconv) -> BasisSetMetadata

Reorder basis functions according to the spherical harmonics convention `shconv`. The input
basis is assumed to be in the wiki (standard) convention.
"""
function convert_basisset_shconv(basisset::BasisSetMetadata, shconv::SHConvention)
    in_basis = basisset.basis
    out_basis = convert_speciesdict_shconv(in_basis, in_basis, shconv)
    return BasisSetMetadata(out_basis, basisset.atom2species)
end

"""
    lift_arbitrary_degeneracy(basis)
    lift_arbitrary_degeneracy(basis_z)

When multiple subblocks share the same `(z, n, l)` (arbitrary degeneracy), add a `:dgen`
label to the `extras` field of each basis function to distinguish them.
"""
function lift_arbitrary_degeneracy(
    basis::Base.ImmutableDict{ChemicalSpecies,B}
) where {K,V,E<:Base.ImmutableDict{K,V},T<:BasisMetadata{E},B<:Vector{T}}
    # Check if there are any arbitrary degeneracies
    if !is_arbitrary_degeneracy(basis)
        return basis
    end

    # Construct new key-value types because a new key-value pair will be added
    Kₒᵤₜ = typejoin(K, Symbol)
    Vₒᵤₜ = typejoin(V, Symbol)

    newbasis = Dict{ChemicalSpecies,Vector{BasisMetadata{Base.ImmutableDict{Kₒᵤₜ,Vₒᵤₜ}}}}()
    for z in keys(basis)
        newbasis[z] = lift_arbitrary_degeneracy(basis[z])
    end

    return Base.ImmutableDict(pairs(newbasis)...)
end

# Split the arbitrarily degenerate basis functions into separate subsets
# - Iteratively search the first m=(-l:l) basis functions
# - Works assuming that the basis functions are in separate subblocks
#   and that each subblock is full (all ms for a given l present)
# Add additional labels to those subsets to extras (`:dgen = 1`) 
function lift_arbitrary_degeneracy(
    basis_z::AbstractVector{T}
) where {K,V,E<:Base.ImmutableDict{K,V},T<:BasisMetadata{E}}

    # Construct new key-value types because a new key-value pair will be added
    Kₒᵤₜ = typejoin(K, Symbol)
    Vₒᵤₜ = typejoin(V, Symbol)

    basis_extras_pairs =
        convert.(
            Vector{Pair{Kₒᵤₜ,Vₒᵤₜ}},
            collect.(pairs.(getfield.(basis_z, :extras))),
        )

    basis_per_subblock = get_basis_per_subblock(basis_z)
    for b in basis_per_subblock
        shared_ib_list = findall(b2 -> share_subblock(b, b2), basis_z)
        nb_subblock = length(shared_ib_list)
        nb_froml = (2b.l + 1)

        # If arbitrary degeneracy is present
        if nb_subblock > nb_froml
            n_subblocks = div(nb_subblock, nb_froml)
            for i_subblock in 1:n_subblocks

                # Find subblock basis indices in an arbitrarily degenerate larger subblock
                subblock_ib_list = shared_ib_list[[
                    findfirst(b2 -> isequal(m, b2.m), basis_z[shared_ib_list])
                    for m in ((-b.l):(b.l))
                ]]

                for ib in subblock_ib_list
                    # Add additional field in extras to the found basis functions
                    push!(basis_extras_pairs[ib], :dgen => Symbol(i_subblock))
                    # Remove the found basis functions from shared_ib_list
                    popat!(shared_ib_list, findfirst(isequal(ib), shared_ib_list))
                end
            end
        end
    end

    # Construct basis functions with the new extra metadata
    newbasis_z = [
        BasisMetadata(
            b.z, b.n, b.l, b.m,
            Base.ImmutableDict(p...),
        )
        for (b, p) in zip(basis_z, basis_extras_pairs)
    ]
    return newbasis_z
end

"""
    is_arbitrary_degeneracy(basisset) -> Bool
    is_arbitrary_degeneracy(basis_z) -> Bool

Check whether any species has multiple subblocks sharing the same `(n, l)` quantum numbers.
"""
is_arbitrary_degeneracy(basisset::BasisSetMetadata) = is_arbitrary_degeneracy(basisset.basis)

function is_arbitrary_degeneracy(basis::SpeciesAnyDict)
    return any(is_arbitrary_degeneracy.(values(basis)))
end

# Check if there are any duplicates in the basis
# - Filter all basis functions to only a single basis function per subblock
# - Check how many basis functions are equal to that basis function
# - Check if that number agrees with the angular momentum of that basis function
#   (check if not more, can be fewer if not all basis functions are used)
function is_arbitrary_degeneracy(basis_z::AbstractVector{<:BasisMetadata})
    basis_per_subblock = get_basis_per_subblock(basis_z)
    for b in basis_per_subblock
        ib_all = findall(b2 -> share_subblock(b, b2), basis_z)
        nb_subblock = length(ib_all)
        nb_froml = (2b.l + 1)
        if nb_subblock > nb_froml
            return true
        end
    end
    return false
end

function get_basis_per_subblock(basis_z::AbstractVector{<:BasisMetadata})
    basis_per_subblock_repeated = map(basis_z) do b
        ib_first = findfirst(b2 -> share_subblock(b, b2), basis_z)
        basis_z[ib_first]
    end
    return unique(basis_per_subblock_repeated)
end

"""
    reduce_basisset(basisset, subbasis; inverted=false) -> BasisSetMetadata
    reduce_basisset(basisset, subbasis_masks) -> BasisSetMetadata

Return a new `BasisSetMetadata` containing only the basis functions selected by `subbasis`
(a vector of `BasisMetadata` to keep). If `inverted=true`, keep the complement instead.
The second form accepts precomputed per-species `BitVector` masks.
"""
function reduce_basisset(
    basisset::BasisSetMetadata, subbasis::AbstractVector{<:BasisMetadata};
    inverted=false,
)
    subbasis_masks = get_subbasis_masks(basisset, subbasis; inverted=inverted)
    return reduce_basisset(basisset, subbasis_masks)
end

function reduce_basisset(
    basisset::BasisSetMetadata, subbasis_masks::SpeciesDictionary
)
    unique_species_keys = Indices(keys(basisset.basis))
    subbasis_dict = map(unique_species_keys) do z
        basis_z = basis_species(basisset, z)
        basis_z[subbasis_masks[z]]
    end
    return BasisSetMetadata(
        Base.ImmutableDict(pairs(subbasis_dict)...), basisset.atom2species
    )
end

"""
    get_subbasis_masks(basisset, subbasis; inverted=false) -> Dictionary{ChemicalSpecies, BitVector}

Compute per-species boolean masks indicating which basis functions belong to `subbasis`.
If `inverted=true`, the masks select the complement (i.e. everything *not* in `subbasis`).
"""
function get_subbasis_masks(
    basisset::BasisSetMetadata, subbasis::AbstractVector{<:BasisMetadata};
    inverted=false,
)
    unique_species_keys = Indices(keys(basisset.basis))

    return map(unique_species_keys) do z
        basis_z = basis_species(basisset, z)
        subbasis_z = filter(orbital -> isequal(orbital.z, z), subbasis)

        # Throw a warning if some orbitals in subbasis are not present in basisset
        # (which could happen if the orbitals were supplied incorrectly in the input file)
        if !all(subbasis_z .∈ Ref(basis_z))
            @warn "Some orbitals in the subbasis are not present in the full basis"
        end

        inverted ? basis_z .∉ Ref(subbasis_z) : basis_z .∈ Ref(subbasis_z)
    end
end

"""
    get_dense_subbasis_mask(basisset, subbasis; inverted=false) -> BitVector
    get_dense_subbasis_mask(basisset, subbasis_masks) -> BitVector

Flatten per-species subbasis masks into a single `BitVector` over all atoms (following
`atom2species` ordering), suitable for indexing into dense matrices.
"""
function get_dense_subbasis_mask(
    basisset::BasisSetMetadata, subbasis::AbstractVector{<:BasisMetadata};
    inverted=false,
)
    subbasis_masks = get_subbasis_masks(basisset, subbasis; inverted=inverted)
    return get_dense_subbasis_mask(basisset, subbasis_masks)
end

function get_dense_subbasis_mask(
    basisset::BasisSetMetadata, subbasis_masks::Dictionary{ChemicalSpecies,BitVector}
)
    return reduce(vcat, (subbasis_masks[z] for z in basisset.atom2species))
end

"""
    get_atom2nbasis(basisset) -> Vector{Int}

Return the number of basis functions per atom.
"""
function get_atom2nbasis(basisset::BasisSetMetadata)
    return get_atom2nbasis(basisset.basis, basisset.atom2species)
end
get_atom2nbasis(basis, atom2species) = [length(basis[z]) for z in atom2species]

"""
    get_atom2basis(basisset) -> Vector{UnitRange{Int}}

Return the global basis function index range for each atom.
"""
function get_atom2basis(basisset::BasisSetMetadata)
    return get_atom2basis(basisset.basis, basisset.atom2species)
end

function get_atom2basis(basis, atom2species)
    atom2nbasis = get_atom2nbasis(basis, atom2species)
    return get_atom2basis(atom2nbasis)
end

function get_atom2basis(atom2nbasis)
    n_atoms = length(atom2nbasis)
    atom2basis_start = cumsum(atom2nbasis) .- atom2nbasis .+ 1
    return [
        atom2basis_start[iat]:(atom2basis_start[iat] + atom2nbasis[iat] - 1)
        for iat in 1:n_atoms
    ]
end

"""
    get_basis2atom(basisset) -> Vector{Int}

Return the atom index for each global basis function index.
"""
function get_basis2atom(basisset::BasisSetMetadata)
    return get_basis2atom(basisset.basis, basisset.atom2species)
end

function get_basis2atom(basis, atom2species)
    atom2nbasis = get_atom2nbasis(basis, atom2species)
    return get_basis2atom(atom2nbasis)
end

function get_basis2atom(atom2nbasis)
    return reduce(vcat, [fill(iat, nbasis) for (iat, nbasis) in enumerate(atom2nbasis)])
end

"""
    get_species2nbasis(basisset) -> ImmutableDict{ChemicalSpecies, Int}

Return the number of basis functions per chemical species.
"""
function get_species2nbasis(basisset::BasisSetMetadata)
    d = Base.ImmutableDict(
        (
            z => length(basis_species(basisset, z))
            for z in keys(basisset.basis)
        )...,
    )
    return d
end

"""
    get_atom2offset(basisset) -> Vector{Int}

Return the cumulative basis function offset for each atom (number of basis functions
belonging to atoms with lower index).
"""
function get_atom2offset(basisset::BasisSetMetadata)
    atom2basis = get_atom2basis(basisset)
    return get_atom2offset(atom2basis)
end

get_atom2offset(atom2basis) = [first(interval) - 1 for interval in atom2basis]

function precompute_shifts(basisset::BasisSetMetadata, shconv::SHConvention)
    return precompute_shifts(basisset.basis, shconv)
end

function precompute_orders(basisset::BasisSetMetadata, shconv::SHConvention)
    return precompute_orders(basisset.basis, shconv)
end

function precompute_shphases(
    basisset::BasisSetMetadata, shconv::SHConvention, ::Val{D}=Val(1), ::Type{T}=Float64
) where {D,T}
    return precompute_shphases(basisset.basis, shconv, Val(D), T)
end

function precompute_signed_perm_matrices(
    basisset::BasisSetMetadata, shconv::SHConvention, ::Type{T}=Float64
) where {T}
    return precompute_signed_perm_matrices(basisset.basis, shconv, T)
end
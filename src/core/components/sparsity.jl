"""
    AbstractSparsity

Abstract supertype for all sparsity patterns. All concrete subtypes store a `hermitian`
flag (accessible via `op_hermicity`) and, for real-space types, an `images` vector of
lattice translation vectors.
"""
abstract type AbstractSparsity end
op_hermicity(sparsity::AbstractSparsity) = sparsity.hermitian
op_images(sparsity::AbstractSparsity) = sparsity.images

struct CSCRealSparsity <: AbstractSparsity
    rowval::Vector{Int}
    colcellptr::Array{Int,3}
    images::Vector{SVector{3,Int}}
    hermitian::Bool
end

struct CSCRecipSparsity <: AbstractSparsity
    rowval::Vector{Int}
    colptr::Vector{Int}
    hermitian::Bool
end

struct BlockRealSparsity <: AbstractSparsity
    ij2images::Dictionary{NTuple{2,Int},Vector{SVector{3,Int}}}
    images::Vector{SVector{3,Int}}
    hermitian::Bool
end

struct BlockRecipSparsity <: AbstractSparsity
    ij::Vector{NTuple{2,Int}}
    hermitian::Bool
end

struct DenseRealSparsity <: AbstractSparsity
    images::Vector{SVector{3,Int}}
    hermitian::Bool
end

struct DenseRecipSparsity <: AbstractSparsity
    hermitian::Bool
end

### SPARSITY FROM NEIGHBOURLIST ###

"""
    build_sparsity(::Type{S}, atoms, radii; hermitian=false)

Build a sparsity pattern of type `S` from an atomic system and per-species interaction
radii. Constructs a neighbour list and retains pairs within the sum of their radii.
"""
function build_sparsity(
    ::Type{S}, atoms::AbstractSystem, radii::SpeciesAnyDict;
    hermitian=false,
) where {S<:AbstractSparsity}
    @debug "Computing sparsity manually from radii"
    maxedges = get_maxedges(radii)
    nlist_cut = maximum(values(maxedges))
    nlist = PairList(atoms, nlist_cut)
    return build_sparsity(S, atoms, nlist, maxedges; hermitian=hermitian)
end

function build_sparsity(
    ::Type{<:BlockRealSparsity},
    atoms::AbstractSystem,
    nlist::PairList,
    maxedges::SpeciesPairAnyDict;
    hermitian=false,
)
    unique_ij = unique(zip(nlist.i, nlist.j))
    if hermitian
        unique_ij = [(i, j) for (i, j) in unique_ij if j ≥ i]
    end

    ij2images = Dictionary(unique_ij, [SVector{3,Int}[] for ij in unique_ij])

    maxedges_nounits = Base.ImmutableDict(
        (
            key => ustrip(value)
            for (key, value) in pairs(maxedges)
        )...,
    )

    for (pair, image) in zip(pairs(nlist), nlist.S)
        i, j = pair[1], pair[2]
        hermitian && i > j && continue

        # Sift edges which are beyond the cutoff value
        r = norm(pair[3])
        if r ≤ maxedges_nounits[(species(atoms, i), species(atoms, j))]
            push!(ij2images[(i, j)], image)
        end
    end

    # Add self-interactions for R = [0, 0, 0]
    # This results in non-hermitian on-site ij2images even when hermitian is true.
    # Before having seen element-wise sparsity, there is no way to tell whether
    # any on-site blocks will only contain data in the upper diagonal.
    # Therefore, to be consistent with the rest of the code, I'm assuming that
    # on-site ij2images always contain mirror images even when they are
    # not strictly required e.g. during sparsity conversion from CSC to block-wise
    for (iat, _) in enumerate(atoms)
        ii = (iat, iat)

        # Add empty array if the key doesn't exist yet
        # (would be in the case of no i self-interactions across periodic boundaries)
        get!(Vector{SVector{3,Int}}, ij2images, ii)
        push!(ij2images[ii], SA[0, 0, 0])
    end

    return BlockRealSparsity(ij2images, unique(nlist.S), hermitian)
end

function get_maxedges(radii::SpeciesAnyDict{T}) where {T<:Number}
    angstrom = unit(T) == NoUnits ? u"Å" : one(T)
    maxedges = Dict(
        (zi, zj) => (i_cut + j_cut) * angstrom
        for (zi, i_cut) in pairs(radii)
        for (zj, j_cut) in pairs(radii)
    )
    return maxedges
end

### SPARSITY CONVERSIONS ###

"""
    convert_sparsity(::Type{Sₒᵤₜ}, in_sparsity, basisset; hermitian=false)

Convert a sparsity pattern to type `Sₒᵤₜ`, optionally changing hermicity. When input and
output types match, only the hermicity is adjusted. Cross-type conversions (e.g.
`CSCRealSparsity → BlockRealSparsity`, `BlockRealSparsity → DenseRecipSparsity`) are also
supported.
"""
function convert_sparsity(
    ::Type{S}, in_sparsity::S, ::BasisSetMetadata;
    hermitian=false,
) where {S<:AbstractSparsity}
    if hermitian
        return convert_to_hermitian(in_sparsity)
    else
        return convert_to_nonhermitian(in_sparsity)
    end
end

function convert_sparsity(
    ::Type{Sₒᵤₜ},
    in_sparsity::CSCRealSparsity,
    basisset::BasisSetMetadata;
    hermitian=false,
) where {Sₒᵤₜ<:BlockRealSparsity}
    basis2atom = get_basis2atom(basisset)
    if in_sparsity.hermitian
        out_sparsity = BlockRealSparsity(
            in_sparsity.colcellptr, in_sparsity.rowval, in_sparsity.images, basis2atom
        )
    else
        throw(error("Currently, converting non-hermitian CSCRealSparsity to
        BlockRealSparsity is not supported."))
    end
    # Pass through `convert_sparsity` to convert hermicity
    return convert_sparsity(Sₒᵤₜ, out_sparsity, basisset; hermitian=hermitian)
end

# NOTE:
# - Tried using LittleSet for a 2 atom system and it was ~10 times slower.
#   For more atoms where unit cells are larger LittleSet might outperform regular Set.
#   One could choose the type of set to use based on the number of images inside `images`
# - When hermitian = true, ij2images keys (ij pairs) are already in the upper
#   triangle as required
# - Currently this function only works when the CSC sparsity is hermitian
function BlockRealSparsity(
    colcellptr::Array{T,3}, rowval::Vector{T}, images, basis2atom::Vector{T}
) where {T<:Integer}
    natoms = maximum(basis2atom)
    ij2images_set = [Set{SVector{3,Int}}() for _ in 1:natoms, _ in 1:natoms]

    # Could define csc iteration, but i_atom_col would have to be evaluated in inner loop
    for i_cell in axes(colcellptr, 2)
        image = images[i_cell]

        for i_basis_col in axes(colcellptr, 3)
            i_atom_col = basis2atom[i_basis_col]

            i_index_first = colcellptr[1, i_cell, i_basis_col]
            i_index_last = colcellptr[2, i_cell, i_basis_col]

            for i_index in i_index_first:i_index_last
                i_basis_row = rowval[i_index]
                i_atom_row = basis2atom[i_basis_row]

                push!(ij2images_set[i_atom_row, i_atom_col], image)
            end
        end
    end

    ij_nonempty_cartesian = findall(!isempty, ij2images_set)
    ij2images = Dictionary(
        Tuple.(ij_nonempty_cartesian),
        collect.(ij2images_set[ij_nonempty_cartesian]),
    )

    # Make on-site atom pair sparsity non-hermitian even when hermitian == true.
    make_onsite_nonhermitian!(ij2images)

    return BlockRealSparsity(ij2images, images, false)
end

function convert_sparsity(
    ::Type{<:DenseRecipSparsity}, in_sparsity::BlockRealSparsity, ::BasisSetMetadata;
    hermitian=false,
)
    return DenseRecipSparsity(hermitian)
end

### SPARSITY HERMICITY CONVERSIONS ###

function convert_to_nonhermitian(sparsity::DenseRecipSparsity)
    !op_hermicity(sparsity) && return sparsity

    return DenseRecipSparsity(false)
end

function convert_to_nonhermitian(sparsity::BlockRealSparsity)
    !op_hermicity(sparsity) && return sparsity

    # Shallow copy
    ij2images_nonherm = copy(sparsity.ij2images)

    upper_ij = filter(keys(sparsity.ij2images)) do ij
        i, j = ij
        i < j
    end

    # Insert ji pairs with -R images to the copy
    for (i, j) in upper_ij
        insert!(ij2images_nonherm, (j, i), -sparsity.ij2images[(i, j)])
    end

    return BlockRealSparsity(ij2images_nonherm, sparsity.images, false)
end

function convert_to_hermitian(sparsity::DenseRecipSparsity)
    op_hermicity(sparsity) && return sparsity

    return DenseRecipSparsity(true)
end

function convert_to_hermitian(sparsity::BlockRealSparsity)
    op_hermicity(sparsity) && return sparsity

    # Upper plus diagonal
    upperd_ij = filter(keys(sparsity.ij2images)) do ij
        i, j = ij
        i ≤ j
    end
    ij2images_herm = getindices(sparsity.ij2images, upperd_ij)

    return BlockRealSparsity(ij2images_herm, sparsity.images, true)
end

function make_onsite_nonhermitian!(ij2images)
    onsite_keys = filter(keys(ij2images)) do ij
        i, j = ij
        i == j
    end

    for ii in onsite_keys
        ij2images[ii] = union(ij2images[ii], -ij2images[ii])
    end
    return ij2images
end

### UTILITY FUNCTIONS ###

"""
    get_iglobal2ilocal(sparsity::BlockRealSparsity)

Return per-pair index maps from global image indices to local image indices in `ij2images`.
"""
function get_iglobal2ilocal(sparsity::BlockRealSparsity)
    return get_iexternal2ilocal(sparsity.images, sparsity)
end

function get_iexternal2ilocal(images, sparsity::BlockRealSparsity)
    return map(sparsity.ij2images) do images_local
        indexin(images, images_local)
    end
end

function get_iexternal2imlocal(images, sparsity::BlockRealSparsity)
    return map(sparsity.ij2images) do images_local
        indexin(images, -images_local)
    end
end

# local ⊆ global
function get_ilocal2iglobal(sparsity::BlockRealSparsity)
    ilocal2iexternal = get_ilocal2iexternal(sparsity.images, sparsity)
    return map(ilocal2iexternal) do ilocal2iexternal_ij
        convert(Vector{Int}, ilocal2iexternal_ij)
    end
end

function get_ilocal2iexternal(images, sparsity::BlockRealSparsity)
    return map(sparsity.ij2images) do images_local
        indexin(images_local, images)
    end
end

function get_imlocal2iexternal(images, sparsity::BlockRealSparsity)
    return map(sparsity.ij2images) do images_local
        indexin(-images_local, images)
    end
end

# Assuming sparsity.ij2images[(i, j)] for i == j will contain
# redundant images even when sparsity.hermitian == true
# ⟹ returned values are never of type Nothing
function get_onsite_ilocal2imlocal(sparsity::BlockRealSparsity)
    onsite_keys = filter(keys(sparsity.ij2images)) do ij
        i, j = ij
        i == j
    end

    # ij keys
    ilocal2imlocal = map(view(sparsity.ij2images, onsite_keys)) do images_local
        indexin(images_local, -images_local)
    end

    # i keys
    return Dictionary(
        values(map(first, keys(ilocal2imlocal))),
        values(ilocal2imlocal),
    )
end

# For a given iR in ij2images[ij] that corresponds to some image R,
# return miR in ij2images[ji] that corresponds to -R.
# Valid only if non-hermitian, otherwise all (iat, jat) pairs where
# iat ≠ jat will yield Nothing arrays
function get_ilocal2imlocal(sparsity::BlockRealSparsity)
    @argcheck !op_hermicity(sparsity)
    return map(keys(sparsity.ij2images)) do ij
        i, j = ij
        images_ijR = sparsity.ij2images[(i, j)]
        images_jimR = sparsity.ij2images[(j, i)]
        indexin(images_ijR, -images_jimR)
    end
end

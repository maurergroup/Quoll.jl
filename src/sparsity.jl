abstract type AbstractSparsity end

struct RealCSCSparsity <: AbstractSparsity
    rowval::Vector{Int}
    colcellptr::Array{Int, 3}
    images::Vector{SVector{3, Int}}
    hermitian::Bool
end

# In reciprocal space we could only have unique (i, j) pairs from `ij2images`
# Assuming 'U' hermicity
struct RealBlockSparsity <: AbstractSparsity
    ij2images::Dictionary{NTuple{2, Int}, Vector{SVector{3, Int}}}
    images::Vector{SVector{3, Int}}
    hermitian::Bool
end

### FROM NEIGHBOURLIST ###

# assumes fractional coordinates of all atoms are [0.0, 1.0)
# if T is unitless assumes angstrom
function RealBlockSparsity(atoms::AbstractSystem, radii::SpeciesAnyDict; hermitian = true)
    maxedges = get_maxedges(radii)
    nlist_cut = maximum(values(maxedges))
    nlist = PairList(atoms, nlist_cut)

    return RealBlockSparsity(atoms, nlist, maxedges; hermitian = hermitian)
end

function RealBlockSparsity(atoms::AbstractSystem, nlist::PairList, maxedges::SpeciesPairAnyDict; hermitian = true)
    unique_ij = unique(zip(nlist.i, nlist.j))
    if hermitian
        unique_ij = [(i, j) for (i, j) in unique_ij if j ≥ i]
    end

    ij2images = Dictionary(unique_ij, [SVector{3, Int}[] for ij in unique_ij])

    maxedges_nounits = Base.ImmutableDict((
        key => ustrip(value)
        for (key, value) in pairs(maxedges)
    )...)

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
        haskey(ij2images, ii) || insert!(ij2images, ii, SVector{3, Int}[])
        push!(ij2images[ii], SA[0, 0, 0])
    end

    return RealBlockSparsity(ij2images, unique(nlist.S), hermitian)
end

function get_maxedges(radii::SpeciesAnyDict{T}) where T<:Number
    angstrom = unit(T) == NoUnits ? u"Å" : one(T)
    maxedges = Dict(
        (zi, zj) => (i_cut + j_cut) * angstrom
        for (zi, i_cut) in pairs(radii)
        for (zj, j_cut) in pairs(radii)
    )
    return maxedges
end

### UTILITY FUNCTIONS ###

# A map from image indices in `images` to local image indices in `ij2images` for a given ij
function get_iglobal2ilocal(sparsity::RealBlockSparsity)
    return get_iexternal2ilocal(sparsity.images, sparsity)
end

function get_iexternal2ilocal(images, sparsity::RealBlockSparsity)
    return map(sparsity.ij2images) do images_local
        indexin(images, images_local)
    end
end

function get_iexternal2imlocal(images, sparsity::RealBlockSparsity)
    return map(sparsity.ij2images) do images_local
        indexin(images, -images_local)
    end
end

# Assuming sparsity.ij2images[(i, j)] for i == j will contain
# redundant images even when sparsity.hermitian == true
# => returned values are never of type Nothing
function get_onsite_ilocal2imlocal(sparsity::RealBlockSparsity)
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
        values(ilocal2imlocal)
    )
end

# For a given iR in ij2images[ij] that corresponds to some image R,
# return miR in ij2images[ji] that corresponds to -R.
# Valid only if non-hermitian, otherwise all (iat, jat) pairs where
# iat ≠ jat will yield Nothing arrays
function get_ilocal2imlocal(sparsity::RealBlockSparsity)
    @argcheck !sparsity.hermitian
    return map(keys(sparsity.ij2images)) do ij
        i, j = ij
        images_ijR = sparsity.ij2images[(i, j)]
        images_jimR = -sparsity.ij2images[(j, i)]
        indexin(images_ijR, images_jimR)
    end
end

### HERMICITY CONVERSIONS ###

# TODO: Not applicable anymore, we always assume it's type 1 sparsity
# Good that it's not inplace because otherwise the sparsity
# of the operator might not match its data.
# The function might look more sophisticated than it needs to be,
# but it might make more sense if one considers an example
# Non-hermitian ij2images:
# (i, j) => [R1, R2]
# (j, i) => [-R1, -R2]
# Hermitian ij2images type 1
# (i, j) => [R1, R2]
# <(j, i) => ...>: not present
# Hermitian ij2images type 2
# (i, j) => [R1]
# (j, i) => [-R2]
# Both possibilities are valid hermitian sparsities.
# 
# If converting from type 1 to nonhermitian sparsity then we will have new keys in ij2images.
# Does that affect how we convert other properties during operator -> operator conversion?
# This matters only for conversions where image order and (i, j) order matters in other operator
# fields, e.g. 3rd axis in keydata in RealBSparseOperator contains some order of images
# and by changing the sparsity the order would not be the same anymore.
# This would make RealBSparseOperator -> RealBSparseOperator conversion a bit tricky
# but it's not something that would have to be done normally.
function convert_to_nonhermitian(sparsity::RealBlockSparsity)
    !sparsity.hermitian && return sparsity

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

    return RealBlockSparsity(ij2images_nonherm, sparsity.images, false)
end

function convert_to_hermitian(sparsity::RealBlockSparsity)
    sparsity.hermitian && return sparsity

    # Upper plus diagonal
    upperd_ij = filter(keys(sparsity.ij2images)) do ij
        i, j = ij
        i ≤ j
    end
    ij2images_herm = getindices(sparsity.ij2images, upperd_ij)

    return RealBlockSparsity(ij2images_herm, sparsity.images, true)
end

function make_onsite_hermitian!(ij2images)
    onsite_keys = filter(keys(ij2images)) do ij
        i, j = ij
        i == j
    end

    for ii in onsite_keys
        ij2images[ii] = union(ij2images[ii], -ij2images[ii])
    end
    return ij2images
end

### SPARSITY CONVERSIONS ###

# Handles in_sparsity == out_sparsity conversions, only hermicity needs to be changed in such cases
function convert_sparsity(in_sparsity::T, ::Type{T}; hermitian = false) where T<:AbstractSparsity
    if hermitian
        return convert_to_hermitian(in_sparsity)
    else
        return convert_to_nonhermitian(in_sparsity)
    end
end

function convert_sparsity(in_sparsity::RealCSCSparsity, basisset::BasisSetMetadata, out_sparsity_type::Type{RealBlockSparsity};
    hermitian = true)
    basis2atom = get_basis2atom(basisset)
    out_sparsity = RealBlockSparsity(in_sparsity.colcellptr, in_sparsity.rowval, in_sparsity.images, basis2atom)
    return convert_sparsity(out_sparsity, out_sparsity_type, hermitian = hermitian)
end

# N.B. Tried using LittleSet for a 2 atom system and it was ~10 times slower.
# For more atoms where unit cells are larger LittleSet might outperform regular Set.
# One could choose the type of set to use based on the number of images inside `images`
function RealBlockSparsity(colcellptr::Array{T, 3}, rowval::Vector{T}, images, basis2atom::Vector{T}) where T<:Integer
    natoms = maximum(basis2atom)
    ij2images_set = [Set{SVector{3, Int}}() for _ in 1:natoms, _ in 1:natoms]

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
        collect.(ij2images_set[ij_nonempty_cartesian])
    )

    # Make on-site blocks non-hermitian even when hermitian == true
    make_onsite_hermitian!(ij2images)
    
    return RealBlockSparsity(ij2images, images, true)
end
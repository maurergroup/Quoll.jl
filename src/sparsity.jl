using NeighbourLists: PairList, pairs

abstract type AbstractSparsity end

struct RealCSCSparsity <: AbstractSparsity
    rowval::Vector{Int}
    colcellptr::Array{Int, 3}
    images::Vector{SVector{3, Int}}
    hermitian::Bool
end

# TODO: this makes sense only in real space
# In reciprocal space we could only have unique (i, j) pairs from `ij2images`
# TODO: hermitian could be either `U` or `L`
struct RealBlockSparsity <: AbstractSparsity
    ij2images::Dictionary{Tuple{Int, Int}, Vector{SVector{3, Int}}}
    images::Vector{SVector{3, Int}}
    hermitian::Bool
end

function get_maxedges(radii::Dict{ChemicalSpecies, T}) where T<:Number
    angstrom = unit(T) == NoUnits ? u"Å" : one(T)
    maxedges = Dict(
        (zi, zj) => (i_cut + j_cut) * angstrom
        for (zi, i_cut) in radii
        for (zj, j_cut) in radii
    )
    return maxedges
end

# assumes fractional coordinates of all atoms are [0.0, 1.0)
# if T is unitless assumes angstrom
function RealBlockSparsity(atoms::AbstractSystem, radii::Dict{ChemicalSpecies}; hermitian = true)
    maxedges = get_maxedges(radii)
    nlist_cut = maximum(values(maxedges))
    nlist = PairList(atoms, nlist_cut)

    return RealBlockSparsity(atoms, nlist, maxedges; hermitian = hermitian)
end

function RealBlockSparsity(atoms::AbstractSystem, nlist::PairList, maxedges::Dict{Tuple{ChemicalSpecies, ChemicalSpecies}}; hermitian = true)
    # Dictionary (i, j) => images
    unique_ij = unique(zip(nlist.i, nlist.j))
    if hermitian
        unique_ij = [(i, j) for (i, j) in unique_ij if j ≥ i]
    end

    ij2images = Dictionary(unique_ij, [SVector{3, Int}[] for ij in unique_ij])

    for (pair, image) in zip(pairs(nlist), nlist.S)
        i, j = pair[1], pair[2]
        hermitian && i > j && continue

        # Sift edges which are beyond the cutoff value
        r = norm(pair[3]) * u"Å"
        if r ≤ maxedges[(species(atoms, i), species(atoms, j))]
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

function RealBlockSparsity(csc_sparsity::RealCSCSparsity, basisset::BasisSetMetadata)
    return RealBlockSparsity(csc_sparsity.colcellptr, csc_sparsity.rowval, csc_sparsity.images, basisset.basis2atom, csc_sparsity.hermitian)
end

# TODO: Does this always return fully non-hermitian (iat, iat)? It doesn't but make sure it does
function RealBlockSparsity(colcellptr, rowval, images, basis2atom, hermitian)
    ij2images = Dictionary{Tuple{Int, Int}, Vector{SVector{3, Int}}}()

    # TODO: could define csc iteration, but i_atom_col would have to be evaluated in inner loop
    for i_cell in axes(colcellptr, 2)
        image = images[i_cell]

        for i_basis_col in axes(colcellptr, 3)
            i_atom_col = basis2atom[i_basis_col]

            i_index_first = colcellptr[1, i_cell, i_basis_col]
            i_index_last = colcellptr[2, i_cell, i_basis_col]

            for i_index in i_index_first:i_index_last
                i_basis_row = rowval[i_index]
                i_atom_row = basis2atom[i_basis_row]

                # If ij doesn't exist in ij2images create an empty array
                ij = (i_atom_row, i_atom_col)
                haskey(ij2images, ij) || insert!(ij2images, ij, SVector{3, Int}[])

                # If the image doesn't exist in ij2images[ij] then push the image inside
                image ∈ ij2images[ij] || push!(ij2images[ij], image)
            end
        end
    end
    return RealBlockSparsity(ij2images, images, hermitian)
end

# A map from image indices in `images` to local image indices in `ij2images` for a given ij
function get_iglobal2ilocal(sparsity::RealBlockSparsity)
    return Dictionary(
        keys(sparsity.ij2images),
        [indexin(sparsity.images, images_local) for images_local in sparsity.ij2images],
    )
end

function get_onsite_ilocal2milocal(sparsity::RealBlockSparsity)
    return Dictionary(
        [iat for (iat, jat) in keys(sparsity.ij2images) if iat == jat],
        [indexin(images_local, -images_local)
        for ((iat, jat), images_local) in pairs(sparsity.ij2images)
        if iat == jat]
    )
end

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

    ij_list = Tuple{Int, Int}[]
    images_list = Vector{SVector{3, Int}}[]

    unique_sorted_ij = filter(keys(sparsity.ij2images)) do ij
        i, j = ij
        i ≤ j
    end

    for (i, j) in unique_sorted_ij
        if i ≠ j
            ij_images = union(
                get(() -> SVector{3, Int}[], sparsity.ij2images, (i, j)),
                (-1 * get(() -> SVector{3, Int}[], sparsity.ij2images, (j, i)))
            )
            ji_images = -1 * ij_images

            push!(ij_list, (i, j), (j, i))
            push!(images_list, ij_images, ji_images)
        else
            ii_images = get(() -> SVector{3, Int}[], sparsity.ij2images, (i, i))
            
            push!(ij_list, (i, i))
            push!(images_list, ii_images)
        end
    end
    nonhermitian_ij2images = Dictionary(ij_list, images_list)

    return RealBlockSparsity(nonhermitian_ij2images, sparsity.images, false)
end

# Converting to "type 1" hermitian sparsity (for more information see `convert_to_nonhermitian`)
# Same as for `convert_to_nonhermitian`, some of the keys will be deleted here.
function convert_to_hermitian(sparsity::RealBlockSparsity)
    sparsity.hermitian && return sparsity

    unique_sorted_ij = filter(keys(sparsity.ij2images)) do ij
        i, j = ij
        i ≤ j
    end
    hermitian_ij2images = Dictionary(view(sparsity.ij2images, unique_sorted_ij))

    return RealBlockSparsity(hermitian_ij2images, sparsity.images, true)
end
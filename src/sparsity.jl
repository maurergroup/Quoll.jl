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
    for (i, _) in enumerate(atoms)
        ii = (i, i)

        # Add empty array if the key doesn't exist yet
        # (would be in the case of no i self-interactions across periodic boundaries)
        ii ∈ keys(ij2images) || insert!(ij2images, ii, SVector{3, Int}[])
        push!(ij2images[ii], SA[0, 0, 0])
    end

    return RealBlockSparsity(ij2images, nlist.S, hermitian)
end

function RealBlockSparsity(csc_sparsity::RealCSCSparsity, basisset::BasisSetMetadata)
    return RealBlockSparsity(csc_sparsity.colcellptr, csc_sparsity.rowval, csc_sparsity.images, basisset.basis2atom)
end

function RealBlockSparsity(colcellptr, rowval, images, basis2atom)
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
                ij ∈ keys(ij2images) || insert!(ij2images, ij, SVector{3, Int}[])

                # If the image doesn't exist in ij2images[ij] then push the image inside
                image ∈ ij2images[ij] || push!(ij2images[ij], image)
            end
        end
    end
    return RealBlockSparsity(ij2images, images, true)
end

# A map from image indices in `images` to local image indices in `ij2images` for a given ij
function get_iglobal2ilocal(sparsity::RealBlockSparsity)
    return Dictionary(
        keys(sparsity.ij2images),
        [indexin(sparsity.images, images_local) for images_local in sparsity.ij2images],
    )
end
using NeighbourLists: PairList, pairs

abstract type AbstractSparsity end

struct RealCSCSparsity <: AbstractSparsity
    rowval::Vector{Int}
    colcellptr::Array{Int, 3}
    images::Vector{SVector{3, Int}}
end

struct BlockSparsity <: AbstractSparsity
    ij2images::Dictionary{Tuple{Int, Int}, Vector{SVector{3, Int}}}
    images::Vector{SVector{3, Int}}
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
function BlockSparsity(atoms::AbstractSystem, radii::Dict{ChemicalSpecies})
    maxedges = get_maxedges(radii)
    nlist_cut = maximum(values(maxedges))
    nlist = PairList(atoms, nlist_cut)

    return BlockSparsity(nlist, maxedges)
end

function BlockSparsity(atoms::AbstractSystem, nlist::PairList, maxedges::Dict{Tuple{ChemicalSpecies, ChemicalSpecies}})
    # Dictionary (i, j) => images
    unique_ij = unique(zip(nlist.i, nlist.j))
    ij2images = Dictionary(unique_ij, [SVector{3, Int}[] for ij in unique_ij])

    for (pair, image) in zip(pairs(nlist), nlist.S)
        iat, jat = pair[1], pair[2]
        ij = (iat, jat)

        # Sift edges which are beyond the cutoff value
        r = norm(pair[3])
        if r <= maxedges[(species(atoms, iat), species(atoms, jat))]
            push!(ij2images[ij], image)
        end
    end

    # Add self-interactions for R = [0, 0, 0]
    for (iat, _) in enumerate(atoms)
        ii = (iat, iat)
        push!(ij2images[ii], SA[0, 0, 0])
    end

    return BlockSparsity(ij2images, nlist.S)
end
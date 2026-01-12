@enum Spin begin
    ⬆ = 1
    ⬇ = -1
end

struct SpinsMetadata{S<:SpeciesAnyDict}
    spins::S
    soc::Bool
end

function spins_species(spins::SpinsMetadata, z::ChemicalSpecies)
    return spins.spins[z]
end

function SpinsMetadata(
    source::AbstractSource, kind::OperatorKind, basisset::BasisSetMetadata
)
    return SpinsMetadata(source, kind, basisset.basis)
end

function SpinsMetadata(
    source::AbstractSource, kind::Hamiltonian, basis::SpeciesAnyDict
)
    return SpinsMetadata(source, Val(kind.spin), basis)
end

function SpinsMetadata(
    ::AbstractSource, ::Val{:up}, basis::SpeciesAnyDict
)
    spins_dict = Base.ImmutableDict(
        (
            z => fill(⬆, length(atom_basis))
            for (z, atom_basis) in pairs(basis)
        )...,
    )
    return SpinsMetadata(spins_dict, false)
end

function SpinsMetadata(
    ::AbstractSource, ::Val{:down}, basis::SpeciesAnyDict
)
    spins_dict = Base.ImmutableDict(
        (
            z => fill(⬇, length(atom_basis))
            for (z, atom_basis) in pairs(basis)
        )...,
    )
    return SpinsMetadata(spins_dict, false)
end

function SpinsMetadata(
    ::AbstractSource, ::Val{:soc}, basis::SpeciesAnyDict
)
    spins_dict = Base.ImmutableDict(
        (
            z => vcat(fill(⬆, length(atom_basis))..., fill(⬇, length(atom_basis))...)
            for (z, atom_basis) in pairs(basis)
        )...,
    )
    return SpinsMetadata(spins_dict, true)
end

function convert_spins_shconv(
    spins::SpinsMetadata, basisset::BasisSetMetadata, shconv::SHConvention
)
    in_spins = spins.spins
    in_basis = basisset.basis
    out_spins = convert_speciesdict_shconv(in_spins, in_basis, shconv)
    return SpinsMetadata(out_spins, spins.soc)
end

function reduce_spins(
    spins::SpinsMetadata, basisset::BasisSetMetadata, subbasis::Vector{<:BasisMetadata};
    inverted=false,
)
    subbasis_masks = get_subbasis_masks(basisset, subbasis; inverted=inverted)
    return reduce_spins(spins, subbasis_masks)
end

# Make a subspins based on the masks for each chemical species
function reduce_spins(
    spins::SpinsMetadata, subbasis_masks::SpeciesDictionary{BitVector}
)
    subspins_dict = getindex.(spins.spins, subbasis_masks)
    return SpinsMetadata(subspins_dict, spins.soc)
end

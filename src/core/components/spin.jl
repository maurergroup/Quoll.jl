"""
    Spin

Enum with values `⬆` (+1) and `⬇` (−1) representing spin-up and spin-down channels.
"""
@enum Spin begin
    ⬆ = 1
    ⬇ = -1
end

"""
    SpinsMetadata{S}

Per-species spin assignment for each basis function.

Fields:
- `spins::S` — species-keyed dictionary mapping each species to a vector of [`Spin`](@ref)
  values (one per basis function). For SOC, each species' vector is doubled.
- `soc::Bool` — `true` if this represents spin-orbit coupling (basis doubled per species).
"""
struct SpinsMetadata{S<:SpeciesAnyDict}
    spins::S
    soc::Bool
end

function spins_species(spins::SpinsMetadata, z::ChemicalSpecies)
    return spins.spins[z]
end

"""
    SpinsMetadata(source, kind, basisset)

Construct `SpinsMetadata` from an operator kind's `:spin` tag. Dispatches on the spin tag
value (`:up`, `:down`, `:soc`) to build the appropriate per-species spin vectors.
"""
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
            z => vcat(
                fill(⬆, round(Int, length(atom_basis) / 2)),
                fill(⬇, round(Int, length(atom_basis) / 2)),
            )
            for (z, atom_basis) in pairs(basis)
        )...,
    )
    return SpinsMetadata(spins_dict, true)
end

"""
    convert_spins_shconv(spins, basisset, shconv)

Reorder spin vectors to match a new SH convention, preserving the spin-basis correspondence.
"""
function convert_spins_shconv(
    spins::SpinsMetadata, basisset::BasisSetMetadata, shconv::SHConvention
)
    in_spins = spins.spins
    in_basis = basisset.basis
    out_spins = convert_speciesdict_shconv(in_spins, in_basis, shconv)
    return SpinsMetadata(out_spins, spins.soc)
end

"""
    reduce_spins(spins, basisset, subbasis; inverted=false)

Restrict spin metadata to a sub-basis, keeping only the spin entries corresponding to the
retained basis functions.
"""
function reduce_spins(
    spins::SpinsMetadata, basisset::BasisSetMetadata, subbasis::Vector{<:BasisMetadata};
    inverted=false,
)
    subbasis_masks = get_subbasis_masks(basisset, subbasis; inverted=inverted)
    return reduce_spins(spins, subbasis_masks)
end

# Make a subspins based on the masks for each chemical species
function reduce_spins(
    spins::SpinsMetadata{<:Base.ImmutableDict}, subbasis_masks::SpeciesAnyDict{BitVector}
)
    subspins_dict = Base.ImmutableDict((
        z => getindex(spins.spins[z], subbasis_masks[z])
        for z in keys(subbasis_masks)
    )...)
    return SpinsMetadata(subspins_dict, spins.soc)
end

"""
    convert_spins_source(in_spins, out_basisset, in_source, out_source)

Make final changes due to the source change (e.g. reorder up and down spins if the two
sources don't agree). This often might leave the spins unchanged.
"""
function convert_spins_source end
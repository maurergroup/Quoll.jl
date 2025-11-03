### TYPE ALIASES ###

const AbstractAnyDict{K, V} = Union{<:AbstractDict{K, V}, <:AbstractDictionary{K, V}}

const SpeciesPairDict{T, N} = AbstractAnyDict{NTuple{2, ChemicalSpecies}, <:AbstractArray{T, N}}
const SpeciesDict{T, N} = AbstractAnyDict{ChemicalSpecies, <:AbstractArray{T, N}}
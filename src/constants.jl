### TYPE ALIASES ###

const AbstractAnyDict{K, V} = Union{<:AbstractDict{K, V}, <:AbstractDictionary{K, V}}
const SpeciesPairAnyDict{V} = AbstractAnyDict{NTuple{2, ChemicalSpecies}, V}
const SpeciesAnyDict{V}  = AbstractAnyDict{ChemicalSpecies, V}
const AtomPairAnyDict{V} = AbstractAnyDict{NTuple{2, Int}, V}

const SpeciesDictionary{V} = AbstractDictionary{ChemicalSpecies, V}
const SpeciesDict{V} = AbstractDict{ChemicalSpecies, V}

const ThreeDimSliceView{T} = SubArray{
    T, 3, Array{T, 3},
    Tuple{
        Base.Slice{Base.OneTo{Int}},
        Base.Slice{Base.OneTo{Int}},
        UnitRange{Int}
    }, true
}

### UNITS ###
# Internally UnitfulAtomic is used to convert between units if needed,
# except during reading and writing. Units defined here are used instead
# to make sure conversions between different codes return a consistent result

const HARTREE_CODATA_2002 = 27.2113845 # eV
const BOHR_RADIUS_CODATA_2002 = 0.52917721 # Å
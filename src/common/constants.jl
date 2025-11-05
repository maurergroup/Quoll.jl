### TYPE ALIASES ###

const AbstractAnyDict{K, V} = Union{<:AbstractDict{K, V}, <:AbstractDictionary{K, V}}

const SpeciesPairAnyDict{V} = AbstractAnyDict{NTuple{2, ChemicalSpecies}, V}
const SpeciesAnyDict{V}  = AbstractAnyDict{ChemicalSpecies, V}
# const SpeciesDict{V} = AbstractDict{ChemicalSpecies, V}
# const SpeciesDictionary{V} = AbstractDictionary{ChemicalSpecies, V}

const AtomPairAnyDict{V} = AbstractAnyDict{NTuple{2, Int}, V}
# const AtomPairKeyedArray{T, N, AT, KT} = AbstractArray{<:KeyedArray{T, N, AT, KT}, 2}

const ThreeDimSliceView{T} = SubArray{
    T, 3, Array{T, 3},
    Tuple{
        Base.Slice{Base.OneTo{Int}},
        Base.Slice{Base.OneTo{Int}},
        UnitRange{Int}
    }, true
}
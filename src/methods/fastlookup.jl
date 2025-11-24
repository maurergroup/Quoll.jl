function convert_to_atomarray(d::SpeciesAnyDict{V}, atom2species) where V
    natoms = length(atom2species)
    sentinel = get_sentinel_value(V)
    return [
        get(d, atom2species[i], sentinel)
        for i in 1:natoms
    ]
end

function convert_to_atomarray(d::SpeciesPairAnyDict{V}, atom2species) where V
    natoms = length(atom2species)
    sentinel = get_sentinel_value(V)
    return [
        get(d, (atom2species[i], atom2species[j]), sentinel)
        for i in 1:natoms, j in 1:natoms
    ]
end

function convert_to_atomarray(d::AtomPairAnyDict{V}, natoms) where V
    sentinel = get_sentinel_value(V)
    return [
        get(d, (i, j), sentinel)
        for i in 1:natoms, j in 1:natoms
    ]
end

get_sentinel_value(::Type{V}) where V = missing

get_sentinel_value(::Type{A}) where A<:AbstractArray = A(undef, Tuple(fill(0, ndims(A)))...)

get_sentinel_value(::Type{A}) where A<:KeyedArray = missing

# This is for keydata.data.data, but I found that having Union{missing, T}
# type instability doesn't lead to measurable performance loss,
# so at least for keydata it seems like it's not really worth it
function get_sentinel_value(::Type{<:ThreeDimSliceView{T}}) where T
    dummy_array = zeros(T, 1, 1, 1)
    return view(dummy_array, :, :, 0:-1)
end

# Could write a method to construct ij2offset::Matrix{UnitRange{Int}}
# for ultimate access speed if z1z2_ij2offset was ever needed in a hot loop
# and dictionary access becomes a bottleneck
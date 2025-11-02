using Base
using OffsetArrays

# Conversion from wiki to a given format
struct SHConversion{T}
    orders::T
    shifts::T
    phases::T
    lmax::Int
end

# E.g.
# [[1], [3,  1,  2], [1,  4,  5,  3,  2], ...]
# [[1], [1, -1,  1], [1,  1, -1, -1,  1], ...]
function SHConversion(orders, phases)
    @argcheck length(phases) == length(orders)
    lmax = length(phases) - 1

    # change from 1 start to -m start
    _orders = tuple((
        offset_static(order .- (l + 1), l)
        for (l, order) in zip(0:lmax, orders)
    )...)
    _shifts = compute_shifts(_orders)

    _phases = tuple((
        offset_static(phase, l)
        for (l, phase) in zip(0:lmax, phases)
    )...)

    return SHConversion(_orders, _shifts, _phases, lmax)
end

function compute_shifts(orders)
    return tuple((
        offset_static(OffsetArrays.no_offset_view(order) .- collect(-l:l), l)
        for (l, order) in zip(0:(length(orders) - 1), orders)
    )...)
end

function Base.inv(shconv::SHConversion)
    lmax = shconv.lmax
    inv_orders = tuple((
        inv_order(shconv.orders, l)
        for l in 0:lmax
    )...)
    inv_shifts = compute_shifts(inv_orders)

    inv_phases = tuple((
        offset_static(shconv.phases[l + 1][inv_order(shconv.orders, l)], l)
        for l in 0:lmax
    )...)
    return SHConversion(inv_orders, inv_shifts, inv_phases, lmax)
end

function inv_order(orders, l)
    return OffsetArray(SVector{2l + 1, Int}(invperm(OffsetArrays.no_offset_view(orders[l + 1]) .+ (l + 1)) .- (l + 1)), -l:l)
end

# Function which converts the vector into a static vector and offsets it to start
# from -m (instead of 1)
# (no_offset_view does nothing to normal vectors, I use this here because
# sometimes the vector will already be an offset vector)
function offset_static(vec, l)
    return OffsetArray(SVector{2l + 1, Int}(OffsetArrays.no_offset_view(vec)), -l:l)
end

# SHConvention is not callable and the result doesn't return ComposedFunction,
# so not sure if it's bad to overwrite this or not
function Base.:(∘)(shconv2::T1, shconv1::T2) where {T1<:SHConversion, T2<:SHConversion}
    if shconv1.lmax != shconv2.lmax
        @debug "lmax of two SH conventions do not match, setting lmax to the lower value"
        shconv1.lmax > shconv2.lmax ? (lmax = shconv2.lmax) : (lmax = shconv1.lmax)
    else
        lmax = shconv1.lmax
    end

    orders = tuple((
        offset_static(shconv1.orders[l + 1][shconv2.orders[l + 1]], l)
        for l in 0:lmax
    )...)
    shifts = compute_shifts(orders)

    phases = tuple((
        offset_static(shconv2.phases[l + 1] .* shconv1.phases[l + 1][shconv2.orders[l + 1]], l)
        for l in 0:lmax
    )...)

    return SHConversion(orders, shifts, phases, lmax)
end

get_shift(orbital::BasisMetadata, shconv::SHConversion) = get_shift(orbital.l, orbital.m, shconv)
get_phase(orbital::BasisMetadata, shconv::SHConversion) = get_phase(orbital.l, orbital.m, shconv)

get_shift(l, m, shconv::SHConversion) = shconv.orders[l + 1][m] - m
get_phase(l, m, shconv::SHConversion) = shconv.phases[l + 1][m]

function precompute_shifts(basisset::BasisSetMetadata, shconv::SHConversion)
    shifts = Dictionary{ChemicalSpecies, Vector{Int}}()
    for z in get_unique_species(basisset)
        insert!(shifts, z, get_shift.(basis_species(basisset, z), Ref(shconv)))
    end
    return shifts
end

# Assuming indexing starts at 1
function precompute_orders(basisset::BasisSetMetadata, shconv::SHConversion)
    shifts = precompute_shifts(basisset, shconv)
    orders = Dictionary{ChemicalSpecies, Vector{Int}}()
    for z in get_unique_species(basisset)
        nbasis = length(basis_species(basisset, z))
        insert!(orders, z, collect(1:nbasis) .+ shifts[z])
    end
    return orders
end

function precompute_phases(basisset::BasisSetMetadata, shconv::SHConversion; DIM = 1)
    phases = Dictionary{ChemicalSpecies, Vector{Int}}()
    for z in get_unique_species(basisset)
        insert!(phases, z, get_phase.(basis_species(basisset, z), Ref(shconv)))
    end
    precompute_phases!(phases, Val(DIM))
end

precompute_phases!(phases, ::Val{1}) = phases

function precompute_phases!(phases, ::Val{2})
    species_tuples = tuple.(keys(phases).values)

    combine_tuples(t1, t2) = tuple(t1..., t2...)
    outer_product = combine_tuples.(species_tuples, reshape(species_tuples, 1, :))

    species_pair_keys = Indices(reshape(outer_product, :))

    return map(species_pair_keys) do key
        z1, z2 = key
        phases[z1] * phases[z2]'
    end
end
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

function get_shift(l, m, shconv::SHConversion)
    return shconv.orders[l + 1][m] - m
end

function get_phase(l, m, shconv::SHConversion)
    return shconv.phases[l + 1][m]
end

function get_shiftphase(basis::BasisMetadata, shconv::SHConversion)
    return get_shift(basis.l, basis.m, shconv), get_phase(basis.l, basis.m, shconv)
end

function precompute_shiftphases(basisset::BasisSetMetadata, shconv::SHConversion)
    shifts = Dictionary{ChemicalSpecies, Vector{Int}}()
    phases = Dictionary{ChemicalSpecies, Vector{Int}}()

    for z in unique(basisset.atom2species)
        n_basis = length(basis_species(basisset, z))
        insert!(shifts, z, Vector{Int}(undef, n_basis))
        insert!(phases, z, Vector{Int}(undef, n_basis))

        for (ib, orbital) in enumerate(basisset.basis[z])
            shifts[z][ib], phases[z][ib] = get_shiftphase(orbital, shconv)
        end
    end

    return shifts, phases
end

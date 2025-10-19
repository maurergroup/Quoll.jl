using Base

# Conversion with respect to wiki
struct SHConversion{T}
    orders::T
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
    _phases = tuple((
        offset_static(phase, l)
        for (l, phase) in zip(0:lmax, phases)
    )...)

    return SHConversion(_orders, _phases, lmax)
end

function Base.inv(shconv::SHConversion)
    lmax = shconv.lmax
    inv_orders = tuple((
        inv_order(shconv.orders, l)
        for l in 0:lmax
    )...)
    inv_phases = tuple((
        offset_static(shconv.phases[l + 1][inv_order(shconv.orders, l)], l)
        for l in 0:lmax
    )...)
    return SHConversion(inv_orders, inv_phases, lmax)
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

function Base.:(∘)(shconv2::T, shconv1::T) where T<:SHConversion
    if shconv1.lmax != shconv2.lmax
        @warn "lmax of two SH conventions do not match, setting lmax to the lower value"
        shconv1.lmax > shconv2.lmax ? (lmax = shconv1.lmax) : (lmax = shconv2.lmax)
    else
        lmax = shconv1.lmax
    end

    orders = tuple((
        offset_static(shconv1.orders[l + 1][shconv2.orders[l + 1]], l)
        for l in 0:lmax
    )...)
    phases = tuple((
        offset_static(shconv2.phases[l + 1] .* shconv1.phases[l + 1][shconv2.orders[l + 1]], l)
        for l in 0:lmax
    )...)

    return SHConversion(orders, phases, lmax)
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

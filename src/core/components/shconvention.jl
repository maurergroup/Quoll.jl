# Convention defined as a conversion from wiki to a given format
struct SHConvention{lmax,T}
    orders::T
    shifts::T
    phases::T

    function SHConvention{lmax}(orders::T, shifts::T, phases::T) where {lmax,T}
        return new{lmax,T}(orders, shifts, phases)
    end
end

# E.g.
# [[1], [3,  1,  2], [1,  4,  5,  3,  2], ...]
# [[1], [1, -1,  1], [1,  1, -1, -1,  1], ...]
function SHConvention(orders, phases)
    @argcheck length(phases) == length(orders)
    lmax = length(phases) - 1

    _orders = tuple((make_static(order, Val(l)) for (l, order) in zip(0:lmax, orders))...)
    _phases = tuple((make_static(phase, Val(l)) for (l, phase) in zip(0:lmax, phases))...)

    _shifts = compute_shifts(_orders)

    return SHConvention{lmax}(_orders, _shifts, _phases)
end

function isidentity(shconv::SHConvention)
    for (orders_l, phases_l) in zip(shconv.orders, shconv.phases)
        if !issorted(orders_l) || any(phases_l .< 0)
            return false
        end
    end
    return true
end

function compute_shifts(orders)
    lmax = length(orders) - 1
    return tuple(
        (
            make_static(order .- collect(1:(2l + 1)), Val(l))
            for (l, order) in zip(0:lmax, orders)
        )...,
    )
end

function Base.inv(shconv::SHConvention{lmax}) where {lmax}
    inv_orders = tuple(
        (
            make_static(invperm(order), Val(l))
            for (l, order) in zip(0:lmax, shconv.orders)
        )...,
    )
    inv_shifts = compute_shifts(inv_orders)

    inv_phases = tuple(
        (
            make_static(phase[invperm(order)], Val(l))
            for (l, phase, order) in zip(0:lmax, shconv.phases, shconv.orders)
        )...,
    )

    return SHConvention{lmax}(inv_orders, inv_shifts, inv_phases)
end

function make_static(vec, ::Val{l}) where {l}
    return SVector{2l + 1,Int}(vec)
end

# SHConvention is not callable and the result doesn't return ComposedFunction,
# so not sure if it's a good idea to overload the composition operator
function Base.:(∘)(
    shconv2::T1, shconv1::T2
) where {T1<:SHConvention{lmax1},T2<:SHConvention{lmax2}} where {lmax1,lmax2}
    if lmax1 != lmax2
        @debug "lmax of two SH conventions do not match, setting lmax to the lower value"
        lmax = min(lmax1, lmax2)
    else
        lmax = lmax1
    end

    orders = tuple(
        (
            make_static(order1[order2], Val(l))
            for (l, order1, order2) in zip(0:lmax, shconv1.orders, shconv2.orders)
        )...,
    )
    shifts = compute_shifts(orders)

    phases = tuple(
        (
            make_static(phase2 .* phase1[order2], Val(l))
            for (l, phase1, phase2, order2)
            in
            zip(0:lmax, shconv1.phases, shconv2.phases, shconv2.orders)
        )...,
    )

    # return SHConvention(orders, shifts, phases, lmax)
    return SHConvention{lmax}(orders, shifts, phases)
end

# ib_sub: basis function index in a subblock
get_order(l, ib_sub, shconv::SHConvention) = shconv.orders[l + 1][ib_sub]
get_shift(l, ib_sub, shconv::SHConvention) = shconv.shifts[l + 1][ib_sub]
get_phase(l, ib_sub, shconv::SHConvention) = shconv.phases[l + 1][ib_sub]

function precompute_shifts(basis::SpeciesAnyDict, shconv::SHConvention)
    shifts = Dictionary{ChemicalSpecies,Vector{Int}}()
    for z in keys(basis)
        l_list = map(b -> b.l, basis[z])
        ib_sub_list = get_indices_in_subblock(basis[z])
        insert!(shifts, z, get_shift.(l_list, ib_sub_list, Ref(shconv)))
    end
    return Base.ImmutableDict(pairs(shifts)...)
end

# Orders for the whole basis set, contrary to `get_order` which computes orders
# only in a given subblock (which is why shifts are used instead of `get_order` directly)
function precompute_orders(basis::SpeciesAnyDict, shconv::SHConvention)
    orders = Dictionary{ChemicalSpecies,Vector{Int}}()
    for z in keys(basis)
        l_list = map(b -> b.l, basis[z])
        ib_sub_list = get_indices_in_subblock(basis[z])

        shifts = get_shift.(l_list, ib_sub_list, Ref(shconv))
        insert!(orders, z, shifts .+ collect(eachindex(basis[z])))
    end
    return Base.ImmutableDict(pairs(orders)...)
end

function precompute_shphases(
    basis::SpeciesAnyDict, shconv::SHConvention,
    ::Val{D}=Val(1), ::Type{T}=Float64,
) where {D,T}
    phases = Dictionary{ChemicalSpecies,Vector{T}}()
    for z in keys(basis)
        l_list = map(b -> b.l, basis[z])
        ib_sub_list = get_indices_in_subblock(basis[z])
        insert!(phases, z, get_phase.(l_list, ib_sub_list, Ref(shconv)))
    end
    return precompute_shphases!(phases, Val(D))
end

precompute_shphases!(phases::AbstractDictionary, ::Val{1}) =
    Base.ImmutableDict(pairs(phases)...)

function precompute_shphases!(phases::AbstractDictionary, ::Val{2})
    species_tuples = tuple.(keys(phases).values)

    combine_tuples(t1, t2) = tuple(t1..., t2...)
    outer_product = combine_tuples.(species_tuples, reshape(species_tuples, 1, :))

    species_pair_keys = Indices(reshape(outer_product, :))
    phases_pairs = map(species_pair_keys) do key
        z1, z2 = key
        phases[z1] * phases[z2]'
    end

    return Base.ImmutableDict(pairs(phases_pairs)...)
end

# SH conversion for a matrix M_z₁z₂ can be performed as P_z₁ * M_z₁z₂ * (P_z₂)ᵀ
function precompute_signed_perm_matrices(
    basis::SpeciesAnyDict, shconv::SHConvention, ::Type{T}=Float64
) where {T}
    signed_perms = Dictionary{ChemicalSpecies,Matrix{T}}()
    for z in keys(basis)
        insert!(signed_perms, z, compute_signed_perm_matrix(basis[z], shconv, T))
    end
    return signed_perms
end

function compute_signed_perm_matrix(
    basis_z, shconv::SHConvention, ::Type{T}=Float64
) where {T}
    perm_matrix = compute_perm_matrix(basis_z, shconv, T)
    l_list = map(b -> b.l, basis_z)
    ib_sub_list = get_indices_in_subblock(basis_z)
    phases = get_phase.(l_list, ib_sub_list, Ref(shconv))
    return Matrix{T}(Diagonal(phases) * perm_matrix)
end

function compute_perm_matrix(
    basis_z, shconv::SHConvention, ::Type{T}=Float64
) where {T}
    nbasis = length(basis_z)
    l_list = map(b -> b.l, basis_z)
    ib_sub_list = get_indices_in_subblock(basis_z)

    shifts = get_shift.(l_list, ib_sub_list, Ref(shconv))
    orders = shifts .+ collect(eachindex(basis_z))

    return Matrix{T}(I(nbasis)[orders, :])
end

function convert_speciesdict_shconv(
    in_d::SpeciesAnyDict, basis::SpeciesAnyDict, shconv::SHConvention
)
    orders = precompute_orders(basis, shconv)
    out_d = deepcopy(in_d)
    for z in keys(out_d)
        permute!(out_d[z], orders[z])
    end
    return out_d
end

function convert_speciesdict_shconv!(
    in_d::SpeciesAnyDict, basis::SpeciesAnyDict, shconv::SHConvention
)
    orders = precompute_orders(basis, shconv)
    for z in keys(in_d)
        permute!(in_d[z], orders[z])
    end
end
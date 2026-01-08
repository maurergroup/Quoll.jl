struct OperatorKind{K}
    tags::UnorderedDictionary{Symbol,Symbol}

    function OperatorKind{K}(tags) where {K}
        @argcheck K isa Symbol
        return new{K}(tags)
    end
end

OperatorKind{K}(pairs::Vararg{P}) where {K,P<:Pair} =
    OperatorKind{K}(UnorderedDictionary(dictionary(pairs)))
OperatorKind{K}(; kwargs...) where {K} =
    OperatorKind{K}(UnorderedDictionary(dictionary(kwargs)))

function Base.hash(op::OperatorKind{K}, h::UInt) where {K}
    return Base.hash(
        getfield(op, :tags),
        Base.hash(K, Base.hash(:OperatorKind, h)),
    )
end

function Base.:(==)(op1::OperatorKind{K₁}, op2::OperatorKind{K₂}) where {K₁,K₂}
    return K₁ == K₂ && getfield(op1, :tags) == getfield(op2, :tags)
end

function Base.isequal(op1::OperatorKind{K₁}, op2::OperatorKind{K₂}) where {K₁,K₂}
    return isequal(K₁, K₂) && isequal(getfield(op1, :tags), getfield(op2, :tags))
end

Base.haskey(op::OperatorKind, key::Symbol) = haskey(op.tags, key)
Base.getindex(op::OperatorKind, key::Symbol) = op.tags[key]

function Base.getproperty(op::OperatorKind, name::Symbol)
    name == :tags && return getfield(op, :tags)
    return if haskey(op, name)
        getfield(op, :tags)[name]
    else
        throw(error("type $(typeof(op)) has no field $name"))
    end
end

function get_tags(op::OperatorKind; excluded_keys=nothing)
    if isnothing(excluded_keys)
        return op.tags
    else
        indices_filtered = filter(keys(op.tags)) do key
            !any(isequal.(key, excluded_keys))
        end
        return getindices(op.tags, indices_filtered)
    end
end

# Group operatorkinds in a list based on their equality of tags
function get_operator_groups(
    op_list::AbstractVector{<:OperatorKind}; op_filter=op -> true, excluded_keys=nothing
)
    # Optionally sift through operatorkinds using the filter function
    op_list_filtered = filter(op_filter, op_list)

    # Optionally filter tags to make the grouping on a subset of tags instead
    tags_filtered = get_tags.(op_list_filtered; excluded_keys=excluded_keys)

    # Get connectivity matrix
    C = isdictequal.(reshape(tags_filtered, :, 1), reshape(tags_filtered, 1, :))

    # Find group indices in the matrix
    group_indices = Vector{Int}[]
    for col in eachcol(C)
        push!(group_indices, findall(col))
    end
    group_indices = unique(group_indices)

    return getindex.(Ref(op_list_filtered), group_indices)
end

function statictuple(kind::OperatorKind{K}) where {K}
    pairtypes = Tuple(Pair(Val(pair[1]), Val(pair[2])) for pair in pairs(kind.tags))
    return tuple(Val(K), pairtypes...)
end

const Hamiltonian = OperatorKind{:Hamiltonian}
const Overlap = OperatorKind{:Overlap}

function allows_symmetry(operatorkinds)
    return all(allows_symmetry, operatorkinds)
end

function allows_symmetry(operatorkind::OperatorKind)
    return (:spin ∉ keys(operatorkind.tags)) || isequal(operatorkind.spin, :none)
end

### PARSING ###

function get_operatorkinds end

const SOURCE_TAGS = (:ref, :pred)
const SPIN_TAGS = (:none, :soc, :up, :down)

function sift_tags(all_tags; requested_tags=nothing, removed_tags=nothing)
    @argcheck isnothing(requested_tags) || isnothing(removed_tags)

    if isnothing(requested_tags) && isnothing(removed_tags)
        return all_tags

    elseif !isnothing(requested_tags)
        if !all(requested_tags .∈ Ref(all_tags))
            throw(error("Not all requested tags are available"))
        else
            return Tuple(requested_tags)
        end

    elseif !isnothing(removed_tags)
        return Tuple(tag for tag in all_tags if tag ∉ removed_tags)
    end
end

get_operatorkinds(::Val{:h}) = vcat(
    get_operatorkinds(Val(:hnospin)),
    get_operatorkinds(Val(:hspin)),
)
get_operatorkinds(::Val{:hnospin}) = [
    Hamiltonian(; source=source)
    for source in SOURCE_TAGS
]
get_operatorkinds(::Val{:hspin}) = [
    Hamiltonian(; source=source, spin=spin)
    for source in SOURCE_TAGS
    for spin in SPIN_TAGS
]

get_operatorkinds(::Val{:s}) = vcat(
    get_operatorkinds(Val(:snospin)),
    get_operatorkinds(Val(:sspin)),
)
get_operatorkinds(::Val{:snospin}) = [
    Overlap(; source=source)
    for source in SOURCE_TAGS
]
get_operatorkinds(::Val{:sspin}) = [
    Overlap(; source=source, spin=spin)
    for source in SOURCE_TAGS
    for spin in SPIN_TAGS
]

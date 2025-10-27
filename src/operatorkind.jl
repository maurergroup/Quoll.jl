using Base
using ArgCheck
using Dictionaries
using AutoHashEquals

@auto_hash_equals struct OperatorKind{K}
    tags::Dictionary{Symbol, Symbol}

    function OperatorKind{K}(tags) where K
        @argcheck K isa Symbol
        new{K}(tags)
    end
end

function OperatorKind{K}(pairs::Vararg{<:Pair}) where K
    return OperatorKind{K}(dictionary(pairs))
end

function OperatorKind{K}(; kwargs...) where K
    return OperatorKind{K}(dictionary(kwargs))
end

Base.getindex(op::OperatorKind, key::Symbol) = op.tags[key]

function Base.getproperty(op::OperatorKind, name::Symbol) 
    name == :tags && return getfield(op, :tags)
    haskey(getfield(op, :tags), name) ? op.tags[name] : throw(error("type $(typeof(op)) has no field $name"))
end

const Hamiltonian = OperatorKind{:Hamiltonian}
const Overlap = OperatorKind{:Overlap}

# Technically eveything below is for parsing only
# TODO: I should probably consider modifying the `operators` field in the input file
# to explicitly accept tags so that we don't need to define `get_operatorkinds` for every case.
# Also I'm not sure if keeping this code here is the best idea

function get_operatorkinds end

const SOURCE_TAGS = (:ref, :pred)
const SPIN_TAGS = (:none, :soc, :up, :down)

function sift_tags(all_tags; requested_tags = nothing, removed_tags = nothing)
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

get_operatorkinds(::Val{:h}) = [
    Hamiltonian(source = source, spin = spin)
    for source in SOURCE_TAGS
    for spin in SPIN_TAGS
]
get_operatorkinds(::Val{:href}) = [
    Hamiltonian(source = source, spin = spin)
    for source in sift_tags(SOURCE_TAGS, requested_tags = (:ref,))
    for spin in SPIN_TAGS
]
get_operatorkinds(::Val{:hpred}) = [
    Hamiltonian(source = source, spin = spin)
    for source in sift_tags(SOURCE_TAGS, requested_tags = (:pred,))
    for spin in SPIN_TAGS
]
get_operatorkinds(::Val{:hsoc}) = [
    Hamiltonian(source = source, spin = spin)
    for source in SOURCE_TAGS
    for spin in sift_tags(SPIN_TAGS, requested_tags = (:soc,))
]
get_operatorkinds(::Val{:hpolarised}) = [
    Hamiltonian(source = source, spin = spin)
    for source in SOURCE_TAGS
    for spin in sift_tags(SPIN_TAGS, requested_tags = (:up, :down))
]

get_operatorkinds(::Val{:s}) = [
    Overlap(source = source)
    for source in SOURCE_TAGS
]
get_operatorkinds(::Val{:sref}) = [
    Overlap(source = source)
    for source in sift_tags(SOURCE_TAGS, requested_tags = (:ref,))
]
get_operatorkinds(::Val{:spred}) = [
    Overlap(source = source)
    for source in sift_tags(SOURCE_TAGS, requested_tags = (:pred,))
]

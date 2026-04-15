### SOURCE ###

"""
    AbstractSource

Abstract supertype for format sources. Concrete subtypes (e.g. `CanonicalSource`,
`DeepHSource`, `FHIaimsSource`) carry format-specific configuration such as the default SH
convention and data layout. Used to dispatch format-specific methods throughout the pipeline.
"""
abstract type AbstractSource end

"""
    namedtuple(source::AbstractSource)

Convert a source's fields to a `NamedTuple`. Useful for forwarding source properties as
keyword arguments.
"""
function namedtuple(source::AbstractSource)
    return NamedTuple(key => getproperty(source, key) for key in propertynames(source))
end

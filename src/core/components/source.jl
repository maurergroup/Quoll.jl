### SOURCE ###

abstract type AbstractSource end

function namedtuple(source::AbstractSource)
    return NamedTuple(key => getproperty(source, key) for key in propertynames(source))
end

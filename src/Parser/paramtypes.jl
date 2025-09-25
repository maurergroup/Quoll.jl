using Configurations

@option struct InputParams
    format::String
    directory::Vector{String}
end

@option struct QuollParams
    input::InputParams
end

function Configurations.from_dict(
    ::Type{InputParams},
    ::OptionField{:directory},
    ::Type{Vector{String}},
    s,
)
    return isa(s, String) ? [s] : s
end

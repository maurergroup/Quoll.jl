using Configurations

@option struct InputParams
    format::String
    directory::Vector{String}
end

@option struct QuollParams
    input::InputParams
end

Configurations.from_dict(
    ::Type{InputParams},
    ::OptionField{:directory},
    ::Type{Vector{String}},
    s,
) = isa(s, String) ? [s] : s

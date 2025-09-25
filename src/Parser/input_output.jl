using Configurations
using ArgCheck
using ..Utils

@option struct InputParams
    format::String
    directory::Vector{String}
end

@option struct OutputParams
    format::String
    directory::String
end

function Configurations.from_dict(
    ::Type{InputParams},
    ::OptionField{:directory},
    ::Type{Vector{String}},
    dirs,
)
    # Ensure `dirs` is of type Vector{String}
    @argcheck isa(dirs, String) || isa(dirs, Vector{String})
    dirs = isa(dirs, String) ? [dirs] : dirs

    # Convert `dirs` to absolute paths
    dirs = abspath.(dirs)

    # Walk the directories to obtain all the paths
    return vcat(Utils.find_leafdirs.(dirs)...)
end

function Configurations.from_dict(
    ::Type{OutputParams},
    ::OptionField{:directory},
    ::Type{String},
    dir,
)
    # Convert `dir` to absolute path
    dir = abspath(dir)
    return dir
end

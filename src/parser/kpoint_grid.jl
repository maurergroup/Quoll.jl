@option struct KPointGridParams <: AbstractQuollParams
    grid::Maybe{SVector{3, Int}} = nothing
    density::Float64 = 10.0
    shift::SVector{3, Bool} = SA[false, false, false]
    kpoints::Maybe{Vector{SVector{4, Float64}}} = nothing

    function KPointGridParams(grid, density, shift, kpoints)

        # Check if supplied arguments are appropriate
        !isnothing(grid) && @argcheck all(grid .> 0)
        @argcheck density > 0

        new(grid, density, shift, kpoints)
    end
end

function Configurations.from_dict(
    ::Type{KPointGridParams},
    ::OptionField{:shift},
    ::Type{SVector{3, Bool}},
    shift,
)
    @argcheck isa(shift, Bool) || isa(shift, Vector)

    return isa(shift, Bool) ? SVector{3, Bool}([shift for _ in 1:3]) : SVector{3, Bool}(shift)
end

function Configurations.from_dict(
    ::Type{KPointGridParams},
    ::OptionField{:kpoints},
    ::Type{Vector{SVector{4, Float64}}},
    kpoints,
)
    @argcheck isa(kpoints, String) || isa(kpoints, Vector)
    if kpoints isa String
        @argcheck ispath(kpoints)
        return eachrow(readdlm(kpoints))

    elseif kpoints isa Vector
        return [SVector{4, Float64}(kpoint) for kpoint in kpoints]
    end
end


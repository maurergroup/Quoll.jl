using Configurations
using ArgCheck

@option struct KPointGridParams <: AbstractQuollParams
    grid::Union{SVector{3, Int}, Nothing} = nothing
    density::Float64 = 10.0

    function KPointGridParams(grid, density)
        grid !== nothing && @argcheck all(grid .> 0)
        @argcheck density > 0
        new(grid, density)
    end
end

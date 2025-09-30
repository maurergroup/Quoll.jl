using Configurations
using ArgCheck

@option struct KPointGridParams <: AbstractQuollParams
    grid::Union{Vector{Int}, Nothing} = nothing
    density::Union{Vector{Float64}, Nothing} = nothing

    function KPointGridParams(grid, density)
        @argcheck !(grid === nothing && density === nothing) "Either k-point grid or k-point density must be supplied"
        @argcheck !(grid !== nothing && density !== nothing) "Both k-point grid and k-point density cannot be supplied"
        @argcheck grid .> 0
        @argcheck density .> 0
        new(grid, density)
    end
end

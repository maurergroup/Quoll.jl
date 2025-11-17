@option struct SymmetryParams
    time_reversal::Maybe{Bool} = nothing
    crystal_symmetry::Maybe{Bool} = nothing
    space_group::Maybe{Int} = nothing
    
    function SymmetryParams(time_reversal, crystal_symmetry, space_group)

        # Check whether space group is appropriate
        !isnothing(space_group) && @argcheck space_group ∈ 1:230

        new(time_reversal, crystal_symmetry, space_group)
    end
end
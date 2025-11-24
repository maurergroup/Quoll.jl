@option struct SymmetryParams
    time_reversal::Maybe{Bool} = nothing
    crystal_symmetry::Maybe{Bool} = nothing
    space_group::Maybe{Int} = nothing
    symprec::Float64 = 1e-5
    
    function SymmetryParams(time_reversal, crystal_symmetry, space_group, symprec)

        # Check whether space group is appropriate
        !isnothing(space_group) && @argcheck space_group ∈ 1:230

        new(time_reversal, crystal_symmetry, space_group, symprec)
    end
end
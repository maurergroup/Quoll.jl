function parse_inputfile(input_filepath::AbstractString)
    params = from_toml(QuollParams, input_filepath)

    input_dirs_eachroot = find_leafdirs.(params.input.directory)
    input_dirs = vcat(input_dirs_eachroot...)

    output_dirs = get_output_dirs(
        params.output.directory,
        params.input.directory,
        input_dirs_eachroot,
    )

    return params, input_dirs, output_dirs
end

### VALIDATION METHODS ###

function search_clashes(basis_projection, error_metrics)
    @argcheck !(!isnothing(basis_projection) && error_metrics.mae)
    @argcheck !(!isnothing(basis_projection) && error_metrics.eigenvalue_error)
    @argcheck !(!isnothing(basis_projection) && error_metrics.el_entropy_error)
end

# Check if the operators are sufficient for the requested job
function validate_operatorkinds(
    operatorkinds, output, basis_projection, postprocessing, error_metrics
)

    # Collect operator kinds into groups containing Sref, Spred, Href, Hpred
    # with the rest of the keys matching in a group
    op_filter = op -> op isa Overlap || op isa Hamiltonian
    groups = get_operator_groups(
        operatorkinds; op_filter=op_filter, excluded_keys=[:source]
    )

    for group in groups
        Sref = !isnothing(findfirst(op -> op isa Overlap && op.source == :ref, group))
        Spred = !isnothing(findfirst(op -> op isa Overlap && op.source == :pred, group))
        Href = !isnothing(findfirst(op -> op isa Hamiltonian && op.source == :ref, group))
        Hpred = !isnothing(findfirst(op -> op isa Hamiltonian && op.source == :pred, group))

        #! format: off
        !isnothing(output)           && (@argcheck Sref || Spred || Href || Hpred)
        !isnothing(basis_projection) && (@argcheck Sref)

        postprocessing.dos            && (@argcheck (Sref && Href) || (Spred && Hpred))
        postprocessing.fermi_level    && (@argcheck (Sref && Href) || (Spred && Hpred))
        postprocessing.band_structure && (@argcheck (Sref && Href) || (Spred && Hpred))

        error_metrics.mae              && (@argcheck Href && Hpred)
        error_metrics.eigenvalue_error && (@argcheck (Sref && Href) || (Spred && Hpred))
        error_metrics.el_entropy_error && (@argcheck (Sref && Href) || (Spred && Hpred))
        #! format: on
    end

    return true
end

# Find operatorkinds and check if they are compatible with the supplied parameters
function find_operatorkinds(dir::AbstractString, params::QuollParams)
    found_operatorkinds = find_operatorkinds(dir, params.input.format)
    if !all(params.input.operators .∈ Ref(found_operatorkinds))
        operatorkinds = params.input.operators ∩ found_operatorkinds

        # Validate whether the calculation can still continue
        # even though not all operatorkinds were found
        validate_operatorkinds(
            operatorkinds,
            params.output,
            params.basis_projection,
            params.postprocessing,
            params.error_metrics,
        )
    else
        operatorkinds = params.input.operators
    end

    # Alternatively could use `@warn` if the function enters the if block above
    # but currently it would be entered in almost all cases due to conflicting
    # default values of operatorkinds
    @debug "Intersection between found and requested operators:" operatorkinds
    return operatorkinds
end

function requires_kpoint_grid(params)
    requires = false
    requires |= !isnothing(params.basis_projection)
    requires |= params.postprocessing.fermi_level
    requires |= params.postprocessing.dos
    requires |= params.error_metrics.eigenvalue_error
    requires |= params.error_metrics.el_entropy_error
    return requires
end

function validate_symmetry(input, basis_projection, symmetry)
    validate_time_reversal(input, symmetry)
    validate_crystal_symmetry(input, basis_projection, symmetry)
    return true
end

function validate_time_reversal(input, symmetry)
    allowed = allows_time_reversal(input.operators)
    supplied = symmetry.time_reversal

    @argcheck !(isequal(supplied, true) && isequal(allowed, false))
    return true
end

function validate_crystal_symmetry(input, basis_projection, symmetry)
    allowed = allows_crystal_symmetry(input.operators, basis_projection)
    supplied = symmetry.crystal_symmetry

    @argcheck !(isequal(supplied, true) && isequal(allowed, false))
    return true
end

allows_time_reversal(operatorkinds, ::QuollParams) = allows_time_reversal(operatorkinds)

function allows_time_reversal(operatorkinds)
    allows = true
    allows &= allows_symmetry(operatorkinds)
    return allows
end

allows_crystal_symmetry(operatorkinds, params::QuollParams) =
    allows_crystal_symmetry(operatorkinds, params.basis_projection)

function allows_crystal_symmetry(operatorkinds, basis_projection)
    allows = true
    allows &= isnothing(basis_projection)
    allows &= allows_symmetry(operatorkinds)
    return allows
end

### CONVENIENCE METHODS ###

function get_time_reversal(operatorkinds, params::QuollParams)
    if !isnothing(params.symmetry.time_reversal)
        return params.symmetry.time_reversal
    else
        return allows_time_reversal(operatorkinds, params)
    end
end

function get_crystal_symmetry(operatorkinds, params::QuollParams)
    if !isnothing(params.symmetry.crystal_symmetry)
        return params.symmetry.crystal_symmetry
    else
        return allows_crystal_symmetry(operatorkinds, params)
    end
end

function construct_kgrid(atoms::AbstractSystem, operatorkinds, params::QuollParams)
    return construct_kgrid(
        atoms;
        density=params.kpoint_grid.density,
        mesh=params.kpoint_grid.mesh,
        shift=params.kpoint_grid.shift,
        time_reversal=get_time_reversal(operatorkinds, params),
        crystal_symmetry=get_crystal_symmetry(operatorkinds, params),
        symprec=params.symmetry.symprec,
    )
end

function perform_core_projection(
    operators, kgrid::KGrid, my_ikpoints, comm::MPI.Comm, params::QuollParams
)
    # TODO: in the future if the projection method is not just a singleton
    # but also depends on parameters, those parameters can be parsed in as well
    # and fed into the constructor `method()` here
    return perform_core_projection(
        operators,
        params.basis_projection.projected_basis,
        kgrid,
        my_ikpoints,
        comm;
        method=params.basis_projection.method(),
    )
end

### FILE HANDLING METHODS ###

# Assuming all supplied paths are absolute
function get_output_dirs(output_root, input_roots, input_eachroot)
    # Including basename if multiple roots were supplied
    # to keep individual trees separate
    include_input_root_basename = length(input_eachroot) > 1

    output_dirs = String[]
    for (abs_input_root, abs_input_leaves) in zip(input_roots, input_eachroot)
        rel_leaves = relpath.(abs_input_leaves, abs_input_root)

        if include_input_root_basename
            input_root_basename = basename(abs_input_root)
            rel_leaves = joinpath.(input_root_basename, rel_leaves)
        end

        push!(output_dirs, normpath.(joinpath.(output_root, rel_leaves))...)
    end

    return output_dirs
end

"""
    find_leafdirs(rootdir::T) where T

Find leaf directories (directories that do not contain
directories themselves) going down from `rootdir` directory.
Absolute paths are returned.
"""
function find_leafdirs(rootdir::T) where {T}
    leafdirs = T[]
    for (root, dirs, _) in walkdir(rootdir)
        if isempty(dirs)
            push!(leafdirs, normpath(joinpath(rootdir, root)))
        end
    end
    return leafdirs
end

# Basis projection
# - requires k-point grid
# - incompatible with error metrics

# Postprocessing
# - most of them require k-point grid (except band structure)
# - should store eigenvalues, fermi level and pass them to error metrics if required

# TODO: DOS requires fermi level calculation
# What happens in dos == true but fermi_level == false?
# In that case I could just ignore the fermi level parameter

# Error metrics
# - some require k-point grid and eigenvalues (e.g. eigenvalue error)
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
function validate_operatorkinds(operatorkinds, output, basis_projection, postprocessing, error_metrics)
    Sref = Overlap(source = :ref) ∈ operatorkinds
    Spred = Overlap(source = :pred) ∈ operatorkinds

    spins = [h.spin for h in operatorkinds if h isa Hamiltonian]

    for spin in spins
        Href = Hamiltonian(source = :ref, spin = spin) ∈ operatorkinds
        Hpred = Hamiltonian(source = :pred, spin = spin) ∈ operatorkinds

        output !== nothing && (@argcheck Sref || Spred || Href || Hpred)
        basis_projection !== nothing && (@argcheck Sref)

        postprocessing.dos && (@argcheck (Sref && Href) || (Spred && Hpred))
        postprocessing.fermi_level && (@argcheck (Sref && Href) || (Spred && Hpred))
        postprocessing.band_structure && (@argcheck (Sref && Href) || (Spred && Hpred))

        error_metrics.mae && (@argcheck Href && Hpred)
        error_metrics.eigenvalue_error && (@argcheck (Sref && Href) || (Spred && Hpred))
        error_metrics.el_entropy_error && (@argcheck (Sref && Href) || (Spred && Hpred))
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

        push!(output_dirs, joinpath.(output_root, rel_leaves)...)
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
            push!(leafdirs, joinpath(rootdir, root))
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
# Check if the operators are sufficient for requested job
function validate_operatorkinds(operatorkinds, output, basis_projection, postprocessing, error_metrics)
    Sref = Overlap(:ref) ∈ operatorkinds
    Spred = Overlap(:pred) ∈ operatorkinds

    spins = [h.spin for h in operatorkinds if h isa Hamiltonian]

    for spin in spins
        Href = Hamiltonian(:ref, spin) ∈ operatorkinds
        Hpred = Hamiltonian(:pred, spin) ∈ operatorkinds

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

# [h.spin for h in operatorkinds if h isa Hamiltonian]

# Find operatorkinds and check if they are compatible with the supplied parameters
function find_operatorkinds(dir::AbstractString, params::QuollParams)
    found_operatorkinds = find_operatorkinds(dir, params.input.format)
    if !all(params.input.operators .∈ Ref(found_operatorkinds))
        operatorkinds = collect(intersect(Set(params.input.operators), Set(found_operatorkinds)))

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
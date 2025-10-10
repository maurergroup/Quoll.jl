# Check if the operators are sufficient for requested job
function validate_operatorkinds(operatorkinds, output, basis_projection, postprocessing, error_metrics)
    Sref = Overlap(:ref) ∈ operatorkinds
    Spred = Overlap(:pred) ∈ operatorkinds
    Href = Hamiltonian(:ref) ∈ operatorkinds
    Hpred = Hamiltonian(:pred) ∈ operatorkinds

    output !== nothing && (@argcheck Sref || Spred || Href || Hpred)
    basis_projection !== nothing && (@argcheck Sref)

    postprocessing.dos && (@argcheck (Sref && Href) || (Spred && Hpred))
    postprocessing.fermi_level && (@argcheck (Sref && Href) || (Spred && Hpred))
    postprocessing.band_structure && (@argcheck (Sref && Href) || (Spred && Hpred))

    error_metrics.mae && (@argcheck Href && Hpred)
    error_metrics.eigenvalue_error && (@argcheck (Sref && Href) || (Spred && Hpred))
    error_metrics.el_entropy_error && (@argcheck (Sref && Href) || (Spred && Hpred))

    return true
end

# Find operatorkinds and check if they are compatible with the supplied parameters
function find_operatorkinds(dir::AbstractString, params::QuollParams)
    found_operatorkinds = find_operatorkinds(dir, params.input.format)
    if !all(params.input.operators .∈ Ref(found_operatorkinds))
        @warn "Could not find all the requested operators"
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
    return operatorkinds
end
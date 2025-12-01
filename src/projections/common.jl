# TODO: Would be quite different for a non-periodic system (no need for a k-point grid)
function perform_core_projection(operators, projected_basis, kgrid::KGrid, my_ikpoints, comm::MPI.Comm,
    method::AbstractBasisProjection; recip_format::Type{O} = DenseOperator) where O<:AbstractOperator

    my_kpoints = view(kgrid.kpoints, my_ikpoints)
    my_weights = view(kgrid.weights, my_ikpoints)
    my_kpoints_float = convert(Vector{SVector{3, Float64}}, my_kpoints)

    # Precompute phases
    # - Needs to be done separately for each operator because in theory the sparsity can differ
    phases_list = [
        collect(eachcol(precompute_phases(my_kpoints_float, get_sparsity(operator).images)))
        for operator in operators
    ]

    # Construct core and valence masks
    # - Again we don't assume that each operator can use the same mask
    c_mask_list = [
        get_dense_subbasis_mask(get_basisset(operator), projected_basis, inverted = false)
        for operator in operators
    ]
    v_mask_list = [
        get_dense_subbasis_mask(get_basisset(operator), projected_basis, inverted = true)
        for operator in operators
    ]

    # Create empty valence-only operators which will contain local contributions
    # and will be synchronised across MPI tasks at the end
    # (for now keeping the same sparsity as in the input operators)
    v_operators = [
        build_operator(
            typeof(operator),
            get_kind(operator),
            construct_subbasis_metadata(
                get_metadata(operator),
                projected_basis,
                inverted = true,
            ),
            uninit = false,
            value = zero(get_float(operator))
        )
        for operator in operators
    ]

    for (ik, (kpoint, weight)) in enumerate(zip(my_kpoints, my_weights))

        # Perform fourier transform, mapping the input operators to the k-point
        recip_operators = [
            fourier_transform(
                operator, tuple(kpoint), tuple(phases[ik]), O,
                float = get_float(operator),
            )
            for (operator, phases) in zip(operators, phases_list)
        ]

        # Obtain projected valence-only operator data at the k-point
        recip_v_operator_data_list = compute_valence_operator_data(
            recip_operators,
            c_mask_list,
            v_mask_list,
            method,
        )
        
        # Construct reciprocal valence-only operators at the k-point using the computed data
        recip_v_metadata_list = [
            construct_subbasis_metadata(
                get_metadata(recip_operator),
                get_basisset(v_operator),
                get_spins(v_operator),
            )
            for (recip_operator, v_operator) in zip(recip_operators, v_operators)
        ]
        recip_v_operators = [
            build_operator(
                typeof(operator), get_kind(operator), data, metadata,
            )
            for (operator, data, metadata)
            in zip(recip_operators, recip_v_operator_data_list, recip_v_metadata_list)
        ]

        # Back-transform the reciprocal valence-only operators to real space
        for (v_operator, recip_v_operator, phases) in zip(v_operators, recip_v_operators, phases_list)
            inv_fourier_transform_data!(v_operator, recip_v_operator, phases[ik], weight)
        end

    end

    # Synchronise valence-only real-space operators across MPI tasks that worked on different k-points
    for v_operator in v_operators
        synchronise_data!(v_operator, comm)
    end

    return v_operators
end

function core_valence_partition(operator::DenseOperator, core_mask::BitVector, valence_mask::BitVector)
    O = get_data(operator)
    O₁₁ = view(O, core_mask, core_mask)
    O₁₂ = view(O, core_mask, valence_mask)
    O₂₂ = view(O, valence_mask, valence_mask)
    return O₁₁, O₁₂, O₂₂
end

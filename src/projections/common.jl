# TODO: Would be quite different for a non-periodic system (no need for a k-point grid)
# TODO: Need to check whether symmetries of KGrid are appropriate
function perform_core_projection(
    operators, projected_basis, kgrid::KGrid, my_ikpoints, comm::MPI.Comm;
    method::AbstractBasisProjection=LaikovCore(),
    recip_operator_type::Type{OPₖ}=Operator,
    recip_format::Type{Mₖ}=CanonicalDenseNoSpinRecipMetadata,
) where {OPₖ<:AbstractOperator,Mₖ<:AbstractMetadata}
    @argcheck !(kgrid.crystal_symmetry)

    my_kpoints = view(kgrid.kpoints, my_ikpoints)
    my_weights = view(kgrid.weights, my_ikpoints)
    my_kpoints_float = convert(Vector{SVector{3,Float64}}, my_kpoints)

    # Precompute phases
    # - Needs to be done separately for each operator because in theory the sparsity can differ
    phases_list = [
        collect(
            eachcol(precompute_phases(my_kpoints_float, op_images(op_sparsity(operator))))
        )
        for operator in operators
    ]

    # Construct core and valence masks
    # - Again we don't assume that each operator can use the same mask
    c_dense_mask_list = [
        get_dense_subbasis_mask(op_basisset(operator), projected_basis; inverted=false)
        for operator in operators
    ]
    v_dense_mask_list = [
        get_dense_subbasis_mask(op_basisset(operator), projected_basis; inverted=true)
        for operator in operators
    ]

    # Create empty valence-only operators which will contain local contributions
    # and will be synchronised across MPI tasks at the end
    # (keeping the same sparsity as in the input operators)
    v_operators = [
        build_operator(
            typeof(operator),
            op_metadata(operator);
            subbasis=projected_basis,
            inverted=true,
            initialised=true,
            value=zero(op_float_type(operator)),
        )
        for operator in operators
    ]

    for (ik, (kpoint, weight)) in enumerate(zip(my_kpoints, my_weights))

        # Perform fourier transform, mapping the input operators to the k-point
        recip_operators = [
            fourier_transform(
                OPₖ, Mₖ, operator, kpoint, phases[ik]; out_shconv=op_shconv(operator)
            )
            for (operator, phases) in zip(operators, phases_list)
        ]

        # Obtain projected valence-only operator data at the k-point
        recip_v_data_list = compute_valence_operator_data(
            recip_operators,
            c_dense_mask_list,
            v_dense_mask_list,
            method,
        )

        # Construct reciprocal valence-only operators at the k-point using the computed data
        recip_v_metadata_list = [
            convert_metadata(
                typeof(op_metadata(recip_operator)),
                op_metadata(recip_operator);
                subbasis=projected_basis,
                inverted=true,
            )
            for recip_operator in recip_operators
        ]
        recip_v_operators = [
            build_operator(OPₖ, recip_v_metadata, recip_v_data)
            for (recip_v_metadata, recip_v_data) in
            zip(recip_v_metadata_list, recip_v_data_list)
        ]

        # Back-transform the reciprocal valence-only operators to real space
        for (v_operator, recip_v_operator, phases) in
            zip(v_operators, recip_v_operators, phases_list)
            inv_fourier_transform_data!(v_operator, recip_v_operator, phases[ik], weight)
        end
    end

    # Synchronise valence-only real-space operators across MPI tasks
    # that worked on different k-points
    for v_operator in v_operators
        synchronise_data!(v_operator, comm)
    end

    return v_operators
end

function core_valence_partition(operator::AbstractOperator, core_mask, valence_mask)
    data = op_data(operator)
    return core_valence_partition(data, core_mask, valence_mask)
end

function core_valence_partition(
    data::DenseRecipData, core_mask::BitVector, valence_mask::BitVector
)
    O = unwrap_data(data)
    O₁₁ = view(O, core_mask, core_mask)
    O₁₂ = view(O, core_mask, valence_mask)
    O₂₂ = view(O, valence_mask, valence_mask)
    return O₁₁, O₁₂, O₂₂
end
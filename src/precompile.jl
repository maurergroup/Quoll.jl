# Precompile workloads exercising the most common code paths from the regression tests.
# The goal is to drive compilation of the core pipeline components (load, convert, Fourier
# transform, write) for the metadata/operator variants used across regression tests so that
# the `quoll` app and typical library entry points incur minimal latency on first use.
#
# Paths covered (mirroring the regression tests):
#   - FHIaimsCSC (NoSpin) -> CanonicalBlock -> DeepHBlock (with and without radii)
#   - DeepHBlock (NoSpin) round-trip: write -> read -> convert
#   - DeepHBlock (Spin) -> CanonicalBlock (Spin) -> CanonicalDenseRecip (Spin) via fourier_transform
#   - CanonicalBlock -> CanonicalDenseRecip via fourier_transform (NoSpin)
#   - Core projection primitives: LaikovCore and FC99V compute_valence_data,
#     inv_fourier_transform_data!, subbasis mask/projection machinery.
#   - TOML parsing to QuollParams.
#
# Concrete metadata variants are instantiated in-memory using the same shape as
# `test/TestFixtures.jl`, so no artifact data is needed at precompile time.

using PrecompileTools
using AtomsBase
using StaticArrays
using Unitful
using Dictionaries

@setup_workload begin
    _pc_H  = ChemicalSpecies(:H)
    _pc_Li = ChemicalSpecies(:Li)

    _pc_cellvecs = SA[
        SA[5.0, 0.0, 0.0],
        SA[0.0, 5.0, 0.0],
        SA[0.0, 0.0, 5.0],
    ]u"Å"

    _pc_atoms = periodic_system(
        [
            _pc_H  => SA[0.0, 0.0, 0.0]u"Å",
            _pc_Li => SA[2.5, 2.5, 2.5]u"Å",
        ],
        _pc_cellvecs,
    )

    # Wiki-ordered basis (1 s-orbital + 1 p-shell per species == 4 functions each)
    _pc_basis_dict = Base.ImmutableDict(
        _pc_H => [
            BasisMetadata(_pc_H, 1, 0,  0, nothing),
            BasisMetadata(_pc_H, 2, 1, -1, nothing),
            BasisMetadata(_pc_H, 2, 1,  0, nothing),
            BasisMetadata(_pc_H, 2, 1,  1, nothing),
        ],
        _pc_Li => [
            BasisMetadata(_pc_Li, 2, 0,  0, nothing),
            BasisMetadata(_pc_Li, 2, 1, -1, nothing),
            BasisMetadata(_pc_Li, 2, 1,  0, nothing),
            BasisMetadata(_pc_Li, 2, 1,  1, nothing),
        ],
    )
    _pc_atom2species = species(_pc_atoms, :)
    _pc_basisset = BasisSetMetadata(_pc_basis_dict, _pc_atom2species)

    # Non-hermitian block sparsity (like DeepH storage): full 2x2 atom pairs
    _pc_block_sparsity_nonherm = BlockRealSparsity(
        dictionary([
            (1, 1) => [SA[0, 0, 0]],
            (1, 2) => [SA[0, 0, 0]],
            (2, 1) => [SA[0, 0, 0]],
            (2, 2) => [SA[0, 0, 0]],
        ]),
        [SA[0, 0, 0]],
        false,
    )

    # Hermitian block sparsity (upper + diagonal). Used for core projection fixtures.
    _pc_block_sparsity_herm = BlockRealSparsity(
        dictionary([
            (1, 1) => [SA[0, 0, 0]],
            (1, 2) => [SA[0, 0, 0]],
            (2, 2) => [SA[0, 0, 0]],
        ]),
        [SA[0, 0, 0]],
        true,
    )

    # Small hermitian CSC sparsity matching the basis (one entry per column = diagonal only)
    _pc_nbasis = 4 + 4
    _pc_csc_rowval     = collect(1:_pc_nbasis)
    _pc_csc_colcellptr = Array{Int,3}(undef, 2, 1, _pc_nbasis)
    for _ib in 1:_pc_nbasis
        _pc_csc_colcellptr[1, 1, _ib] = _ib
        _pc_csc_colcellptr[2, 1, _ib] = _ib
    end
    _pc_csc_sparsity = CSCRealSparsity(
        _pc_csc_rowval, _pc_csc_colcellptr, [SA[0, 0, 0]], true,
    )

    _pc_kind_h     = Hamiltonian(; source=:ref)
    _pc_kind_s     = Overlap(; source=:ref)
    _pc_kind_h_soc = Hamiltonian(; source=:ref, spin=:soc)

    _pc_canon_source = CanonicalSource()
    _pc_deeph_source = DeepHSource()
    _pc_fhi_source   = FHIaimsSource()

    _pc_canon_shconv = default_shconv(_pc_canon_source)
    _pc_deeph_shconv = default_shconv(_pc_deeph_source)
    _pc_fhi_shconv   = default_shconv(_pc_fhi_source)

    _pc_deeph_basisset = convert_basisset_shconv(_pc_basisset, _pc_deeph_shconv)
    _pc_fhi_basisset   = convert_basisset_shconv(_pc_basisset, _pc_fhi_shconv)

    # SOC-doubled DeepH basis (DeepH stores both spins concatenated when soc)
    _pc_deeph_soc_basis_dict = Base.ImmutableDict(
        (
            z => vcat(_pc_deeph_basisset.basis[z], _pc_deeph_basisset.basis[z])
            for z in keys(_pc_deeph_basisset.basis)
        )...,
    )
    _pc_deeph_soc_basisset = BasisSetMetadata(
        _pc_deeph_soc_basis_dict, _pc_atom2species,
    )

    _pc_radii = Dict{ChemicalSpecies,LengthAngstrom}(
        _pc_H  => 3.0u"Å",
        _pc_Li => 3.0u"Å",
    )

    _pc_core_basis = [
        BasisMetadata(_pc_H,  1, 0, 0, nothing),
        BasisMetadata(_pc_Li, 2, 0, 0, nothing),
    ]

    @compile_workload begin

        ### ---------- CANONICAL BLOCK REAL (NoSpin) ---------- ###

        _pc_canon_basic = BasicMetadataContainer(
            _pc_kind_h, _pc_canon_source, _pc_block_sparsity_nonherm,
            _pc_basisset, _pc_canon_shconv, _pc_atoms,
        )
        _pc_canon_meta = RealMetadata(_pc_canon_basic)

        _pc_canon_op_plain = build_operator(Operator,      _pc_canon_meta)
        _pc_canon_op_keyed = build_operator(KeyedOperator, _pc_canon_meta)

        ### ---------- FOURIER TRANSFORM (NoSpin) ---------- ###

        _pc_kpt = SA[0.0, 0.0, 0.0]
        fourier_transform(
            Operator, CanonicalDenseRecipMetadata, _pc_canon_op_keyed, _pc_kpt,
        )
        fourier_transform(
            KeyedOperator, CanonicalDenseRecipMetadata, _pc_canon_op_keyed, _pc_kpt,
        )

        ### ---------- DEEPH BLOCK REAL (NoSpin) ---------- ###

        _pc_deeph_basic = BasicMetadataContainer(
            _pc_kind_h, _pc_deeph_source, _pc_block_sparsity_nonherm,
            _pc_deeph_basisset, _pc_deeph_shconv, _pc_atoms,
        )
        _pc_deeph_meta = RealMetadata(_pc_deeph_basic)
        _pc_deeph_op   = build_operator(Operator, _pc_deeph_meta)

        # DeepH -> Canonical
        _pc_canon_from_deeph = convert_operator(
            KeyedOperator, CanonicalBlockRealMetadata, _pc_deeph_op,
        )

        # Canonical -> DeepH (hermitian input, non-hermitian output — the silicon carbide path)
        _pc_canon_basic_herm = BasicMetadataContainer(
            _pc_kind_h, _pc_canon_source, _pc_block_sparsity_herm,
            _pc_basisset, _pc_canon_shconv, _pc_atoms,
        )
        _pc_canon_herm_keyed = build_operator(
            KeyedOperator, RealMetadata(_pc_canon_basic_herm),
        )
        _pc_deeph_from_canon = convert_operator(
            Operator, DeepHBlockRealMetadata, _pc_canon_herm_keyed;
            hermitian=false,
        )

        ### ---------- DEEPH WRITE / READ ROUND-TRIP ---------- ###

        mktempdir() do _pc_dir
            write_operators(
                DeepHBlockNoSpinRealMetadata, _pc_dir, [_pc_deeph_from_canon],
            )
            load_operator(
                Operator, DeepHBlockNoSpinRealMetadata, _pc_dir,
                Hamiltonian(; source=:ref),
            )
            find_operatorkinds(DeepHBlockNoSpinRealMetadata, _pc_dir)
        end

        ### ---------- FHI-AIMS CSC -> CANONICAL ---------- ###

        _pc_fhi_basic = BasicMetadataContainer(
            _pc_kind_h, _pc_fhi_source, _pc_csc_sparsity,
            _pc_fhi_basisset, _pc_fhi_shconv, _pc_atoms,
        )
        _pc_fhi_meta = RealMetadata(_pc_fhi_basic)
        _pc_fhi_data = wrap_data(typeof(_pc_fhi_meta), zeros(Float64, _pc_nbasis))
        _pc_fhi_op   = Operator(_pc_fhi_meta, _pc_fhi_data)

        # Path without radii (sparsity inherited from FHI-aims CSC indices)
        convert_operator(
            KeyedOperator, CanonicalBlockRealMetadata, _pc_fhi_op,
        )

        # Path with radii (sparsity rebuilt from neighbour list — the Au / H2O regression tests)
        convert_operator(
            KeyedOperator, CanonicalBlockRealMetadata, _pc_fhi_op;
            radii=_pc_radii,
        )

        ### ---------- SPIN PATH: DEEPH SOC -> CANONICAL -> DENSE RECIP ---------- ###

        _pc_deeph_soc_spins = SpinsMetadata(
            _pc_deeph_source, _pc_kind_h_soc, _pc_deeph_soc_basisset,
        )
        _pc_deeph_soc_basic = BasicMetadataContainer(
            _pc_kind_h_soc, _pc_deeph_source, _pc_block_sparsity_nonherm,
            _pc_deeph_soc_basisset, _pc_deeph_shconv, _pc_atoms,
        )
        _pc_deeph_spin_meta = SpinRealMetadata(_pc_deeph_soc_basic, _pc_deeph_soc_spins)
        _pc_deeph_spin_op = build_operator(
            Operator, _pc_deeph_spin_meta; value=zero(ComplexF64),
        )
        _pc_canon_spin = convert_operator(
            KeyedOperator, CanonicalBlockRealMetadata, _pc_deeph_spin_op,
        )
        fourier_transform(
            KeyedOperator, CanonicalDenseRecipMetadata, _pc_canon_spin,
            SA[0.0, 0.5, 0.5],
        )

        ### ---------- CORE PROJECTION PRIMITIVES ---------- ###

        # Build hermitian real-space canonical H/S to feed the projection machinery
        _pc_canon_basic_s_herm = BasicMetadataContainer(
            _pc_kind_s, _pc_canon_source, _pc_block_sparsity_herm,
            _pc_basisset, _pc_canon_shconv, _pc_atoms,
        )
        _pc_canon_h_real = build_operator(
            KeyedOperator, RealMetadata(_pc_canon_basic_herm),
        )
        _pc_canon_s_real = build_operator(
            KeyedOperator, RealMetadata(_pc_canon_basic_s_herm),
        )

        # Exercise subbasis machinery (drives convert_metadata/convert_spins subbasis paths)
        get_dense_subbasis_mask(_pc_basisset, _pc_core_basis; inverted=false)
        _pc_v_mask = get_dense_subbasis_mask(_pc_basisset, _pc_core_basis; inverted=true)
        _pc_c_mask = get_dense_subbasis_mask(_pc_basisset, _pc_core_basis; inverted=false)

        # Valence-only metadata/operator construction (subbasis + inverted)
        build_operator(
            KeyedOperator, RealMetadata(_pc_canon_basic_herm);
            subbasis=_pc_core_basis, inverted=true,
        )

        # Fourier-transform to dense recip (the intermediate used by perform_core_projection)
        _pc_recip_h = fourier_transform(
            Operator, CanonicalDenseNoSpinRecipMetadata, _pc_canon_h_real, _pc_kpt,
        )[1]
        _pc_recip_s = fourier_transform(
            Operator, CanonicalDenseNoSpinRecipMetadata, _pc_canon_s_real, _pc_kpt,
        )[1]

        # Seed the overlap with an identity block so LaikovCore's `inv(S₁₁)` succeeds.
        _pc_s_body = unwrap_data(op_data(_pc_recip_s))
        for _i in axes(_pc_s_body, 1)
            _pc_s_body[_i, _i] = one(ComplexF64)
        end

        Projections.compute_valence_data(
            [_pc_recip_h, _pc_recip_s],
            [_pc_c_mask, _pc_c_mask],
            [_pc_v_mask, _pc_v_mask],
            LaikovCore(),
        )
        Projections.compute_valence_data(
            [_pc_recip_h, _pc_recip_s],
            [_pc_c_mask, _pc_c_mask],
            [_pc_v_mask, _pc_v_mask],
            FC99V(),
        )

        # Exercise inverse Fourier transform (real-space accumulation loop)
        _pc_v_meta = convert_metadata(
            typeof(op_metadata(_pc_canon_h_real)),
            op_metadata(_pc_canon_h_real);
            subbasis=_pc_core_basis, inverted=true,
        )
        _pc_v_h_real = build_operator(
            KeyedOperator, _pc_v_meta; value=zero(Float64),
        )
        _pc_phases = precompute_phases([_pc_kpt], op_images(op_sparsity(_pc_canon_h_real)))
        _pc_recip_v_meta = convert_metadata(
            typeof(op_metadata(_pc_recip_h)),
            op_metadata(_pc_recip_h);
            subbasis=_pc_core_basis, inverted=true,
        )
        _pc_recip_v = build_operator(
            Operator, _pc_recip_v_meta; value=zero(ComplexF64),
        )
        inv_fourier_transform_data!(
            _pc_v_h_real, _pc_recip_v, _pc_phases[:, 1], 1.0,
        )

        ### ---------- TOML PARSING ---------- ###

        mktempdir() do _pc_dir
            _pc_input_dir = joinpath(_pc_dir, "input")
            mkpath(_pc_input_dir)
            _pc_toml_path = joinpath(_pc_dir, "quoll.toml")
            open(_pc_toml_path, "w") do _io
                write(_io,
                    """
                    [input]
                    format = "DeepH"
                    directory = "$(_pc_input_dir)"
                    operators = ["H", "S"]

                    [output]
                    format = "DeepH"
                    directory = "$(joinpath(_pc_dir, "out"))"
                    hermitian = false
                    """,
                )
            end
            Parser.parse_inputfile(_pc_toml_path)
        end
    end
end

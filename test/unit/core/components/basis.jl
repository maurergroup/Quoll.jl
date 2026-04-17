using Quoll
using Dictionaries
using AtomsBase

using Test

const basisset = let a = ChemicalSpecies(:H1), b = ChemicalSpecies(:H2)
    Quoll.BasisSetMetadata(
        dictionary([
            a => [
                Quoll.BasisMetadata(a, 1, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, -1, nothing)
                Quoll.BasisMetadata(a, 2, 1, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, 1, nothing)
            ],
            b => [
                Quoll.BasisMetadata(b, 2, 0, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, -1, nothing)
                Quoll.BasisMetadata(b, 2, 1, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, 1, nothing)
            ],
        ]),
        [a, b, b, a],
    )
end

const basisset_degen = let a = ChemicalSpecies(:H1), b = ChemicalSpecies(:H2)
    Quoll.BasisSetMetadata(
        Base.ImmutableDict(
            a => [
                Quoll.BasisMetadata(a, 1, 0, 0,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(a, 2, 0, 0,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(a, 2, 1, -1, Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(a, 2, 1, 0,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(a, 2, 1, 1,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(a, 2, 1, -1, Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(a, 2, 1, 0,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(a, 2, 1, 1,  Base.ImmutableDict(:type => :atomic))
            ],
            b => [
                Quoll.BasisMetadata(b, 2, 0, 0,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(b, 2, 0, 0,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(b, 2, 1, -1, Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(b, 2, 1, 0,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(b, 2, 1, 1,  Base.ImmutableDict(:type => :atomic))
            ],
        ),
        [a, b, b, a],
    )
end

const basisset_soc = let a = ChemicalSpecies(:H1), b = ChemicalSpecies(:H2)
    Quoll.BasisSetMetadata(
        dictionary([
            a => [
                Quoll.BasisMetadata(a, 1, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, -1, nothing)
                Quoll.BasisMetadata(a, 2, 1, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, 1, nothing)
                Quoll.BasisMetadata(a, 1, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, -1, nothing)
                Quoll.BasisMetadata(a, 2, 1, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, 1, nothing)
            ],
            b => [
                Quoll.BasisMetadata(b, 2, 0, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, -1, nothing)
                Quoll.BasisMetadata(b, 2, 1, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, 1, nothing)
                Quoll.BasisMetadata(b, 2, 0, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, -1, nothing)
                Quoll.BasisMetadata(b, 2, 1, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, 1, nothing)
            ],
        ]),
        [a, b, b, a],
    )
end

@testset "get_angular_momenta" begin
    ref_ells_a = [0, 0, 1]
    ref_ells_b = [0, 1]
    @test Quoll.get_angular_momenta(Quoll.basis_species(basisset, ChemicalSpecies(:H1))) ==
        ref_ells_a
    @test Quoll.get_angular_momenta(Quoll.basis_species(basisset, ChemicalSpecies(:H2))) ==
        ref_ells_b
end

@testset "get_unique_species" begin
    ref_unique_species = [ChemicalSpecies(:H1), ChemicalSpecies(:H2)]
    @test ref_unique_species == Quoll.get_unique_species(basisset)
end

@testset "get_atom2nbasis" begin
    ref_atom2nbasis = [5, 4, 4, 5]
    @test ref_atom2nbasis == Quoll.get_atom2nbasis(basisset)
end

@testset "get_atom2basis" begin
    ref_atom2basis = [1:5, 6:9, 10:13, 14:18]
    @test ref_atom2basis == Quoll.get_atom2basis(basisset)
end

@testset "get_basis2atom" begin
    ref_basis2atom = [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 4]
    @test ref_basis2atom == Quoll.get_basis2atom(basisset)
end

@testset "get_species2nbasis" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)

    ref_species2nbasis = Base.ImmutableDict(a => 5, b => 4)
    species2nbasis = Quoll.get_species2nbasis(basisset)

    @test keys(species2nbasis) == keys(ref_species2nbasis)
    for z in keys(ref_species2nbasis)
        @test species2nbasis[z] == ref_species2nbasis[z]
    end
end

@testset "basis_atom" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)

    basis_atom3 = [
        Quoll.BasisMetadata(b, 2, 0, 0, nothing)
        Quoll.BasisMetadata(b, 2, 1, -1, nothing)
        Quoll.BasisMetadata(b, 2, 1, 0, nothing)
        Quoll.BasisMetadata(b, 2, 1, 1, nothing)
    ]
    @test Quoll.basis_atom(basisset, 3) == basis_atom3

    basis_atom4 = [
        Quoll.BasisMetadata(a, 1, 0, 0, nothing)
        Quoll.BasisMetadata(a, 2, 0, 0, nothing)
        Quoll.BasisMetadata(a, 2, 1, -1, nothing)
        Quoll.BasisMetadata(a, 2, 1, 0, nothing)
        Quoll.BasisMetadata(a, 2, 1, 1, nothing)
    ]
    @test Quoll.basis_atom(basisset, 4) == basis_atom4
end

@testset "basis_species" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)

    basis_species_b = [
        Quoll.BasisMetadata(b, 2, 0, 0, nothing)
        Quoll.BasisMetadata(b, 2, 1, -1, nothing)
        Quoll.BasisMetadata(b, 2, 1, 0, nothing)
        Quoll.BasisMetadata(b, 2, 1, 1, nothing)
    ]
    @test Quoll.basis_species(basisset, b) == basis_species_b

    basis_species_a = [
        Quoll.BasisMetadata(a, 1, 0, 0, nothing)
        Quoll.BasisMetadata(a, 2, 0, 0, nothing)
        Quoll.BasisMetadata(a, 2, 1, -1, nothing)
        Quoll.BasisMetadata(a, 2, 1, 0, nothing)
        Quoll.BasisMetadata(a, 2, 1, 1, nothing)
    ]
    @test Quoll.basis_species(basisset, a) == basis_species_a
end

@testset "get_atom2offset" begin
    ref_atom2offset = [0, 5, 9, 13]
    @test ref_atom2offset == Quoll.get_atom2offset(basisset)
end

@testset "get_subbasis_masks" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)
    subbasis = [
        Quoll.BasisMetadata(a, 2, 0, 0, nothing)
        Quoll.BasisMetadata(b, 2, 0, 0, nothing)
        Quoll.BasisMetadata(b, 2, 1, 1, nothing)
        Quoll.BasisMetadata(a, 1, 0, 0, nothing)
    ]

    @testset "No spin" begin
        subbasis_masks_ref = dictionary([
            a => [1, 1, 0, 0, 0],
            b => [1, 0, 0, 1],
        ])
        @test Quoll.get_subbasis_masks(basisset, subbasis) == subbasis_masks_ref
    end

    # SOC case: for a given basis function in the subbasis, we want both up and down spin
    # basis functions to evaluate to true in the mask
    @testset "SOC" begin
        subbasis_masks_ref = dictionary([
            a => [1, 1, 0, 0, 0, 1, 1, 0, 0, 0],
            b => [1, 0, 0, 1, 1, 0, 0, 1],
        ])
        @test Quoll.get_subbasis_masks(basisset_soc, subbasis) == subbasis_masks_ref
    end
end

@testset "get_dense_subbasis_mask" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)
    subbasis = [
        Quoll.BasisMetadata(a, 2, 0, 0, nothing)
        Quoll.BasisMetadata(b, 2, 0, 0, nothing)
        Quoll.BasisMetadata(b, 2, 1, 1, nothing)
        Quoll.BasisMetadata(a, 1, 0, 0, nothing)
    ]
    dense_subbasis_mask_ref = [
        1, 1, 0, 0, 0,
        1, 0, 0, 1,
        1, 0, 0, 1,
        1, 1, 0, 0, 0,
    ]
    @test Quoll.get_dense_subbasis_mask(basisset, subbasis) == dense_subbasis_mask_ref
end

@testset "share_subblock" begin
    a = ChemicalSpecies(:H1)
    b1 = Quoll.BasisMetadata(a, 2, 0, 0, nothing)
    b2 = Quoll.BasisMetadata(a, 2, 1, -1, nothing)
    b3 = Quoll.BasisMetadata(a, 2, 1, 0, nothing)
    @test !Quoll.share_subblock(b1, b2)
    @test !Quoll.share_subblock(b1, b3)
    @test Quoll.share_subblock(b2, b3)
end

@testset "get_subblock_ranges" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)
    ranges_ref = dictionary([
        a => [1:1, 2:2, 3:5],
        b => [1:1, 2:4]
    ])
    @test Quoll.get_subblock_ranges(basisset) == ranges_ref
end

@testset "get_indices_in_subblock" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)
    indices_ref = dictionary([
        a => [1, 1, 1, 2, 3],
        b => [1, 1, 2, 3],
    ])
    @test Quoll.get_indices_in_subblock(basisset) == indices_ref
end

@testset "is_arbitrary_degeneracy" begin
    @test Quoll.is_arbitrary_degeneracy(basisset_degen.basis)
end

@testset "lift_arbitrary_degeneracy" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)
    basisset_degen_lifted_ref = Quoll.BasisSetMetadata(
        Base.ImmutableDict(
            a => [
                Quoll.BasisMetadata(a, 1, 0, 0,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(a, 2, 0, 0,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(a, 2, 1, -1, Base.ImmutableDict(:type => :atomic, :dgen => Symbol(1)))
                Quoll.BasisMetadata(a, 2, 1, 0,  Base.ImmutableDict(:type => :atomic, :dgen => Symbol(1)))
                Quoll.BasisMetadata(a, 2, 1, 1,  Base.ImmutableDict(:type => :atomic, :dgen => Symbol(1)))
                Quoll.BasisMetadata(a, 2, 1, -1, Base.ImmutableDict(:type => :atomic, :dgen => Symbol(2)))
                Quoll.BasisMetadata(a, 2, 1, 0,  Base.ImmutableDict(:type => :atomic, :dgen => Symbol(2)))
                Quoll.BasisMetadata(a, 2, 1, 1,  Base.ImmutableDict(:type => :atomic, :dgen => Symbol(2)))
            ],
            b => [
                Quoll.BasisMetadata(b, 2, 0, 0,  Base.ImmutableDict(:type => :atomic, :dgen => Symbol(1)))
                Quoll.BasisMetadata(b, 2, 0, 0,  Base.ImmutableDict(:type => :atomic, :dgen => Symbol(2)))
                Quoll.BasisMetadata(b, 2, 1, -1, Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(b, 2, 1, 0,  Base.ImmutableDict(:type => :atomic))
                Quoll.BasisMetadata(b, 2, 1, 1,  Base.ImmutableDict(:type => :atomic))
            ],
        ),
        [a, b, b, a],
    )
    basis_degen_lifted_ref = basisset_degen_lifted_ref.basis
    basis_degen_lifted = Quoll.lift_arbitrary_degeneracy(basisset_degen.basis)

    @test keys(basis_degen_lifted) == keys(basis_degen_lifted_ref)
    for z in keys(basis_degen_lifted_ref)
        bfuncs = basis_degen_lifted[z]
        bfuncs_ref = basis_degen_lifted_ref[z]
        @test bfuncs == bfuncs_ref
    end
end

@testset "reduce_basisset" begin
    a = ChemicalSpecies(:H1)
    b = ChemicalSpecies(:H2)
    subbasis = [
        Quoll.BasisMetadata(a, 2, 1, -1, nothing),
        Quoll.BasisMetadata(a, 2, 1, 0, nothing),
        Quoll.BasisMetadata(a, 2, 1, 1, nothing),
        Quoll.BasisMetadata(b, 2, 0, 0, nothing)
    ]
    basisset_reduced_ref = Quoll.BasisSetMetadata(
        Base.ImmutableDict(
            a => [
                Quoll.BasisMetadata(a, 2, 1, -1, nothing)
                Quoll.BasisMetadata(a, 2, 1, 0, nothing)
                Quoll.BasisMetadata(a, 2, 1, 1, nothing)
            ],
            b => [
                Quoll.BasisMetadata(b, 2, 0, 0, nothing)
            ],
        ),
        [a, b, b, a],
    )
    basisset_reduced = Quoll.reduce_basisset(basisset, subbasis; inverted=false)
    @test keys(basisset_reduced.basis) == keys(basisset_reduced_ref.basis)
    for z in keys(basisset_reduced_ref.basis)
        bfuncs = basisset_reduced.basis[z]
        bfuncs_ref = basisset_reduced_ref.basis[z]
        @test bfuncs == bfuncs_ref
    end

    # Inverted
    basisset_reduced_ref = Quoll.BasisSetMetadata(
        Base.ImmutableDict(
            a => [
                Quoll.BasisMetadata(a, 1, 0, 0, nothing)
                Quoll.BasisMetadata(a, 2, 0, 0, nothing)
            ],
            b => [
                Quoll.BasisMetadata(b, 2, 1, -1, nothing)
                Quoll.BasisMetadata(b, 2, 1, 0, nothing)
                Quoll.BasisMetadata(b, 2, 1, 1, nothing)
            ],
        ),
        [a, b, b, a],
    )
    basisset_reduced = Quoll.reduce_basisset(basisset, subbasis; inverted=true)
    @test keys(basisset_reduced.basis) == keys(basisset_reduced_ref.basis)
    for z in keys(basisset_reduced_ref.basis)
        bfuncs = basisset_reduced.basis[z]
        bfuncs_ref = basisset_reduced_ref.basis[z]
        @test bfuncs == bfuncs_ref
    end
end
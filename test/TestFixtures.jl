module TestFixtures

using Quoll
using AtomsBase
using StaticArrays
using Unitful
using Dictionaries

export make_test_atoms, make_test_basisset, make_test_block_sparsity,
    make_test_dense_recip_sparsity, make_test_spins,
    make_canonical_block_real_metadata, make_deeph_block_real_metadata,
    make_canonical_dense_recip_metadata,
    make_spin_real_metadata, make_spin_recip_metadata

### ATOMS ###

function make_test_atoms()
    cellvecs = SA[
        SA[5.0, 0.0, 0.0],
        SA[0.0, 5.0, 0.0],
        SA[0.0, 0.0, 5.0],
    ]u"Å"

    return periodic_system(
        [
            ChemicalSpecies(:H)  => SA[0.0, 0.0, 0.0]u"Å",
            ChemicalSpecies(:Li) => SA[2.5, 2.5, 2.5]u"Å",
        ],
        cellvecs,
    )
end

### BASIS SET ###

# Wiki-ordered basis: H has 1s + 2p (4 functions), Li has 1s + 2p (4 functions)
function make_test_basisset(atoms)
    H = ChemicalSpecies(:H)
    Li = ChemicalSpecies(:Li)

    basis = Base.ImmutableDict(
        H => [
            Quoll.BasisMetadata(H, 1, 0, 0, nothing),   # 1s
            Quoll.BasisMetadata(H, 2, 1, -1, nothing),   # 2p_{-1}
            Quoll.BasisMetadata(H, 2, 1, 0, nothing),    # 2p_0
            Quoll.BasisMetadata(H, 2, 1, 1, nothing),    # 2p_1
        ],
        Li => [
            Quoll.BasisMetadata(Li, 2, 0, 0, nothing),   # 2s
            Quoll.BasisMetadata(Li, 2, 1, -1, nothing),   # 2p_{-1}
            Quoll.BasisMetadata(Li, 2, 1, 0, nothing),    # 2p_0
            Quoll.BasisMetadata(Li, 2, 1, 1, nothing),    # 2p_1
        ],
    )

    atom2species = species(atoms, :)
    return Quoll.BasisSetMetadata(basis, atom2species)
end

### SPARSITY ###

function make_test_block_sparsity(; hermitian=false)
    if hermitian
        ij2images = dictionary([
            (1, 1) => [SA[0, 0, 0]],
            (1, 2) => [SA[0, 0, 0]],
            (2, 2) => [SA[0, 0, 0]],
        ])
    else
        ij2images = dictionary([
            (1, 1) => [SA[0, 0, 0]],
            (1, 2) => [SA[0, 0, 0]],
            (2, 1) => [SA[0, 0, 0]],
            (2, 2) => [SA[0, 0, 0]],
        ])
    end

    images = [SA[0, 0, 0]]
    return Quoll.BlockRealSparsity(ij2images, images, hermitian)
end

function make_test_dense_recip_sparsity(; hermitian=false)
    return Quoll.DenseRecipSparsity(hermitian)
end

### SPINS ###

function make_test_spins(basisset; spin=:up)
    kind = Quoll.Hamiltonian(; source=:ref, spin=spin)
    source = Quoll.CanonicalSource()
    return Quoll.SpinsMetadata(source, kind, basisset)
end

### METADATA FACTORIES ###

function make_canonical_block_real_metadata(; hermitian=false)
    atoms = make_test_atoms()
    basisset = make_test_basisset(atoms)
    sparsity = make_test_block_sparsity(; hermitian=hermitian)
    source = Quoll.CanonicalSource()
    shconv = Quoll.default_shconv(source)
    kind = Quoll.Hamiltonian(; source=:ref)

    basic = Quoll.BasicMetadataContainer(kind, source, sparsity, basisset, shconv, atoms)
    return Quoll.RealMetadata(basic)
end

function make_deeph_block_real_metadata(; hermitian=false)
    atoms = make_test_atoms()
    basisset = make_test_basisset(atoms)
    sparsity = make_test_block_sparsity(; hermitian=hermitian)
    source = Quoll.DeepHSource()
    shconv = Quoll.default_shconv(source)
    kind = Quoll.Hamiltonian(; source=:ref)

    # Reorder basis to DeepH SH convention
    basisset = Quoll.convert_basisset_shconv(basisset, shconv)

    basic = Quoll.BasicMetadataContainer(kind, source, sparsity, basisset, shconv, atoms)
    return Quoll.RealMetadata(basic)
end

function make_canonical_dense_recip_metadata(;
    hermitian=false, kpoint=SA[0.0, 0.0, 0.0],
)
    atoms = make_test_atoms()
    basisset = make_test_basisset(atoms)
    sparsity = make_test_dense_recip_sparsity(; hermitian=hermitian)
    source = Quoll.CanonicalSource()
    shconv = Quoll.default_shconv(source)
    kind = Quoll.Hamiltonian(; source=:ref)

    basic = Quoll.BasicMetadataContainer(kind, source, sparsity, basisset, shconv, atoms)
    return Quoll.RecipMetadata(basic, kpoint)
end

function make_spin_real_metadata(; hermitian=false, spin=:up)
    atoms = make_test_atoms()
    basisset = make_test_basisset(atoms)
    sparsity = make_test_block_sparsity(; hermitian=hermitian)
    source = Quoll.CanonicalSource()
    shconv = Quoll.default_shconv(source)
    kind = Quoll.Hamiltonian(; source=:ref, spin=spin)

    spins = Quoll.SpinsMetadata(source, kind, basisset)

    basic = Quoll.BasicMetadataContainer(kind, source, sparsity, basisset, shconv, atoms)
    return Quoll.SpinRealMetadata(basic, spins)
end

function make_spin_recip_metadata(;
    hermitian=false, kpoint=SA[0.0, 0.0, 0.0], spin=:up,
)
    atoms = make_test_atoms()
    basisset = make_test_basisset(atoms)
    sparsity = make_test_dense_recip_sparsity(; hermitian=hermitian)
    source = Quoll.CanonicalSource()
    shconv = Quoll.default_shconv(source)
    kind = Quoll.Hamiltonian(; source=:ref, spin=spin)

    spins = Quoll.SpinsMetadata(source, kind, basisset)

    basic = Quoll.BasicMetadataContainer(kind, source, sparsity, basisset, shconv, atoms)
    return Quoll.SpinRecipMetadata(basic, spins, kpoint)
end

end # module

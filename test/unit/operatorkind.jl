using Quoll

using Test

@testset "OperatorKind" begin
    h1 = Hamiltonian(source = :ref, spin = :none)
    @test h1 isa Quoll.OperatorKind{:Hamiltonian}
    @test h1.source == :ref
    @test h1.spin == :none

    h2 = Hamiltonian(Pair(:source, :ref), Pair(:spin, :none))
    @test h2.source == :ref
    @test h2.spin == :none
end

@testset "allows_symmetry" begin
    h1 = Hamiltonian(source = :ref, spin = :none)
    @test Quoll.allows_symmetry(h1) == true

    h2 = Hamiltonian(source = :ref, spin = :soc)
    @test Quoll.allows_symmetry(h2) == false

    s = Overlap(source = :ref)
    operatorkinds = [h2, s]
    @test Quoll.allows_symmetry(operatorkinds) == false
end

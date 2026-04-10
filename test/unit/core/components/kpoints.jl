using Quoll
using LinearAlgebra
using StaticArrays
using Unitful
using AtomsBase

using Test

@testset "check_kpoint_reduction" begin
    #     ┌─────────────────────────────────────────────────┐
    #  0.5│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⡀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⡧⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤⠤│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠂⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     │⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    # -0.5│⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⡇⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀│
    #     └─────────────────────────────────────────────────┘
    #     ⠀-0.5⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀⠀0.5⠀

    is_shift = trues(3)
    mesh = [2, 2, 2]
    grid_address = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [1, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [0, 1, 1],
        [1, 1, 1],
    ]
    k_all_reducible = [[1, 8], [2, 7], [3, 6], [4, 5]]
    #! format: off
    kirreds = [
        [ 1//4,  1//4,  1//4],
        [-1//4,  1//4,  1//4],
        [ 1//4, -1//4,  1//4],
        [-1//4, -1//4,  1//4],
    ]
    #! format: on

    rotations = [SMatrix{3,3}(I(3))]

    time_reversal = true
    @test Quoll.check_kpoint_reduction(
        rotations, time_reversal, is_shift, mesh, grid_address, k_all_reducible, kirreds
    ) == true

    time_reversal = false
    @test Quoll.check_kpoint_reduction(
        rotations, time_reversal, is_shift, mesh, grid_address, k_all_reducible, kirreds
    ) == false
end

@testset "construct_kgrid" begin
    cellvecs = SA[
        SA[10.0, 0.0, 0.0],
        SA[0.0, 10.0, 0.0],
        SA[0.0, 0.0, 10.0],
    ]u"Å"
    atoms = periodic_system(
        [
            ChemicalSpecies(:H1) => SA[0.00, 0.00, 0.00]u"Å",
            ChemicalSpecies(:H2) => SA[1.00, 1.00, 1.00]u"Å",
        ],
        cellvecs,
    )

    mesh = [2, 2, 2]
    shift = trues(3)

    @testset "No symmetries" begin
        kgrid = construct_kgrid(
            atoms; mesh=mesh, shift=shift, time_reversal=false, crystal_symmetry=false
        )
        @test length(kgrid.kpoints) == 8
        @test sum(kgrid.weights) ≈ 1.0
    end

    @testset "Time reversal" begin
        kgrid = construct_kgrid(
            atoms; mesh=mesh, shift=shift, time_reversal=true, crystal_symmetry=false
        )
        @test length(kgrid.kpoints) == 4
        @test sum(kgrid.weights) ≈ 1.0
    end

    @testset "All symmetries" begin
        kgrid = construct_kgrid(
            atoms; mesh=mesh, shift=shift, time_reversal=true, crystal_symmetry=true
        )
        @test length(kgrid.kpoints) == 2
        @test sum(kgrid.weights) ≈ 1.0
    end
end

@testset "get_recip_mesh" begin
    real_lattice = 2π * [
        1.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 1.0;
    ]
    density = 5
    mesh_ref = [5, 5, 5]
    @test Quoll.get_recip_mesh(real_lattice, density) == mesh_ref
end

@testset "precompute_phases" begin
    ks = [[0, 0, 0], [0, 0, 0.5], [0, 0, -0.5]]
    ts = [[0, 0, 0], [0, 0, 1], [0, 0, -1]]
    phases_ref = [
        exp(2π * im * 0.0) exp(2π * im *  0.0) exp(2π * im *  0.0);
        exp(2π * im * 0.0) exp(2π * im *  0.5) exp(2π * im * -0.5);
        exp(2π * im * 0.0) exp(2π * im * -0.5) exp(2π * im *  0.5);
    ]
    @test Quoll.precompute_phases(ks, ts) ≈ phases_ref
end
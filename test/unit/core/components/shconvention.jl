#! format: off
using Quoll
using Dictionaries
using AtomsBase

using Test

function reorder_matrix_type1(A1, basis1, basis2, shconv)
    A2 = zero(A1)
    
    iis1 = Quoll.get_indices_in_subblock(basis1)
    iis2 = Quoll.get_indices_in_subblock(basis2)

    for (j, b2) in enumerate(basis2)
        shift_j = Quoll.get_shift(b2.l, iis2[j], shconv)
        phase_j = Quoll.get_phase(b2.l, iis2[j], shconv)

        for (i, b1) in enumerate(basis1)
            shift_i = Quoll.get_shift(b1.l, iis1[i], shconv)
            phase_i = Quoll.get_phase(b1.l, iis1[i], shconv)

            A2[i, j] = A1[i + shift_i, j + shift_j] * phase_i * phase_j
        end
    end

    return A2
end

function reorder_matrix_type2(A1, basis1, basis2, inv_shconv)
    A2 = zero(A1)

    iis1 = Quoll.get_indices_in_subblock(basis1)
    iis2 = Quoll.get_indices_in_subblock(basis2)

    for (j, b2) in enumerate(basis2)
        shift_j = Quoll.get_shift(b2.l, iis2[j], inv_shconv)
        phase_j = Quoll.get_phase(b2.l, iis2[j], inv_shconv)

        for (i, b1) in enumerate(basis1)
            shift_i = Quoll.get_shift(b1.l, iis1[i], inv_shconv)
            phase_i = Quoll.get_phase(b1.l, iis1[i], inv_shconv)

            # Note shifts on A2 instead of A1 as in type1 conversion
            A2[i + shift_i, j + shift_j] = A1[i, j] * phase_i * phase_j
        end
    end

    return A2
end

@testset "SHConvention" begin
    orders = [[1], [3,  1,  2], [1,  4,  5,  3,  2]]
    shphases = [[1], [1, -1,  1], [1,  1, -1, -1,  1]]
    shconv = Quoll.SHConvention(orders, shphases)

    l = 0
    @test shconv.orders[l + 1][1] == 1
    @test shconv.shifts[l + 1][1] == 1 - (l + 1) - 0
    @test shconv.phases[l + 1][1] == 1

    l = 1
    @test shconv.orders[l + 1][1] == 3
    @test shconv.orders[l + 1][2] == 1
    @test shconv.orders[l + 1][3] == 2
    @test shconv.shifts[l + 1][1] == 3 - (l + 1) + 1
    @test shconv.shifts[l + 1][2] == 1 - (l + 1)
    @test shconv.shifts[l + 1][3] == 2 - (l + 1) - 1
    @test shconv.phases[l + 1][1] == 1
    @test shconv.phases[l + 1][2] == -1
    @test shconv.phases[l + 1][3] == 1

    l = 2
    @test shconv.orders[l + 1][1] == 1
    @test shconv.orders[l + 1][2] == 4
    @test shconv.orders[l + 1][3] == 5
    @test shconv.orders[l + 1][4] == 3
    @test shconv.orders[l + 1][5] == 2
    @test shconv.shifts[l + 1][1] == 1 - (l + 1) + 2
    @test shconv.shifts[l + 1][2] == 4 - (l + 1) + 1
    @test shconv.shifts[l + 1][3] == 5 - (l + 1)
    @test shconv.shifts[l + 1][4] == 3 - (l + 1) - 1
    @test shconv.shifts[l + 1][5] == 2 - (l + 1) - 2
    @test shconv.phases[l + 1][1] == 1
    @test shconv.phases[l + 1][2] == 1
    @test shconv.phases[l + 1][3] == -1
    @test shconv.phases[l + 1][4] == -1
    @test shconv.phases[l + 1][5] == 1
end

@testset "Reorder" begin
    basis = [
        Quoll.BasisMetadata(ChemicalSpecies(:H), 2, 0, 0, nothing),
        Quoll.BasisMetadata(ChemicalSpecies(:H), 2, 1, -1, nothing),
        Quoll.BasisMetadata(ChemicalSpecies(:H), 2, 1, 0, nothing),
        Quoll.BasisMetadata(ChemicalSpecies(:H), 2, 1, 1, nothing),
    ]

    A1 = [
        0 1 2 3;
        1 2 3 4;
        2 3 4 5;
        3 4 5 6;
    ]

    # 1  2  3 start
    # 2 -1  3 end
    orders1 = [[1], [2,  1,  3]]
    shphases1 = [[1], [1, -1,  1]]
    shconv1 = Quoll.SHConvention(orders1, shphases1)

    # 2 -1  3 start
    # 1  3  2 end
    orders2 = [[1], [ 2,  3,  1]]
    shphases2 = [[1], [-1,  1,  1]]
    shconv2 = Quoll.SHConvention(orders2, shphases2)

    @testset "Basic" begin
        ref = [
            0  2 -1  3;
            2  4 -3  5;
           -1 -3  2 -4;
            3  5 -4  6;
        ]
        A2 = reorder_matrix_type1(A1, basis, basis, shconv1)
        @test A2 == ref

        P = Quoll.compute_signed_perm_matrix(basis, shconv1)
        A2 = P * A1 * P'
        @test A2 == ref
    end

    @testset "Combination" begin
        ref = [
            0 1 3 2;
            1 2 4 3;
            3 4 6 5;
            2 3 5 4;
        ]

        A2 = reorder_matrix_type1(A1, basis, basis, shconv1)
        A3 = reorder_matrix_type1(A2, basis, basis, shconv2)
        @test A3 == ref

        A3_direct = reorder_matrix_type1(A1, basis, basis, shconv2 ∘ shconv1)
        @test A3_direct == ref
    end

    @testset "Inverse" begin
        A2 = reorder_matrix_type1(A1, basis, basis, shconv1)
        A3 = reorder_matrix_type1(A2, basis, basis, inv(shconv1))
        @test A1 == A3

        A2type1 = reorder_matrix_type1(A1, basis, basis, shconv1)
        A2type2 = reorder_matrix_type2(A1, basis, basis, inv(shconv1))
        @test A2type1 == A2type2

        identity_shconv = inv(shconv1) ∘ shconv1

        l = 0
        @test identity_shconv.orders[l + 1][1] == 1
        @test identity_shconv.shifts[l + 1][1] == 0
        @test identity_shconv.phases[l + 1][1] == 1

        l = 1
        @test identity_shconv.orders[l + 1][1] == 1
        @test identity_shconv.orders[l + 1][2] == 2
        @test identity_shconv.orders[l + 1][3] == 3
        @test identity_shconv.shifts[l + 1][1] == 0
        @test identity_shconv.shifts[l + 1][2] == 0
        @test identity_shconv.shifts[l + 1][3] == 0
        @test identity_shconv.phases[l + 1][1] == 1
        @test identity_shconv.phases[l + 1][2] == 1
        @test identity_shconv.phases[l + 1][3] == 1
    end

end

@testset "Precompute" begin
    z1 = ChemicalSpecies(:H1)
    z2 = ChemicalSpecies(:H2)
    basis = dictionary([
        z1 => [
            Quoll.BasisMetadata(z1, 2, 0, 0, nothing),
            Quoll.BasisMetadata(z1, 2, 1, -1, nothing),
            Quoll.BasisMetadata(z1, 2, 1, 0, nothing),
            Quoll.BasisMetadata(z1, 2, 1, 1, nothing),
        ],
        z2 => [
            Quoll.BasisMetadata(z2, 1, 0, 0, nothing),
        ]
    ])

    basisset = Quoll.BasisSetMetadata(basis, [z2, z1])

    orders = [[1], [2,  1,  3]]
    shphases = [[1], [1, -1,  1]]
    shconv = Quoll.SHConvention(orders, shphases)

    @testset "precompute_shifts" begin
        ref = dictionary([
            z1 => [0, 1, -1, 0]
            z2 => [0]
        ])
        shifts = Quoll.precompute_shifts(basisset, shconv)
        @test Set(keys(shifts)) == Set(keys(ref))
        for key in keys(ref)
            @test shifts[key] == ref[key]
        end
    end

    @testset "precompute_orders" begin
        ref = dictionary([
            z1 => [1, 3, 2, 4]
            z2 => [1]
        ])
        orders = Quoll.precompute_orders(basisset, shconv)
        @test Set(keys(orders)) == Set(keys(ref))
        for key in keys(ref)
            @test orders[key] == ref[key]
        end
    end

    @testset "precompute_shphases" begin
        ref = dictionary([
            z1 => [1, 1, -1, 1]
            z2 => [1]
        ])
        shphases = Quoll.precompute_shphases(basisset, shconv, Val(1))
        @test Set(keys(shphases)) == Set(keys(ref))
        for key in keys(ref)
            @test shphases[key] == ref[key]
        end

        ref = dictionary([
            (z1, z1) => [
                1  1 -1  1;
                1  1 -1  1;
               -1 -1  1 -1;
                1  1 -1  1;
            ],
            (z1, z2) => [
                1;
                1;
               -1;
                1;;
            ],
            (z2, z1) => [
                1  1 -1  1;
            ],
            (z2, z2) => [
                1;;
            ]
        ])
        shphases = Quoll.precompute_shphases(basisset, shconv, Val(2))
        @test Set(keys(shphases)) == Set(keys(ref))
        for key in keys(ref)
            @test shphases[key] == ref[key]
        end
    end
    
    @testset "precompute_signed_perm_matrices" begin
        ref = dictionary([
            z1 => [
                 1  0  0  0;
                 0  0  1  0;
                 0 -1  0  0;
                 0  0  0  1
            ],
            z2 => [
                 1;;
            ]
        ])

        signed_perm_matrices = Quoll.precompute_signed_perm_matrices(basisset, shconv, Int)
        @test Set(keys(signed_perm_matrices)) == Set(keys(ref))
        for key in keys(ref)
            @test signed_perm_matrices[key] == ref[key]
        end
    end
    
end
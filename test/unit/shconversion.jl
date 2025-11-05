using Quoll
using Dictionaries
using AtomsBase

using Test

function reorder_matrix_type1(A1, basis, shconv)
    A2 = zero(A1)

    # Precompute shifts and phases
    shifts = Quoll.get_shift.(basis, Ref(shconv))
    phases = Quoll.get_phase.(basis, Ref(shconv))

    for j in axes(A2, 2)
        shift_j = shifts[j]
        phase_j = phases[j]

        for i in axes(A2, 1)
            shift_i = shifts[i]
            phase_i = phases[i]

            A2[i, j] = A1[i + shift_i, j + shift_j] * phase_i * phase_j
        end
    end

    return A2
end

function reorder_matrix_type2(A1, basis, inv_shconv)
    A2 = zero(A1)

    # Precompute shifts and phases
    shifts = Quoll.get_shift.(basis, Ref(inv_shconv))
    phases = Quoll.get_phase.(basis, Ref(inv_shconv))

    for j in axes(A2, 2) 
        shift_j = shifts[j]
        phase_j = phases[j]

        for i in axes(A2, 1) 
            shift_i = shifts[i]
            phase_i = phases[i]

            # Note shifts on A2 instead of A1 as in type1 conversion
            A2[i + shift_i, j + shift_j] = A1[i, j] * phase_i * phase_j
        end
    end

    return A2
end

function reorder_matrix_type1_orders(A1, basis, shconv)
    A2 = zero(A1)

    # Precompute orders and phases
    orders = Quoll.get_shift.(basis, Ref(shconv)) .+ collect(1:length(basis))
    phases = Quoll.get_phase.(basis, Ref(shconv))

    for j in axes(A2, 2) 
        order_j = orders[j]
        phase_j = phases[j]

        for i in axes(A2, 1) 
            order_i = orders[i]
            phase_i = phases[i]

            A2[i, j] = A1[order_i, order_j] * phase_i * phase_j
        end
    end

    return A2
end


@testset "SHConversion" begin
    orders = [[1], [3,  1,  2], [1,  4,  5,  3,  2]]
    phases = [[1], [1, -1,  1], [1,  1, -1, -1,  1]]
    shconv = Quoll.SHConversion(orders, phases)

    l = 0
    @test shconv.orders[l + 1][0] == 1 - (l + 1)
    @test shconv.shifts[l + 1][0] == 1 - (l + 1) - 0
    @test shconv.phases[l + 1][0] == 1

    l = 1
    @test shconv.orders[l + 1][-1] == 3 - (l + 1)
    @test shconv.orders[l + 1][0] == 1 - (l + 1)
    @test shconv.orders[l + 1][1] == 2 - (l + 1)
    @test shconv.shifts[l + 1][-1] == 3 - (l + 1) + 1
    @test shconv.shifts[l + 1][0] == 1 - (l + 1)
    @test shconv.shifts[l + 1][1] == 2 - (l + 1) - 1
    @test shconv.phases[l + 1][-1] == 1
    @test shconv.phases[l + 1][0] == -1
    @test shconv.phases[l + 1][1] == 1

    l = 2
    @test shconv.orders[l + 1][-2] == 1 - (l + 1)
    @test shconv.orders[l + 1][-1] == 4 - (l + 1)
    @test shconv.orders[l + 1][0] == 5 - (l + 1)
    @test shconv.orders[l + 1][1] == 3 - (l + 1)
    @test shconv.orders[l + 1][2] == 2 - (l + 1)
    @test shconv.shifts[l + 1][-2] == 1 - (l + 1) + 2
    @test shconv.shifts[l + 1][-1] == 4 - (l + 1) + 1
    @test shconv.shifts[l + 1][0] == 5 - (l + 1)
    @test shconv.shifts[l + 1][1] == 3 - (l + 1) - 1
    @test shconv.shifts[l + 1][2] == 2 - (l + 1) - 2
    @test shconv.phases[l + 1][-2] == 1
    @test shconv.phases[l + 1][-1] == 1
    @test shconv.phases[l + 1][0] == -1
    @test shconv.phases[l + 1][1] == -1
    @test shconv.phases[l + 1][2] == 1
end

@testset "Reorder" begin
    basis = [
        BasisMetadata(ChemicalSpecies(:H), 2, 0, 0, nothing),
        BasisMetadata(ChemicalSpecies(:H), 2, 1, -1, nothing),
        BasisMetadata(ChemicalSpecies(:H), 2, 1, 0, nothing),
        BasisMetadata(ChemicalSpecies(:H), 2, 1, 1, nothing),
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
    phases1 = [[1], [1, -1,  1]]
    shconv1 = Quoll.SHConversion(orders1, phases1)

    # 2 -1  3 start
    # 1  3  2 end
    orders2 = [[1], [ 2,  3,  1]]
    phases2 = [[1], [-1,  1,  1]]
    shconv2 = Quoll.SHConversion(orders2, phases2)

    @testset "Basic" begin
        ref = [
            0  2 -1  3;
            2  4 -3  5;
           -1 -3  2 -4;
            3  5 -4  6;
        ]
        A2 = reorder_matrix_type1(A1, basis, shconv1)
        @test A2 == ref

        A2 = reorder_matrix_type1_orders(A1, basis, shconv1)
        @test A2 == ref
    end

    @testset "Combination" begin
        ref = [
            0 1 3 2;
            1 2 4 3;
            3 4 6 5;
            2 3 5 4;
        ]

        A2 = reorder_matrix_type1(A1, basis, shconv1)
        A3 = reorder_matrix_type1(A2, basis, shconv2)
        @test A3 == ref

        A3_direct = reorder_matrix_type1(A1, basis, shconv2 ∘ shconv1)
        @test A3_direct == ref
    end

    @testset "Inverse" begin
        A2 = reorder_matrix_type1(A1, basis, shconv1)
        A3 = reorder_matrix_type1(A2, basis, inv(shconv1))
        @test A1 == A3

        A2type1 = reorder_matrix_type1(A1, basis, shconv1)
        A2type2 = reorder_matrix_type2(A1, basis, inv(shconv1))
        @test A2type1 == A2type2

        identity_shconv = inv(shconv1) ∘ shconv1

        l = 0
        @test identity_shconv.orders[l + 1][0] == 1 - (l + 1)
        @test identity_shconv.shifts[l + 1][0] == 0
        @test identity_shconv.phases[l + 1][0] == 1

        l = 1
        @test identity_shconv.orders[l + 1][-1] == 1 - (l + 1)
        @test identity_shconv.orders[l + 1][0] == 2 - (l + 1)
        @test identity_shconv.orders[l + 1][1] == 3 - (l + 1)
        @test identity_shconv.shifts[l + 1][-1] == 0
        @test identity_shconv.shifts[l + 1][0] == 0
        @test identity_shconv.shifts[l + 1][1] == 0
        @test identity_shconv.phases[l + 1][-1] == 1
        @test identity_shconv.phases[l + 1][0] == 1
        @test identity_shconv.phases[l + 1][1] == 1
    end

end

@testset "Precompute" begin
    z1 = ChemicalSpecies(:H1)
    z2 = ChemicalSpecies(:H2)
    basis = dictionary([
        z1 => [
            BasisMetadata(z1, 2, 0, 0, nothing),
            BasisMetadata(z1, 2, 1, -1, nothing),
            BasisMetadata(z1, 2, 1, 0, nothing),
            BasisMetadata(z1, 2, 1, 1, nothing),
        ],
        z2 => [
            BasisMetadata(z2, 1, 0, 0, nothing),
        ]
    ])

    basisset = Quoll.BasisSetMetadata(basis, [z2, z1])

    orders = [[1], [2,  1,  3]]
    phases = [[1], [1, -1,  1]]
    shconv = Quoll.SHConversion(orders, phases)

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

    @testset "precompute_phases" begin
        ref = dictionary([
            z1 => [1, 1, -1, 1]
            z2 => [1]
        ])
        phases = Quoll.precompute_phases(basisset, shconv, DIM = Val(1))
        @test Set(keys(phases)) == Set(keys(ref))
        for key in keys(ref)
            @test phases[key] == ref[key]
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
        phases = Quoll.precompute_phases(basisset, shconv, DIM = Val(2))
        @test Set(keys(phases)) == Set(keys(ref))
        for key in keys(ref)
            @test phases[key] == ref[key]
        end
    end
    
end
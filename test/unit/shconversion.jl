using Quoll
using AtomsBase

using Test

function reorder_matrix(A1, basis, shconv)
    A2 = zero(A1)

    # Precompute shifts and phases
    n = size(A1, 2)
    shifts = Vector{Int}(undef, n)
    phases = Vector{Int}(undef, n)

    for j in 1:n
        shifts[j], phases[j] = Quoll.get_shiftphase(basis[j], shconv)
    end

    for j in 1:n
        shift_j = shifts[j]
        phase_j = phases[j]

        for i in 1:n
            shift_i = shifts[i]
            phase_i = phases[i]

            A2[i, j] = A1[i + shift_i, j + shift_j] * phase_i * phase_j
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
    @test shconv.phases[l + 1][0] == 1

    l = 1
    @test shconv.orders[l + 1][-1] == 3 - (l + 1)
    @test shconv.orders[l + 1][0] == 1 - (l + 1)
    @test shconv.orders[l + 1][1] == 2 - (l + 1)
    @test shconv.phases[l + 1][-1] == 1
    @test shconv.phases[l + 1][0] == -1
    @test shconv.phases[l + 1][1] == 1

    l = 2
    @test shconv.orders[l + 1][-2] == 1 - (l + 1)
    @test shconv.orders[l + 1][-1] == 4 - (l + 1)
    @test shconv.orders[l + 1][0] == 5 - (l + 1)
    @test shconv.orders[l + 1][1] == 3 - (l + 1)
    @test shconv.orders[l + 1][2] == 2 - (l + 1)
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
        A2 = reorder_matrix(A1, basis, shconv1)
        @test A2 == [
            0  2 -1  3;
            2  4 -3  5;
           -1 -3  2 -4;
            3  5 -4  6;
        ]
    end

    @testset "Combination" begin
        ref = [
            0 1 3 2;
            1 2 4 3;
            3 4 6 5;
            2 3 5 4;
        ]

        A2 = reorder_matrix(A1, basis, shconv1)
        A3 = reorder_matrix(A2, basis, shconv2)
        @test A3 == ref

        A3_direct = reorder_matrix(A1, basis, shconv2 ∘ shconv1)
        @test A3_direct == ref
    end

    @testset "Inverse" begin
        A2 = reorder_matrix(A1, basis, shconv1)
        A3 = reorder_matrix(A2, basis, inv(shconv1))
        @test A1 == A3

        identity_shconv = inv(shconv1) ∘ shconv1

        l = 0
        @test identity_shconv.orders[l + 1][0] == 1 - (l + 1)
        @test identity_shconv.phases[l + 1][0] == 1

        l = 1
        @test identity_shconv.orders[l + 1][-1] == 1 - (l + 1)
        @test identity_shconv.orders[l + 1][0] == 2 - (l + 1)
        @test identity_shconv.orders[l + 1][1] == 3 - (l + 1)
        @test identity_shconv.phases[l + 1][-1] == 1
        @test identity_shconv.phases[l + 1][0] == 1
        @test identity_shconv.phases[l + 1][1] == 1
    end

end
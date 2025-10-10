using Quoll

using Test

@testset "FHIaimsLAPACKSplit" begin
    split_method = Quoll.FHIaimsLAPACKSplit()

    @testset "N_mpi < N" begin
        N_mpi = 2
        N = 4

        my_rank = 0
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [2, 4]
        my_rank = 1
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [1, 3]
    end

    @testset "N_mpi = N" begin
        N_mpi = 4
        N = 4

        my_rank = 0
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [4]
        my_rank = 1
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [1]
        my_rank = 2
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [2]
        my_rank = 3
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [3]
    end

    @testset "N_mpi > N" begin

        N_mpi = 8
        N = 4

        my_rank = 0
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == Int[]
        my_rank = 1
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [1]
        my_rank = 2
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [2]
        my_rank = 3
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [3]
        my_rank = 4
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [4]
        my_rank = 5
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == Int[]
        my_rank = 6
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == Int[]
        my_rank = 7
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == Int[]

    end

end

@testset "DefaultSplit" begin
    split_method = Quoll.DefaultSplit()

    @testset "N_mpi < N" begin
        N_mpi = 2
        N = 4

        my_rank = 0
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [1, 2]
        my_rank = 1
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [3, 4]
    end

    @testset "N_mpi = N" begin
        N_mpi = 4
        N = 4

        my_rank = 0
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [1]
        my_rank = 1
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [2]
        my_rank = 2
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [3]
        my_rank = 3
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [4]
    end

    @testset "N_mpi > N" begin

        N_mpi = 8
        N = 4

        my_rank = 0
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [1]
        my_rank = 1
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [2]
        my_rank = 2
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [3]
        my_rank = 3
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [4]
        my_rank = 4
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [1]
        my_rank = 5
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [2]
        my_rank = 6
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [3]
        my_rank = 7
        @test Quoll.split_work(N, N_mpi, my_rank, split_method) == [4]

    end

end
# TODO: some notes, move somewhere else
# GLOBAL COMMUNICATOR
# comm
# comm_size (can be obtained from comm)
# my_rank (can be obtained from comm)

# my_configs (can be obtained from (my_rank, comm_size, n_structures))

# SINGLE CONFIG COMMUNICATOR
# (if more than one config per mpi task, then subcomm_size = 1)
# subcomm
# subcomm_size (can be obtained from subcomm)
# my_subrank (can be obtained from subcomm)

# my_kpoints (can be obtained from (my_subrank, subcomm_size, n_kpoints))

# SINGLE KPOINT COMMUNICATOR
# (if more than one k-point per mpi task, then subsubcomm_size = 1)
# - AFAIK we don't need this one because we will not be implementing ScaLAPACK for now
# - Instead of this, it's better to initialise the script with fewer 
# subsubcomm
# subsubcomm_size (can be obtained from subsubcomm)
# my_subsubrank (can be obtained from subsubcomm)

abstract type AbstractWorkSplit end

struct DefaultSplit <: AbstractWorkSplit end
struct FHIaimsLAPACKSplit <: AbstractWorkSplit end

function split_work(N::Integer, comm::MPI.Comm, split_method::AbstractWorkSplit)
    N_comm = MPI.Comm_size(comm)
    my_rank = MPI.Comm_rank(comm)
    return split_work(N, N_comm, my_rank, split_method)
end

function split_work(N::Integer, N_comm::Integer, my_rank::Integer, ::DefaultSplit)
    min_N_pp = fld(N, N_comm)
    remainder = mod(N, N_comm)

    if N ≥ N_comm
        i_begin = min_N_pp * my_rank + min(my_rank, remainder) + 1
        i_end = min_N_pp * (my_rank + 1) + min(my_rank + 1, remainder)
    else
        i_begin = i_end = mod(my_rank, N) + 1
    end

    return collect(i_begin:i_end)
end

function split_work(N::Integer, N_comm::Integer, my_rank::Integer, ::FHIaimsLAPACKSplit)
    if my_rank != 0
        my_offset = my_rank
    else
        my_offset = N_comm
    end

    my_work = Int[]
    for max_N_pp in 1:((N - 1) / N_comm + 1)
        i = (max_N_pp - 1) * N_comm + my_offset
        if i <= N
            push!(my_work, i)
        end
    end

    return my_work
end
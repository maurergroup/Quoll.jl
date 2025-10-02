using Logging
using LoggingExtras
using Dates

const JULIA_QUOLL_MPI_LOGALLRANKS = parse(Bool, get(ENV, "JULIA_QUOLL_MPI_LOGALLRANKS", "false"))

timestamp_rank_logger(logger, rank) = TransformerLogger(logger) do log
    # This approach is not ideal because kwargs are printed into the message
    # tmp = pairs(NamedTuple(vcat(collect(log.kwargs), [:rank => rank])))
    # log = merge(log, (; kwargs = tmp))
    merge(log, (; message = "[$(Dates.format(now(), """yyyy-mm-dd HH:MM:SS"""))|Rank $(rank)] $(log.message)"))
end

filterrank_logger(logger, rank) = EarlyFilteredLogger(logger) do log
    if JULIA_QUOLL_MPI_LOGALLRANKS || rank == 0
        return true
    else
        return false
    end
end

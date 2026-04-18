using Logging
using LoggingExtras
using Dates

timestamp_rank_logger(logger, rank) = TransformerLogger(logger) do log
    # This approach is not ideal because kwargs are printed into the message
    # tmp = pairs(NamedTuple(vcat(collect(log.kwargs), [:rank => rank])))
    # log = merge(log, (; kwargs = tmp))
    merge(log, (; message = "[$(Dates.format(now(), """yyyy-mm-dd HH:MM:SS"""))|Rank $(rank)] $(log.message)"))
end

filterrank_logger(logger, rank; all_ranks::Bool=false) = EarlyFilteredLogger(logger) do log
    if all_ranks || rank == 0
        return true
    else
        return false
    end
end

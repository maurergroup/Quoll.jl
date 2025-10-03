# mpiexecjl -n 8 -t 1 quoll.jl input_file.toml
# (non-MPI version works with julia quoll.jl even when MPI is not properly set up)
# In theory I probably need a separate Project.toml in bin/ because we need
# additional dependencies. Ideally this could be part of [targets], but
# currently only `test` and `build` targets are supported by Pkg.
# Maybe it would make sense to simply have Quoll as a dependency like in docs
# and add extra bits like parsing and logging dependencies. In that case I would
# probably need separate tests as well. However, only units test then would be 
# testing functions in src/ (maybe that's fine?)
using Quoll
using Configurations
using LinearAlgebra
using MPI

### START MPI ###
MPI.Init()
BLAS.set_num_threads(1)

global_comm = MPI.COMM_WORLD
global_rank = MPI.Comm_rank(global_comm)

### SETUP LOGGER ###
# - I think it's fair to simply move logs to stderr because if quoll job
#   is submitted on hpc normally stderr goes to a file like a slurm file.
# - However, it would be nice to make sure logging works as expected
#   for MPI, e.g. could use a Tranformer logger to only show rank 0 messages
#   (could do this with an environment variable instead of input file to allow
#   immediate logging even during input file reading).
# - Also could add timestamps to logging messages.
# - If we keep logger here, maybe we should keep the parser here as well?
#   For example, Configurations.jl is only needed for parsing.
include("logger.jl")
global_logger(
    filterrank_logger(
        timestamp_rank_logger(
            ConsoleLogger(stderr, Logging.Info),
            global_rank
        ),
        global_rank,
    )
)

@info "Parsing QUOLL input file"

input_filepath = abspath(ARGS[1])
params = from_toml(Quoll.Parser.QuollParams, input_filepath)

@info "Splitting configurations across MPI tasks"

N_dirs = length(params.input.directories)
my_idirs = Quoll.MPITools.split_work(N_dirs, global_comm, Quoll.MPITools.DefaultSplit())


for idir in my_idirs
    dir = params.input.directories[idir]
    @info "Starting configuration $(basename(dir))"

    @info "Loading atoms"
    atoms = load_atoms(dir, params.input.format)

end
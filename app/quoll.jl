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

include("logger.jl")

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
params, input_dirs, output_dirs = Quoll.Parser.parse_inputfile(input_filepath)

@info "Splitting configurations across MPI tasks"

ndirs = length(input_dirs)
my_idirs = Quoll.split_work(ndirs, global_comm, Quoll.DefaultSplit())

# Construct an MPI sub-communicator for MPI tasks working on the same atomic configuration
global_color = first(my_idirs)
config_comm = MPI.Comm_split(global_comm, global_color, global_rank)
config_rank = MPI.Comm_rank(config_comm)

for idir in my_idirs
    input_dir = input_dirs[idir]

    @info "Starting configuration $(basename(input_dir))"
    @info "Loading operators"
    
    operatorkinds = find_operatorkinds(input_dir, params)
    operators = load_operators(input_dir, operatorkinds, params.input.format)

    @info "Converting operators into BSparseOperator format"

    operators = convert_operators(operators, BSparseOperator)

    if Quoll.Parser.requires_kpoint_grid(params)
        @info "Initialising k-point grid"

        # Assuming all operators in the directory are related to the same atoms
        atoms = Quoll.get_atoms(first(operators))
        kgrid = get_kgrid(atoms, operatorkinds, params)

        @info "Splitting k-points across MPI tasks"
        # TODO: add an option for different splits in the input file
        nkpoints = Quoll.get_nkpoints(kgrid)
        my_ikpoints = Quoll.split_work(nkpoints, config_comm, Quoll.FHIaimsLAPACKSplit())
    end

    if !isnothing(params.basis_projection)
        @info "Projecting the core states"

        operators = perform_core_projection(operators, kgrid, my_ikpoints, config_comm, params)
    end

    if !isnothing(params.output)
        @info "Converting operators into $(params.output.format) format"

        operators = convert_operators(operators, params.output.format, hermitian = params.output.hermitian)

        @info "Writing operators"

        if global_rank == 0
            write_operators(output_dirs[idir], operators)
        end
    end

end

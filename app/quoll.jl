using Quoll
using Configurations
using LinearAlgebra
using MPI

include("logger.jl")

### START MPI ###
# Set BLAS threads to 1 to prioritise MPI parallelisation
MPI.Init()
BLAS.set_num_threads(1)

global_comm = MPI.COMM_WORLD
global_rank = MPI.Comm_rank(global_comm)

### SETUP LOGGER ###
global_logger(
    filterrank_logger(
        timestamp_rank_logger(
            ConsoleLogger(stderr, Logging.Info),
            global_rank,
        ),
        global_rank,
    ),
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

# TODO: Add options for operator types that are passed into `load_operators` and
# `convert_operator`, this would be important if a particular conversion method needs
# to be selectedm, e.g. the conversion between Operator and KeyedOperator might not be
# implemented

for idir in my_idirs
    input_dir = input_dirs[idir]

    @info "Starting configuration $(basename(input_dir))"
    @info "Loading operators"

    operatorkinds = find_operatorkinds(input_dir, params)
    operators = load_operators(
        Quoll.Operator, params.input.format, input_dir, operatorkinds
    )

    @info "Converting operators to canonical format"
    operators =
        convert_operator.(
            Quoll.KeyedOperator, Quoll.CanonicalBlockRealMetadata, operators;
            radii=params.input.radii,
        )

    if Quoll.Parser.requires_kpoint_grid(params)
        @info "Initialising k-point grid"

        # Assuming all operators in the directory are related to the same atoms
        atoms = Quoll.op_atoms(first(operators))
        kgrid = construct_kgrid(atoms, operatorkinds, params)

        @info "Splitting k-points across MPI tasks"
        # TODO: Add an option to change the k-point split, e.g. would be required
        # if ScaLAPACK was used
        nkpoints = Quoll.get_nkpoints(kgrid)
        my_ikpoints = Quoll.split_work(nkpoints, config_comm, Quoll.FHIaimsLAPACKSplit())
    end

    if !isnothing(params.basis_projection)
        @info "Projecting the core states"

        operators = perform_core_projection(
            operators, kgrid, my_ikpoints, config_comm, params
        )
    end

    if !isnothing(params.output)
        @info "Converting operators into $(params.output.format) format"

        operators =
            convert_operator.(
                Quoll.Operator, params.output.format, operators;
                hermitian=params.output.hermitian,
            )

        @info "Writing operators"

        if config_rank == 0
            write_operators(params.output.format, output_dirs[idir], operators)
        end
    end
end

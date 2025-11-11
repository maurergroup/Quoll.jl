module TestUtils

using Tar
using CodecZlib
using Logging
using LoggingExtras
using Dictionaries
using MPI

export create_temptree
export populate_from_dict!
export get_translation_invariant
export setupteardown_tmp, with_nowarn_logger, mpiexec_quollapp

"""
    create_temptree(files)

Create a temporary tree of paths stored in a collection of relative paths
in `files`. Returns the root and the resulting absolute paths of the created files.
"""
function create_temptree(paths)
    root = mktempdir(; cleanup = true)

    abspaths = String[]
    for p in paths
        push!(abspaths, mkpath(joinpath(root, p)))
    end

    return (root, abspaths)
end

# Setup:
# - create temporary directory /tmp/<...>
# - go to that directory
# - extract each tarball in /tmp/<...>/<tarball_name>
# Teardown:
# - move back to the original directory
function setupteardown_tmp(f; tarballs = [])
    starting_dir = pwd()
    tarballs = abspath.(tarballs)
    try
        # Create temporary directory
        tempdir = mktempdir(; cleanup = true)

        # Move to the temporary directory
        cd(tempdir)

        # Extract any supplied tarballs
        for tpath in tarballs
            bname = match(r"^[^\.]+", basename(tpath)).match
            tempdir_ext = joinpath(tempdir, bname)

            tar = GzipDecompressorStream(open(tpath))
            Tar.extract(tar, tempdir_ext)
        end

        # run the test
        f()

    finally
        cd(starting_dir)
    end
end

function with_nowarn_logger(f)
    error_only_logger = MinLevelLogger(current_logger(), Logging.Error)
    return with_logger(f, error_only_logger)
end

function populate_from_dict!(arr::AbstractArray, d::AbstractDictionary)
    setindex!(arr, collect(values(d)), CartesianIndex.(collect(keys(d))))
end

function populate_from_dict!(v::AbstractVector, d::AbstractDictionary)
    setindex!(v, collect(values(d)), collect(keys(d)))
end

# Precompile = true could alleviate race conditions during precompilation
# if there are any (MPI.jl known issue, but if I understand correctly
# it might have been fixed in Base Julia already)
function mpiexec_quollapp(app, project, inputfile; np, nt = 1, precompile = false)
    cmd_list = Cmd[]
    if precompile
        push!(cmd_list, `$(Base.julia_cmd()) --project=$project -e 'using Pkg; Pkg.instantiate()'`)
        push!(cmd_list, `$(Base.julia_cmd()) --project=$project -e 'using Pkg; Pkg.precompile()'`)
    end
    cmd_mpiexec = `
        $(mpiexec()) -n $np
        $(Base.julia_cmd()) -t $nt
        --startup-file=no
        --project=$project
        $app $inputfile
    `
    push!(cmd_list, cmd_mpiexec)
    return cmd_list
end

# Assuming 3 x N matrix
function get_translation_invariant(coords::AbstractMatrix)
    centroid = get_centroid(coords)
    return coords .- centroid
end

# Assuming 3 x N matrix -> 3 x 1 matrix
function get_centroid(coords::AbstractMatrix)
    norm = 1 / size(coords, 2)
    return norm * sum(coords, dims = 2)
end

end
module TestUtils

using Tar
using CodecZlib
using Logging
using LoggingExtras
using Dictionaries
using MPI

export create_temptree
export populate_from_dict!
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

function mpiexec_quollapp(app, project, inputfile; np, nt = 1)
    `$(mpiexec()) -n $np $(Base.julia_cmd()) -t $nt --startup-file=no --project=$project $app $inputfile`
end

end
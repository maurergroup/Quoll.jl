module TestUtils

using Tar
using CodecZlib
using Logging
using LoggingExtras

export create_temptree, setupteardown_tmp, with_nowarn_logger

"""
    create_temptree(files)

Create a temporary tree of paths stored in a collection of relative paths
in `files`. Returns the root and the resulting absolute paths of the created files.
"""
function create_temptree(paths)
    root = mktempdir(; cleanup = true)

    abspaths = []
    for p in paths
        push!(abspaths, mkpath(joinpath(root, p)))
    end

    return (root, abspaths)
end

# Setup:
# - create temporary directory,
# - go to that directory,
# - copy tarballs into that directory
# - extract tarballs
# Teardown:
# - move back to original directory
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

end
using Tar
using CodecZlib

"""
    create_tempdirs(dirs)

Create a temporary tree of directories stored in a collection of
relative paths in `dirs`. Returns the root and the resulting paths
of the created directories.
"""
function create_tempdirs(dirs)
    root = mktempdir(; cleanup = true)

    paths = []
    for dir in dirs
        push!(paths, mkpath(joinpath(root, dir)))
    end

    return (root, paths)
end

# Setup:
# - create temporary directory,
# - go to that directory,
# - copy tarballs into that directory
# - extract tarballs
# Teardown:
# - move back to original directory
function setupteardown_tarballs(f, args, kwargs, tarballs)
    starting_dir = pwd()
    tarballs = abspath.(tarballs)
    try
        # Create temporary directory
        tempdir = mktempdir(; cleanup = true)

        # Move to the temporary directory
        cd(tempdir)

        # Extract tarballs
        for tpath in tarballs
            bname = match(r"^[^\.]+", basename(tpath)).match
            tempdir_ext = joinpath(tempdir, bname)

            tar = GzipDecompressorStream(open(tpath))
            Tar.extract(tar, tempdir_ext)
        end

        # run the function
        f(args...; kwargs...)

    finally
        cd(starting_dir)
    end
end

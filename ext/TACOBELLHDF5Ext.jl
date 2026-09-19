"""
    TACOBELLHDF5Ext

Package extension that implements [`TACOBELL.ingest_jld2!`](@ref) once the
user loads HDF5.jl (`using HDF5`). Kept out of the main module so that
TACOBELL itself has no hard dependency on HDF5.jl.

JLD2 files are HDF5 containers, so a TNRKit `.jld2` output can be read
directly with HDF5.jl without needing JLD2.jl as a dependency.
"""
module TACOBELLHDF5Ext

using TACOBELL
using HDF5

using TACOBELL: RunParameters, CFTResults, ScalingDimSector, insert_run!, generate_catalog

"""
    ingest_jld2!(filepath, params; db_path, allow_duplicate)

See [`TACOBELL.ingest_jld2!`](@ref) for the full docstring. This method
becomes available once `HDF5` is loaded.
"""
function TACOBELL.ingest_jld2!(
        filepath :: String,
        params   :: RunParameters;
        db_path  :: String  = TACOBELL.DB_PATH(),
        allow_duplicate :: Bool = false,
    )
    HDF5.h5open(filepath, "r") do f
        chi    = Int(HDF5.read(f["chi"]))
        K      = Int(HDF5.read(f["K"]))
        mu0_sq = Float64(HDF5.read(f["μ0"]))
        lambda = Float64(HDF5.read(f["λ"]))
        t      = Float64(HDF5.read(f["t"]))

        filled_params = RunParameters(
            model     = params.model,
            symmetry  = params.symmetry,
            algorithm = params.algorithm,
            chi       = chi,
            K         = K,
            mu0_sq    = mu0_sq,
            lambda    = lambda,
        )

        data_refs = HDF5.read(f["data"])   # Vector of object references
        iters = CFTResults[]

        for (i, ref) in enumerate(data_refs)
            step = i - 1   # 0-based iteration index
            entry = f[ref][]

            norm = Float64(entry["1"])

            cc_ref = entry["2"]["central_charge"]
            cc_re  = Float64(HDF5.read(f[cc_ref]["re"]))

            sd     = entry["2"]["scaling_dimensions"]
            sd_re  = Float64.(HDF5.read(f[sd["data"]]["re"]))

            sectors = _parse_sectors(f, sd["structure"], sd_re)

            push!(iters, CFTResults(
                iteration      = step,
                normalization  = norm,
                central_charge = cc_re,
                sectors        = sectors,
            ))
        end

        insert_run!(filled_params, iters;
                    db_path, allow_duplicate,
                    source_file=basename(filepath),
                    runtime_s=t)
    end
end

# Parse the sector structure out of the JLD2 scaling_dimensions.structure field
function _parse_sectors(f, struct_ref, all_dims::Vector{Float64})
    struct_val = f[struct_ref][]
    kvvec_ref  = struct_val["kvvec"]
    kvvec      = f[kvvec_ref][]          # Vector of object references, one per sector

    sectors = ScalingDimSector[]
    for ref in kvvec
        pair_val = f[ref][]              # named tuple: first=sector label, second=indices
        first_ref  = pair_val["first"]
        second_ref = pair_val["second"]

        label     = f[first_ref][]
        j_ref     = label["j"]
        twice_j   = Int(f[j_ref][]["twice"])
        s         = Int(label["s"])

        indices   = Int.(f[second_ref][])   # 1-based into all_dims
        dims      = isempty(indices) ? Float64[] : sort(all_dims[indices])

        push!(sectors, ScalingDimSector(twice_j=twice_j, s=s, dims=dims))
    end
    return sectors
end

"""
    ingest_directory!(dir; infer_params, db_path, extension=".jld2", allow_duplicate=false)

See [`TACOBELL.ingest_directory!`](@ref) for the full docstring.
"""
function TACOBELL.ingest_directory!(
        dir :: String;
        infer_params,
        db_path     :: String  = TACOBELL.DB_PATH(),
        extension   :: String  = ".jld2",
        allow_duplicate :: Bool = false,
    )
    files = String[]
    for (root, _, fnames) in walkdir(dir)
        for fn in fnames
            endswith(lowercase(fn), lowercase(extension)) && push!(files, joinpath(root, fn))
        end
    end
    sort!(files)
    isempty(files) && @warn "No *$(extension) files found under $dir"

    n_ok = n_skip = n_err = 0
    for (i, filepath) in enumerate(files)
        label = "[$i/$(length(files))] $(basename(filepath))"
        try
            params = infer_params(filepath)
            TACOBELL.ingest_jld2!(filepath, params; db_path, allow_duplicate)
            n_ok += 1
            @info "$label -> inserted"
        catch err
            if err isa ErrorException && startswith(err.msg, "Duplicate:")
                n_skip += 1
                @info "$label -> skipped (already in db)"
            else
                n_err += 1
                @warn "$label -> FAILED" exception=(err, catch_backtrace())
            end
        end
    end

    n_ok > 0 && generate_catalog(; db_path)

    summary = (total=length(files), inserted=n_ok, skipped=n_skip, failed=n_err)
    @info "Bulk ingest complete" summary...
    return summary
end

end # module

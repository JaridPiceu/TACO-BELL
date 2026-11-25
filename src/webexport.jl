# ─────────────────────────────────────────────────────────────────────────────
# Explorer webpage export (web/index.html)
#
#   web/data/index.json        all runs, column-wise, strings dictionary-encoded
#   web/data/runs/<id>.bin     one run's full trajectory, zlib-compressed binary
#
# Run files are what the page fetches lazily for the matched run. They hold
# exactly what the page can display: scaling dimensions quantized to 1e-4
# (the page shows 4 decimals), delta-encoded within each sector (sorted, so
# deltas are small) and byte-shuffled before compression — ~10 KB per run
# instead of ~70 KB of JSON. Layout after inflating (little-endian):
#
#   "TBW1"  u32 nIter  u32 nSec  u32 nDims  u32 labelsLen
#   labels  (labelsLen bytes of UTF-8 JSON: array of charge objects)
#   i32[nIter] iteration     f64[nIter] central charge (NaN = missing)
#   u16[nIter] sectors per iteration
#   u16[nSec]  label index   u16[nSec]  dims per sector
#   i32[nDims] round(Δ·1e4), delta-coded within each sector, stored
#              byte-shuffled (all 1st bytes, then all 2nd bytes, …)
#
# The decoder is `decodeRun` in web/index.html.
# ─────────────────────────────────────────────────────────────────────────────

const WEB_DIMS_SCALE = 10_000

function _zlib_compress(data::Vector{UInt8}; level::Integer = 9)
    bound = ccall((:compressBound, Zlib_jll.libz), Culong, (Culong,), length(data))
    out = Vector{UInt8}(undef, bound)
    outlen = Ref{Culong}(bound)
    rc = ccall((:compress2, Zlib_jll.libz), Cint,
               (Ptr{UInt8}, Ref{Culong}, Ptr{UInt8}, Culong, Cint),
               out, outlen, data, length(data), level)
    rc == 0 || error("zlib compress2 failed (code $rc)")
    return resize!(out, outlen[])
end

function _zlib_decompress(data::Vector{UInt8}, maxlen::Integer)
    out = Vector{UInt8}(undef, maxlen)
    outlen = Ref{Culong}(maxlen)
    rc = ccall((:uncompress, Zlib_jll.libz), Cint,
               (Ptr{UInt8}, Ref{Culong}, Ptr{UInt8}, Culong),
               out, outlen, data, length(data))
    rc == 0 || error("zlib uncompress failed (code $rc)")
    return resize!(out, outlen[])
end

function _web_run_bytes(e::DatabaseEntry)
    iters = sort(e.iterations; by = r -> r.iteration)
    labels = Dict{String,Float64}[]
    label_idx = Dict{Dict{String,Float64},Int}()
    nsec_per_iter = UInt16[]; sec_label = UInt16[]; sec_n = UInt16[]; q = Int32[]
    for r in iters
        push!(nsec_per_iter, length(r.sectors))
        for s in r.sectors
            idx = get!(label_idx, s.charge) do
                push!(labels, s.charge)
                length(labels)
            end
            push!(sec_label, idx - 1)
            push!(sec_n, length(s.dims))
            prev = Int32(0)
            for d in sort(s.dims)
                v = round(Int32, d * WEB_DIMS_SCALE)
                push!(q, v - prev)
                prev = v
            end
        end
    end
    labels_json = Vector{UInt8}(_to_json(labels))

    io = IOBuffer()
    write(io, b"TBW1")
    for n in (length(iters), length(sec_label), length(q), length(labels_json))
        write(io, htol(UInt32(n)))
    end
    write(io, labels_json)
    foreach(r -> write(io, htol(Int32(r.iteration))), iters)
    foreach(r -> write(io, htol(_nan(r.central_charge))), iters)
    foreach(x -> write(io, htol(x)), nsec_per_iter)
    foreach(x -> write(io, htol(x)), sec_label)
    foreach(x -> write(io, htol(x)), sec_n)
    raw = reinterpret(UInt8, htol.(q))
    write(io, vec(permutedims(reshape(raw, 4, :))))   # byte-shuffle
    return _zlib_compress(take!(io))
end

"""
    export_web(; db_path, dir="web/data") -> String

(Re)build everything the Explorer webpage (`web/index.html`) reads from the
whole database: `dir/index.json` (every run's parameters and plateau central
charge, stored column-wise) and one compact `dir/runs/<id>.bin` per run with
its full trajectory. Stale `.bin`/`.json` files already in `dir/runs` are
removed first. Works shard by shard, so memory use stays bounded no matter
how large the database grows. Returns `dir`.

Run it via `julia --project=. scripts/update_web.jl`; the GitHub Pages
workflow does exactly that on every push that changes the database.
"""
function export_web(; db_path::String = DB_PATH(), dir::String = joinpath("web", "data"))
    rundir = joinpath(dir, "runs")
    mkpath(rundir)
    for f in readdir(rundir; join = true)
        (endswith(f, ".bin") || endswith(f, ".json")) && rm(f)
    end

    rows = _load_index(db_path)
    # 12 hex digits of the UUID is plenty to be unique (checked below) and
    # keeps the index small; fall back to the full id if it ever isn't.
    short = [first(replace(r["id"], "-" => ""), 12) for r in rows]
    allunique(short) || (short = [replace(r["id"], "-" => "") for r in rows])
    short_of = Dict(r["id"] => s for (r, s) in zip(rows, short))

    models = sort(unique(r["model"] for r in rows))
    syms   = sort(unique(r["symmetry"] for r in rows))
    algs   = sort(unique(r["algorithm"] for r in rows))
    code(v, x) = findfirst(==(x), v) - 1
    col(k) = [r[k] for r in rows]

    index = Dict{String,Any}(
        "format" => 1,
        "generated_at" => string(now()),
        "models" => models, "symmetries" => syms, "algorithms" => algs,
        "id" => short,
        "model" => [code(models, r["model"]) for r in rows],
        "symmetry" => [code(syms, r["symmetry"]) for r in rows],
        "algorithm" => [code(algs, r["algorithm"]) for r in rows],
        "chi" => col("chi"), "K" => col("K"),
        "mu0_sq" => col("mu0_sq"), "lambda" => col("lambda"),
        "n_iterations" => col("n_iterations"),
        "central_charge" => [isnan(r["central_charge"]) ? missing : r["central_charge"] for r in rows],
    )
    write(joinpath(dir, "index.json"), _to_json(index))

    for path in _shard_files(db_path)
        for e in _read_shard(path)
            write(joinpath(rundir, short_of[e.id] * ".bin"), _web_run_bytes(e))
        end
    end
    @info "Wrote web index ($(length(rows)) runs) and run files to $dir"
    return dir
end

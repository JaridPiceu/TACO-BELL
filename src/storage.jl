# ─────────────────────────────────────────────────────────────────────────────
# On-disk storage: append-only, columnar HDF5 shards
#
#   db/shards/<yyyymmdd-HHMMSS>-<rand>.h5
#
# Every ingest writes one or more *new* shard files and never touches the
# existing ones, so git stores each byte of data exactly once (no rewritten
# index or catalog blobs piling up in history), and adding a batch of runs
# is a single bulk write instead of one file + one full index rewrite per run.
#
# Inside a shard, everything is stored column-wise with "counts" columns that
# link the levels together (CSR-style), which is both compact and trivially
# readable from Python (h5py), MATLAB, R, … without any TACOBELL code:
#
#   /runs/        id, created_at, source_file, model, symmetry, algorithm (strings)
#                 chi, K, n_iterations (Int32)   mu0_sq, lambda, runtime_s (Float64)
#   /iterations/  iteration, n_sectors (Int32)   normalization, central_charge (Float64)
#                 notes_index (Int64, 0-based) + notes_text — only non-empty notes
#   /sectors/     charge (UInt16 index into charge_labels), n_dims (Int32)
#                 charge_labels: "j=1.0;s=2.0"-style strings ("" = trivial sector)
#   /dims         every scaling dimension, flat (Float32 by default, see DIMS_TYPE)
#
# Run r owns the next `n_iterations[r]` iteration rows, iteration i owns the
# next `n_sectors[i]` sector rows, sector s owns the next `n_dims[s]` values
# of /dims. `missing` is stored as NaN. All numeric columns are chunked and
# compressed (shuffle + deflate), which every HDF5 reader supports natively.
# ─────────────────────────────────────────────────────────────────────────────

const SHARD_FORMAT = 1

"""
    TACOBELL.DIMS_TYPE[]

Element type used to store scaling dimensions in newly written shards —
`Float32` by default (7 significant digits; the relative rounding error, ≤6e-8, is
orders of magnitude below the accuracy of any TNR calculation, and it makes
the database ~2.5× smaller than storing full `Float64`). Set
`TACOBELL.DIMS_TYPE[] = Float64` before ingesting to store them losslessly
instead. Each shard records its own type, so both kinds can coexist; reading
always returns `Float64`. Central charges, normalizations and parameters are
always stored as full `Float64`.
"""
const DIMS_TYPE = Ref{DataType}(Float32)

"""
    TACOBELL.SHARD_MAX_DIMS[]

Soft cap on the number of scaling dimensions per shard file (default 12
million ≈ 30–35 MB on disk), keeping every file well under GitHub's 100 MB
per-file limit. A batch larger than this is split over several shards.
"""
const SHARD_MAX_DIMS = Ref(12_000_000)

_shard_dir(db_path) = joinpath(db_path, "shards")

function _shard_files(db_path)
    d = _shard_dir(db_path)
    isdir(d) || return String[]
    return sort!(filter(f -> endswith(f, ".h5"), readdir(d; join = true)))
end

_new_shard_path(db_path) = joinpath(_shard_dir(db_path),
    Dates.format(now(), "yyyymmdd-HHMMSS") * "-" * first(string(uuid4()), 8) * ".h5")

_nan(x) = x === missing ? NaN : Float64(x)
_unnan(x::AbstractFloat) = isnan(x) ? missing : Float64(x)

# Charge dicts ⇄ compact labels. Keys are TensorKit field names (`j`, `s`,
# `n`, `charge`, …), so they never contain the separators.
_charge_key(c::Dict{String,Float64}) =
    join(("$k=$(repr(v))" for (k, v) in sort!(collect(c); by = first)), ";")

function _charge_from_key(s::AbstractString)
    d = Dict{String,Float64}()
    isempty(s) && return d
    for kv in split(s, ';')
        k, v = split(kv, '='; limit = 2)
        d[String(k)] = parse(Float64, v)
    end
    return d
end

function _wcol(parent, name, v::AbstractVector{<:AbstractString})
    if isempty(v)
        create_dataset(parent, name, datatype(String), dataspace((0,)))
    else
        parent[name] = String.(v)
    end
end

function _wcol(parent, name, v::AbstractVector{T}) where {T<:Real}
    n = length(v)
    ds = create_dataset(parent, name, datatype(T), dataspace((n,));
                        chunk = (clamp(n, 1, 1 << 16),), shuffle = true, deflate = 6)
    n > 0 && write(ds, collect(v))
    return ds
end

_rcol(parent, name) = (ds = parent[name]; length(ds) == 0 ? eltype(ds)[] : read(ds))
_rcol_str(parent, name) = (ds = parent[name]; length(ds) == 0 ? String[] : String.(read(ds)))

_ndims_total(e::DatabaseEntry) = sum((length(s.dims) for r in e.iterations for s in r.sectors); init = 0)

"""
Write `entries` to a new shard file at `path` (atomically: written to a
temporary file first, then moved into place).
"""
function _write_shard(path::AbstractString, entries::AbstractVector{DatabaseEntry};
                      dims_type::DataType = DIMS_TYPE[])
    labels = String[]
    label_idx = Dict{String,Int}()

    n_iterations = Int32[]
    it_iteration = Int32[]; it_norm = Float64[]; it_cc = Float64[]; it_nsec = Int32[]
    notes_index = Int64[]; notes_text = String[]
    sec_charge = UInt16[]; sec_ndims = Int32[]
    dims = dims_type[]
    sizehint!(dims, sum(_ndims_total, entries; init = 0))

    for e in entries
        push!(n_iterations, length(e.iterations))
        for r in e.iterations
            isempty(r.notes) || (push!(notes_index, length(it_iteration)); push!(notes_text, r.notes))
            push!(it_iteration, r.iteration)
            push!(it_norm, _nan(r.normalization))
            push!(it_cc, _nan(r.central_charge))
            push!(it_nsec, length(r.sectors))
            for s in r.sectors
                key = _charge_key(s.charge)
                idx = get!(label_idx, key) do
                    push!(labels, key)
                    length(labels)
                end
                idx <= typemax(UInt16) || error("more than $(typemax(UInt16)) distinct charge labels in one shard")
                push!(sec_charge, UInt16(idx - 1))
                push!(sec_ndims, length(s.dims))
                append!(dims, s.dims)
            end
        end
    end

    mkpath(dirname(path))
    tmp = path * ".tmp"
    h5open(tmp, "w") do f
        attributes(f)["tacobell_format"] = SHARD_FORMAT
        attributes(f)["written_at"] = string(now())

        g = create_group(f, "runs")
        _wcol(g, "id",          [e.id for e in entries])
        _wcol(g, "created_at",  [e.created_at for e in entries])
        _wcol(g, "source_file", [e.source_file for e in entries])
        _wcol(g, "model",       [e.params.model for e in entries])
        _wcol(g, "symmetry",    [e.params.symmetry for e in entries])
        _wcol(g, "algorithm",   [e.params.algorithm for e in entries])
        _wcol(g, "chi",         Int32[e.params.chi for e in entries])
        _wcol(g, "K",           Int32[e.params.K for e in entries])
        _wcol(g, "mu0_sq",      Float64[e.params.mu0_sq for e in entries])
        _wcol(g, "lambda",      Float64[e.params.lambda for e in entries])
        _wcol(g, "runtime_s",   Float64[_nan(e.runtime_s) for e in entries])
        _wcol(g, "n_iterations", n_iterations)

        g = create_group(f, "iterations")
        _wcol(g, "iteration",      it_iteration)
        _wcol(g, "normalization",  it_norm)
        _wcol(g, "central_charge", it_cc)
        _wcol(g, "n_sectors",      it_nsec)
        _wcol(g, "notes_index",    notes_index)
        _wcol(g, "notes_text",     notes_text)

        g = create_group(f, "sectors")
        _wcol(g, "charge",        sec_charge)
        _wcol(g, "n_dims",        sec_ndims)
        _wcol(g, "charge_labels", labels)

        _wcol(f, "dims", dims)
    end
    mv(tmp, path; force = true)
    lock(_CACHE_LOCK) do                 # in case an existing shard was rewritten
        delete!(_META_CACHE, path)
        delete!(_DETAIL_CACHE, path)
    end
    return path
end

# ─── Cached readers ──────────────────────────────────────────────────────────
#
# Shards are immutable in normal use, so their decoded metadata is cached in
# memory keyed by (path, mtime, size): the index of every run in the database
# is rebuilt from a handful of small column reads the first time it's needed
# and then reused, and anything rewritten on disk is picked up automatically.

_stamp(path) = (s = stat(path); (s.mtime, s.size))

struct ShardMeta
    stamp      :: Tuple{Float64,Int}
    rows       :: Vector{Dict{String,Any}}   # one index row per run
    iter_start :: Vector{Int}                # run r ↦ iterations iter_start[r]+1 : iter_start[r+1]
end

struct ShardDetail
    stamp          :: Tuple{Float64,Int}
    iteration      :: Vector{Int32}
    normalization  :: Vector{Float64}
    central_charge :: Vector{Float64}
    sec_start      :: Vector{Int}            # iteration i ↦ sectors sec_start[i]+1 : sec_start[i+1]
    charge         :: Vector{UInt16}
    labels         :: Vector{Dict{String,Float64}}
    dim_start      :: Vector{Int}            # sector s ↦ dims dim_start[s]+1 : dim_start[s+1]
    notes          :: Dict{Int,String}       # 1-based iteration row ↦ note
end

const _META_CACHE   = Dict{String,ShardMeta}()
const _DETAIL_CACHE = Dict{String,ShardDetail}()
const _CACHE_LOCK   = ReentrantLock()

_offsets(counts) = (o = zeros(Int, length(counts) + 1); cumsum!(view(o, 2:length(o)), counts); o)

function _shard_meta(path::String)
    st = _stamp(path)
    cached = lock(() -> get(_META_CACHE, path, nothing), _CACHE_LOCK)
    cached !== nothing && cached.stamp == st && return cached

    meta = h5open(path, "r") do f
        g = f["runs"]
        ids   = _rcol_str(g, "id")
        cat   = _rcol_str(g, "created_at")
        model = _rcol_str(g, "model"); sym = _rcol_str(g, "symmetry"); alg = _rcol_str(g, "algorithm")
        chi   = _rcol(g, "chi"); K = _rcol(g, "K")
        mu0   = _rcol(g, "mu0_sq"); lam = _rcol(g, "lambda")
        niter = _rcol(g, "n_iterations")
        it    = f["iterations"]
        iters = _rcol(it, "iteration"); cc = _rcol(it, "central_charge")
        istart = _offsets(niter)

        rows = Vector{Dict{String,Any}}(undef, length(ids))
        for r in eachindex(ids)
            rng = (istart[r] + 1):istart[r + 1]
            traj = Float64[cc[i] for i in rng[sortperm(iters[rng])] if !isnan(cc[i])]
            rows[r] = Dict{String,Any}(
                "id" => ids[r], "created_at" => cat[r],
                "model" => model[r], "symmetry" => sym[r], "algorithm" => alg[r],
                "chi" => Int(chi[r]), "K" => Int(K[r]),
                "mu0_sq" => mu0[r], "lambda" => lam[r],
                "n_iterations" => Int(niter[r]),
                "central_charge" => isempty(traj) ? NaN : plateau_estimate(traj).value,
                "_shard" => path, "_row" => r,
            )
        end
        ShardMeta(st, rows, istart)
    end
    lock(() -> (_META_CACHE[path] = meta), _CACHE_LOCK)
    return meta
end

function _shard_detail(path::String)
    st = _stamp(path)
    cached = lock(() -> get(_DETAIL_CACHE, path, nothing), _CACHE_LOCK)
    cached !== nothing && cached.stamp == st && return cached

    det = h5open(path, "r") do f
        it = f["iterations"]; sc = f["sectors"]
        nidx = _rcol(it, "notes_index"); ntxt = _rcol_str(it, "notes_text")
        ShardDetail(st,
            _rcol(it, "iteration"), _rcol(it, "normalization"), _rcol(it, "central_charge"),
            _offsets(_rcol(it, "n_sectors")),
            _rcol(sc, "charge"),
            [_charge_from_key(s) for s in _rcol_str(sc, "charge_labels")],
            _offsets(_rcol(sc, "n_dims")),
            Dict{Int,String}(Int(i) + 1 => t for (i, t) in zip(nidx, ntxt)),
        )
    end
    lock(() -> (_DETAIL_CACHE[path] = det), _CACHE_LOCK)
    return det
end

"""
Read full `DatabaseEntry`s for the given 1-based run rows of one shard. Few
rows → only their slices of `/dims` are read; many → the whole column once.
"""
function _read_entries(path::String, runs::AbstractVector{<:Integer})
    meta = _shard_meta(path)
    det  = _shard_detail(path)
    isempty(runs) && return DatabaseEntry[]
    h5open(path, "r") do f
        dsd = f["dims"]
        whole = length(runs) > 16 ? Float64.(_rcol(f, "dims")) : nothing
        srcs = _rcol_str(f["runs"], "source_file")
        rts  = _rcol(f["runs"], "runtime_s")
        map(runs) do r
            row = meta.rows[r]
            irng = (meta.iter_start[r] + 1):meta.iter_start[r + 1]
            # this run's contiguous slice of /dims
            d0 = isempty(irng) ? 0 : det.dim_start[det.sec_start[first(irng)] + 1]
            d1 = isempty(irng) ? 0 : det.dim_start[det.sec_start[last(irng) + 1] + 1]
            vals = whole !== nothing ? view(whole, (d0 + 1):d1) :
                   d1 > d0 ? Float64.(dsd[(d0 + 1):d1]) : Float64[]
            iters = map(irng) do i
                srng = (det.sec_start[i] + 1):det.sec_start[i + 1]
                secs = map(srng) do s
                    ScalingDimSector(
                        # shared per label rather than copied per sector:
                        # millions of sectors, only a handful of distinct labels
                        charge = det.labels[det.charge[s] + 1],
                        dims   = collect(vals[(det.dim_start[s] - d0 + 1):(det.dim_start[s + 1] - d0)]),
                    )
                end
                CFTResults(
                    iteration      = Int(det.iteration[i]),
                    normalization  = _unnan(det.normalization[i]),
                    central_charge = _unnan(det.central_charge[i]),
                    sectors        = secs,
                    notes          = get(det.notes, i, ""),
                )
            end
            p = RunParameters(
                model = row["model"], symmetry = row["symmetry"], algorithm = row["algorithm"],
                chi = row["chi"], K = row["K"], mu0_sq = row["mu0_sq"], lambda = row["lambda"],
            )
            DatabaseEntry(id = row["id"], created_at = row["created_at"], source_file = srcs[r],
                          runtime_s = _unnan(rts[r]), params = p, iterations = collect(CFTResults, iters))
        end
    end
end

_shard_dims_type(path) = h5open(f -> eltype(f["dims"]), path, "r")

"""All entries stored in one shard, in stored order."""
_read_shard(path::String) = _read_entries(path, 1:length(_shard_meta(path).rows))

# ─── Database-level index ────────────────────────────────────────────────────

"""
    _load_index(db_path) -> Vector{Dict{String,Any}}

One summary row per run across all shards (parameters, iteration count and
plateau central charge). Derived from the shards on the fly and cached — it
is never written to disk, so it can't go stale or cause merge conflicts.
"""
function _load_index(db_path = DB_PATH())
    rows = Dict{String,Any}[]
    for path in _shard_files(db_path)
        append!(rows, _shard_meta(path).rows)
    end
    return rows
end

_param_key(p::RunParameters) =
    (p.model, p.symmetry, p.algorithm, p.chi, p.K, round(p.mu0_sq; digits = 10), round(p.lambda; digits = 10))
_param_key(row::AbstractDict) =
    (row["model"], row["symmetry"], row["algorithm"], row["chi"], row["K"],
     round(row["mu0_sq"]; digits = 10), round(row["lambda"]; digits = 10))

"""Map from parameter key to the id of an existing run with those parameters."""
_existing_keys(db_path) = Dict(_param_key(r) => r["id"] for r in _load_index(db_path))

"""Load the full entries behind a set of index rows, grouped per shard."""
function _entries_for_rows(rows)
    isempty(rows) && return DatabaseEntry[]
    byshard = Dict{String,Vector{Int}}()
    order = Dict{String,Int}()
    for (k, row) in enumerate(rows)
        push!(get!(byshard, row["_shard"], Int[]), row["_row"])
        order[row["id"]] = k
    end
    out = Vector{DatabaseEntry}(undef, length(rows))
    for (path, rs) in byshard
        for e in _read_entries(path, rs)
            out[order[e.id]] = e
        end
    end
    return out
end

"""Split `entries` into consecutive groups of at most ~`SHARD_MAX_DIMS[]` dims."""
function _chunk_by_dims(entries::AbstractVector{DatabaseEntry}; max_dims::Int = SHARD_MAX_DIMS[])
    chunks = Vector{DatabaseEntry}[]
    cur = DatabaseEntry[]; n = 0
    for e in entries
        k = _ndims_total(e)
        if !isempty(cur) && n + k > max_dims
            push!(chunks, cur); cur = DatabaseEntry[]; n = 0
        end
        push!(cur, e); n += k
    end
    isempty(cur) || push!(chunks, cur)
    return chunks
end

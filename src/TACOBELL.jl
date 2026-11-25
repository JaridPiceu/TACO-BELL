"""
    TACOBELL

Tensor Archive of Conformal Output — Best Ever Lattice Labour.

A lightweight, file-based database for CFT data (central charges, scaling
dimensions, ...) extracted from Tensor Network Renormalization (TNR)
calculations, designed for use alongside
[TNRKit](https://github.com/QuantumKitHub/TNRKit.jl/).

See the package README for a full walkthrough. The core entry points are
[`add_data!`](@ref) (add TNRKit output: a file or a whole folder; or
[`insert_run!`](@ref) for results built by hand),[`query_runs`](@ref)/[`find_closest`](@ref) (look it up), and
[`export_csv`](@ref)/[`export_web`](@ref)/[`generate_catalog`](@ref) (get it
out again, as a CSV file, as data for the Explorer webpage, or as a Markdown
summary for browsing on GitHub).

Data lives in append-only, compressed, columnar HDF5 shard files under
`db/shards/` — see `src/storage.jl` for the exact layout.
"""
module TACOBELL

using TOML, UUIDs, Dates, Printf, Statistics
using HDF5
using Plots
import Zlib_jll

export RunParameters, ScalingDimSector, CFTResults, DatabaseEntry
export insert_run!, insert_runs!, query_runs, load_entry, get_iteration, final_iteration, summarize_db
export find_closest, find_sector, list_algorithms, list_symmetries, rebuild_index!, correct_field!
export ingest_jld2!, ingest_directory!, generate_catalog, export_csv, export_json, export_json_runs
export export_web, compact!, import_toml_runs!
export ingest_jld2_old!, ingest_directory_old!
export add_data!, MissingParameterError
export central_charge_trajectory, plateau_estimate, plateau_central_charge
export plot_central_charge

# ─────────────────────────────────────────────────────────────────────────────
# Types
# ─────────────────────────────────────────────────────────────────────────────

"""
    RunParameters

All input parameters that uniquely identify a TNR calculation run.

Fields
------
- `model`      : e.g. `"phi4_real"`, `"phi4_complex"`
- `symmetry`   : manifest symmetry in the tensor, e.g. `"O(2)"`, `"Z2"`, `"U(1)"`, `"none"`
- `algorithm`  : e.g. `"LoopTNR"`, `"BTRG"`, `"GILT-TNR"`
- `chi`        : bond dimension χ
- `K`          : initial truncation
- `mu0_sq`     : bare mass squared μ₀²
- `lambda`     : quartic coupling λ
"""
Base.@kwdef struct RunParameters
    model     :: String
    symmetry  :: String
    algorithm :: String
    chi       :: Int
    K         :: Int
    mu0_sq    :: Float64
    lambda    :: Float64
end

"""
    ScalingDimSector

The scaling dimensions within one symmetry sector, i.e. one fusion-tree
charge/irrep. Storage is **symmetry-agnostic**: `charge` holds whatever
quantum number(s) TensorKit's own charge/irrep type carries for that
symmetry, keyed by the same field names TensorKit itself uses for that
irrep type — so this works unchanged for `Trivial`, `U(1)`, `Z_N`, `SU(2)`,
`O(2)`, or anything else TensorKit can label a sector with.

- `charge` : `Dict{String,Float64}` of the sector's quantum number(s), e.g.
             `Dict("charge"=>1.0)` for `U1Irrep`, `Dict("n"=>1.0)` for a
             `ZNIrrep`, `Dict("j"=>0.5)` for `SU2Irrep`, or
             `Dict("j"=>1.0, "s"=>2.0)` for `CU1Irrep` (= O(2)). Empty for
             the trivial (no-symmetry) sector.
- `dims`   : vector of scaling dimensions Δ in this sector, sorted ascending

Values that are half-integers (as `j` commonly is) are stored as their
actual numeric value (e.g. `0.5`), not doubled — use [`find_sector`](@ref)
with the real value.
"""
Base.@kwdef struct ScalingDimSector
    charge :: Dict{String, Float64} = Dict{String, Float64}()
    dims   :: Vector{Float64}
end

"""
    CFTResults

CFT data extracted from **one iteration step** of a TNR run.

Fields
------
- `iteration`       : RG step index (0 = initial tensor)
- `normalization`   : log of the per-site tensor norm at this step (not
                      recorded in pre-TNRKit-v0.5 output — `missing` there)
- `central_charge`  : extracted central charge `c`
- `sectors`         : scaling dimensions grouped by symmetry sector
- `notes`           : free-text annotation

Also exposes a derived `.dims` property: all scaling dimensions across every
sector, flattened into one sorted `Vector{Float64}` — e.g. `result.dims[2]`
for "the first excited state overall, regardless of which symmetry sector
it's in", matching the flat `cft.scaling_dimensions[i]` indexing TNRKit
itself provides. Use `.sectors` instead when you need to know *which*
sector a given Δ came from.
"""
Base.@kwdef struct CFTResults
    iteration       :: Int
    normalization   :: Union{Float64, Missing}              = missing
    central_charge  :: Union{Float64, Missing}              = missing
    sectors         :: Vector{ScalingDimSector}             = ScalingDimSector[]
    notes           :: String                               = ""
end

function Base.getproperty(r::CFTResults, name::Symbol)
    if name === :dims
        sectors = getfield(r, :sectors)
        return isempty(sectors) ? Float64[] : sort(vcat((sec.dims for sec in sectors)...))
    end
    return getfield(r, name)
end

Base.propertynames(::CFTResults) = (fieldnames(CFTResults)..., :dims)

"""
    DatabaseEntry

One record: a parameter set, results for every stored iteration, and metadata.
"""
Base.@kwdef struct DatabaseEntry
    id           :: String             = string(uuid4())
    created_at   :: String             = string(now())
    source_file  :: String             = ""          # original JLD2 path if ingested
    runtime_s    :: Union{Float64,Missing} = missing  # wall-clock seconds stored in JLD2
    params       :: RunParameters
    iterations   :: Vector{CFTResults}  # one per stored RG step
end

# ─────────────────────────────────────────────────────────────────────────────
# Pretty-printing — a run's raw data is deeply nested and can be hundreds of
# KB (dozens of symmetry sectors × tens of iterations), so the default
# struct dump is unreadable. These give one-line summaries instead; the full
# data underneath is unchanged and still fully accessible via the fields.
# ─────────────────────────────────────────────────────────────────────────────

function Base.show(io::IO, p::RunParameters)
    print(io, "RunParameters(", p.model, ", ", p.symmetry, ", ", p.algorithm,
          ", χ=", p.chi, ", K=", p.K, ", μ₀²=", p.mu0_sq, ", λ=", p.lambda, ")")
end

# Print a half-integer charge as "1/2" rather than "0.5", matching how
# physicists write spins/charges by hand; falls back to the plain number
# for anything that isn't an integer or half-integer.
function _fmt_charge(v::Real)
    isinteger(v) && return string(Int(v))
    twice = 2v
    isinteger(twice) && return string(Int(round(twice)), "/2")
    return string(v)
end

_charge_str(charge::Dict{String,Float64}) =
    isempty(charge) ? "trivial" :
        join(("$k=$(_fmt_charge(v))" for (k, v) in sort(collect(charge))), ", ")

function Base.show(io::IO, s::ScalingDimSector)
    if isempty(s.dims)
        print(io, "ScalingDimSector(", _charge_str(s.charge), ", 0 dims)")
    else
        print(io, "ScalingDimSector(", _charge_str(s.charge), ", ", length(s.dims),
              " dims, Δ∈[", @sprintf("%.4f", first(s.dims)), ", ", @sprintf("%.4f", last(s.dims)), "])")
    end
end

function Base.show(io::IO, r::CFTResults)
    cc = r.central_charge === missing ? "—" : @sprintf("%.4f", r.central_charge)
    norm = r.normalization === missing ? "—" : @sprintf("%.4f", r.normalization)
    populated = count(s -> !isempty(s.dims), r.sectors)
    print(io, "CFTResults(iteration=", r.iteration, ", norm=", norm,
          ", c=", cc, ", ", populated, "/", length(r.sectors), " sectors populated)")
end

function Base.show(io::IO, e::DatabaseEntry)
    p = e.params
    cc = _best_c(e.iterations)
    cc_str = isnan(cc) ? "—" : @sprintf("%.4f", cc)
    print(io, "DatabaseEntry(", first(e.id, 8), "…, ", p.model, "/", p.symmetry, "/", p.algorithm,
          ", χ=", p.chi, ", K=", p.K, ", μ₀²=", p.mu0_sq, ", λ=", p.lambda,
          ", ", length(e.iterations), " iters, c=", cc_str, ")")
end

# ─────────────────────────────────────────────────────────────────────────────
# Storage (see src/storage.jl for the on-disk shard format)
# ─────────────────────────────────────────────────────────────────────────────

DB_PATH() = joinpath(@__DIR__, "..", "db")

include("storage.jl")

# ─────────────────────────────────────────────────────────────────────────────
# Legacy TOML readers — the database used to be one TOML file per run
# (db/runs/<uuid>.toml). Kept only so old run files can still be imported
# with `import_toml_runs!`.
# ─────────────────────────────────────────────────────────────────────────────

function _sector_from_dict(d)
    if haskey(d, "charge")
        charge = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in d["charge"])
    else
        # Pre-generalization files only ever stored O(2) sectors as (twice_j, s).
        charge = Dict{String,Float64}("j" => d["twice_j"] / 2, "s" => Float64(d["s"]))
    end
    ScalingDimSector(charge=charge, dims=Float64.(d["dims"]))
end

function _results_from_dict(d)
    CFTResults(
        iteration       = d["iteration"],
        normalization   = get(d, "normalization",   missing),
        central_charge  = get(d, "central_charge",  missing),
        sectors         = [_sector_from_dict(s) for s in get(d, "sectors", [])],
        notes           = get(d, "notes", ""),
    )
end

function _params_from_dict(d)
    RunParameters(
        model     = d["model"],
        symmetry  = get(d, "symmetry", "unknown"),
        algorithm = d["algorithm"],
        chi       = d["chi"],
        K         = d["K"],
        mu0_sq    = d["mu0_sq"],
        lambda    = d["lambda"],
    )
end

function _entry_from_dict(d)
    DatabaseEntry(
        id          = d["id"],
        created_at  = d["created_at"],
        source_file = get(d, "source_file", ""),
        runtime_s   = get(d, "runtime_s", missing),
        params      = _params_from_dict(d["params"]),
        iterations  = [_results_from_dict(r) for r in d["iterations"]],
    )
end

function _best_c(iters::Vector{CFTResults})
    # The literal last iteration can be numerically unstable (finite-χ
    # truncation error compounding under the RG flow), so use the flattest
    # plateau rather than trusting it blindly — see plateau_estimate.
    traj = Float64[r.central_charge for r in sort(iters, by = r -> r.iteration) if r.central_charge !== missing]
    isempty(traj) && return NaN
    return plateau_estimate(traj).value
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API — writing
# ─────────────────────────────────────────────────────────────────────────────

"""
    insert_runs!(entries::Vector{DatabaseEntry}; db_path, allow_duplicate=false)

Insert many entries at once, written together as new shard file(s) under
`db/shards/` (split over several files if larger than
[`TACOBELL.SHARD_MAX_DIMS`](@ref)). This is the fast path for bulk data —
one write per batch rather than one per run. Existing shards are never
modified. Returns `entries`.

Raises an error (and writes nothing) if any entry has the same parameter set
as a run already in the database or earlier in `entries`, unless
`allow_duplicate=true`.
"""
function insert_runs!(
        entries :: Vector{DatabaseEntry};
        db_path :: String = DB_PATH(),
        allow_duplicate :: Bool = false,
    )
    isempty(entries) && return entries
    if !allow_duplicate
        seen = _existing_keys(db_path)
        for e in entries
            k = _param_key(e.params)
            haskey(seen, k) && error(
                "Duplicate: entry $(seen[k]) has identical parameters to $(e.params). " *
                "Pass allow_duplicate=true to insert anyway.")
            seen[k] = e.id
        end
    end
    for chunk in _chunk_by_dims(entries)
        _write_shard(_new_shard_path(db_path), chunk)
    end
    return entries
end

"""
    insert_run!(params, iterations; db_path, allow_duplicate, source_file, runtime_s)

Insert a new entry (parameter set + all iteration results) into the database.
Returns the created `DatabaseEntry`.

Raises an error if an identical parameter set is already in the database,
unless `allow_duplicate=true`.

Each call writes one small shard file; when adding many runs, collect them
and use [`insert_runs!`](@ref) (or [`ingest_directory!`](@ref)) instead, or
merge small shards afterwards with [`compact!`](@ref).
"""
function insert_run!(
        params      :: RunParameters,
        iterations  :: Vector{CFTResults};
        db_path     :: String             = DB_PATH(),
        allow_duplicate :: Bool           = false,
        source_file :: String             = "",
        runtime_s   :: Union{Float64,Missing} = missing,
    )
    entry = DatabaseEntry(
        params       = params,
        iterations   = iterations,
        source_file  = source_file,
        runtime_s    = runtime_s,
    )
    insert_runs!([entry]; db_path, allow_duplicate)
    @info "Inserted run $(entry.id) ($(length(iterations)) iterations)"
    return entry
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API — JLD2 ingestion
# ─────────────────────────────────────────────────────────────────────────────

"""
    ingest_jld2!(filepath, model; symmetry="", algorithm="", db_path, allow_duplicate)
    ingest_jld2!(filepath, params::RunParameters; db_path, allow_duplicate)

Read a JLD2 output file produced by TNRKit and insert it into the database.

The simple form only needs `model` (`"phi4_real"` / `"phi4_complex"`) — the
one field that's never stored inside the JLD2:

```julia
using TACOBELL
ingest_jld2!("Com_PD_O2_mu0-2_0_lam1_0_K8_chi16_iter20.jld2", "phi4_complex")
```

Newer TNRKit output also stores `symmetry`/`algorithm` directly, in which
case those are read from the file; pass them as keywords only as a fallback
for older files that predate this. Everything else (χ, K, μ₀², λ,
iterations, normalization, central charge, scaling dimensions) is always
read from the file, regardless of which symmetry it used.

The `RunParameters` form does the same thing but takes a full struct
instead — mainly useful when you already have one on hand, e.g. from
[`ingest_directory!`](@ref)'s `infer_params`.

JLD2 files are HDF5 containers, so this reads them directly with HDF5.jl
rather than depending on JLD2.jl.
"""
function ingest_jld2!(
        filepath :: String,
        model    :: AbstractString;
        symmetry :: AbstractString = "",
        algorithm :: AbstractString = "",
        db_path  :: String  = DB_PATH(),
        allow_duplicate :: Bool = false,
    )
    params = _blank_params(model; symmetry, algorithm)
    ingest_jld2!(filepath, params; db_path, allow_duplicate)
end

function ingest_jld2!(
        filepath :: String,
        params   :: RunParameters;
        db_path  :: String  = DB_PATH(),
        allow_duplicate :: Bool = false,
    )
    entry = _read_jld2(filepath, params)
    insert_runs!([entry]; db_path, allow_duplicate)
    @info "Inserted run $(entry.id) ($(length(entry.iterations)) iterations)"
    return entry
end

# Placeholder parameters: everything except model (and, for old files,
# symmetry/algorithm) is overwritten from the file itself.
_blank_params(model::AbstractString; symmetry::AbstractString = "", algorithm::AbstractString = "") =
    RunParameters(model = String(model), symmetry = String(symmetry), algorithm = String(algorithm),
                  chi = 0, K = 0, mu0_sq = 0.0, lambda = 0.0)

"""Read a current-format TNRKit JLD2 file into a `DatabaseEntry` (not inserted)."""
function _read_jld2(filepath::String, params::RunParameters)
    HDF5.h5open(filepath, "r") do f
        chi    = Int(HDF5.read(f["chi"]))
        K      = Int(HDF5.read(f["K"]))
        mu0_sq = Float64(HDF5.read(f["μ0"]))
        lambda = Float64(HDF5.read(f["λ"]))
        t      = Float64(HDF5.read(f["t"]))

        # Newer TNRKit output stores these directly; prefer them when
        # present and only fall back to `params` for older files that don't.
        symmetry  = haskey(f, "symmetry")  ? String(HDF5.read(f["symmetry"]))  : params.symmetry
        algorithm = haskey(f, "algorithm") ? String(HDF5.read(f["algorithm"])) : params.algorithm

        filled_params = RunParameters(
            model     = params.model,
            symmetry  = symmetry,
            algorithm = algorithm,
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

            norm = Float64(_field(entry, "1"))

            results = _field(entry, "2")
            cc_ref  = _field(results, "central_charge")
            cc_val  = HDF5.read(f[cc_ref])          # complex scalar (Complex or (re,im)-like)
            cc_re   = Float64(_field(cc_val, "re"))

            sd       = _field(results, "scaling_dimensions")
            data_ref = _field(sd, "data")
            sd_arr   = HDF5.read(f[data_ref])       # vector of complex-like values
            sd_re    = Float64[_field(v, "re") for v in sd_arr]

            sectors = _parse_sectors(f, _field(sd, "structure"), sd_re)

            push!(iters, CFTResults(
                iteration      = step,
                normalization  = norm,
                central_charge = cc_re,
                sectors        = sectors,
            ))
        end

        DatabaseEntry(params = filled_params, iterations = iters,
                      source_file = basename(filepath), runtime_s = t)
    end
end

# HDF5.jl reads compound HDF5 records back as `NamedTuple`s (or, for
# built-in compound layouts like complex numbers, as `Complex`) rather than
# Dict-like objects. `getfield` handles both — even fields with a purely
# numeric name (from a plain Julia `Tuple`, e.g. "1", "2") — since NamedTuple
# field names don't have to be valid Julia identifiers. Falls back to plain
# `getindex` for actual Dict/HDF5.Group objects, which have no such field.
function _field(x, key::AbstractString)
    sym = Symbol(key)
    hasfield(typeof(x), sym) && return getfield(x, sym)
    return x[key]
end

# TensorKit charge/irrep types encode a value either as a plain number (e.g.
# ZNIrrep's `n`) or, for half-integer quantum numbers (spins, U(1)/O(2)
# charges — HalfInt from HalfIntegers.jl), as the doubled integer under a
# `twice` field. This detects the latter structurally (by field name, not by
# which symmetry it belongs to) and returns the real value either way.
function _charge_component(x)
    hasfield(typeof(x), :twice) && return getfield(x, :twice) / 2
    return Float64(x)
end

# Parse the sector structure out of the JLD2 scaling_dimensions.structure
# field. Symmetry-agnostic: reads whatever fields the sector's charge/irrep
# label actually has (e.g. `charge` for U1Irrep, `n` for a ZNIrrep, `j`+`s`
# for CU1Irrep/O(2), none at all for Trivial) via reflection, rather than
# assuming a fixed (j, s)-shaped schema.
function _parse_sectors(f, struct_ref, all_dims::Vector{Float64})
    struct_val = f[struct_ref][]
    kvvec_ref  = _field(struct_val, "kvvec")
    kvvec      = f[kvvec_ref][]          # Vector of object references, one per sector

    sectors = ScalingDimSector[]
    for ref in kvvec
        pair_val = f[ref][]              # named tuple: first=sector label, second=indices ref
        # The sector label is stored inline, not as a Reference. A label with
        # no fields (TensorKit's `Trivial`) is zero-sized, and JLD2 drops
        # zero-sized fields entirely, so the pair then has only `second`.
        label      = hasfield(typeof(pair_val), :first) ? getfield(pair_val, :first) : (;)
        second_ref = _field(pair_val, "second")  # indices vector IS a Reference

        charge = Dict{String, Float64}(
            String(name) => _charge_component(getfield(label, name))
            for name in propertynames(label)
        )

        indices = Int.(f[second_ref][])   # 1-based into all_dims
        isempty(indices) && continue   # nothing found in this sector — same information as
                                        # "sector absent", so skip storing it (files have
                                        # dozens of these; they're most of the file size)
        dims = sort(all_dims[indices])

        push!(sectors, ScalingDimSector(charge=charge, dims=dims))
    end
    return sectors
end

"""
    ingest_directory!(dir; infer_params, db_path, extension=".jld2", allow_duplicate=false)

Recursively ingest every matching file under `dir` (hundreds of files is the
intended use case). `infer_params` is a function `filepath -> model`
returning just the model name for that file — the same one thing
[`ingest_jld2!`](@ref) needs — or, for finer control, a full
`RunParameters` (e.g. to force `symmetry`/`algorithm` for older files that
predate those being stored in the JLD2). Pass a closure that always returns
the same value if the whole directory shares one model, or inspect
`filepath` (e.g. a regex on its filename or parent folder) if it doesn't.

A file that is already in the database is skipped and counted, not treated
as an error. A file that fails to read or parse is logged with `@warn` and
skipped, without aborting the rest of the batch. `db/CATALOG.md` is
regenerated once at the end if anything new was inserted.

Returns `(; total, inserted, skipped, failed)`.

# Example
```julia
using TACOBELL

summary = ingest_directory!("data/loop_tnr_runs"; infer_params = _ -> "phi4_complex")
```
"""
function ingest_directory!(
        dir :: String;
        infer_params,
        db_path     :: String  = DB_PATH(),
        extension   :: String  = ".jld2",
        allow_duplicate :: Bool = false,
    )
    reader(fp, p::RunParameters) = _read_jld2(fp, p)
    reader(fp, model::AbstractString) = _read_jld2(fp, _blank_params(model))
    _ingest_directory_generic!(reader, dir; infer_params, db_path, extension, allow_duplicate)
end

# Shared resilience loop behind ingest_directory!/ingest_directory_old! — the
# only thing that differs between the two is how a single file is read, so
# it's a parameter rather than duplicated code. Files are read one by one but
# written in bulk: runs are buffered and flushed as whole shards, so ingesting
# thousands of files costs a handful of writes, not one index rewrite per file.
function _ingest_directory_generic!(
        read_fn, dir :: String;
        infer_params,
        db_path     :: String  = DB_PATH(),
        extension   :: String  = ".jld2",
        allow_duplicate :: Bool = false,
    )
    files = _find_files(dir, extension)
    isempty(files) && @warn "No *$(extension) files found under $dir"

    seen = _existing_keys(db_path)
    buffer = DatabaseEntry[]; buffered_dims = 0
    flush!() = (insert_runs!(buffer; db_path, allow_duplicate = true); empty!(buffer); buffered_dims = 0)

    n_ok = n_skip = n_err = 0
    for (i, filepath) in enumerate(files)
        label = "[$i/$(length(files))] $(basename(filepath))"
        try
            entry = read_fn(filepath, infer_params(filepath))
            k = _param_key(entry.params)
            if !allow_duplicate && haskey(seen, k)
                n_skip += 1
                @info "$label -> skipped (already in db as $(seen[k]))"
                continue
            end
            seen[k] = entry.id
            k_dims = _ndims_total(entry)
            # flush *before* overflowing, so every shard ends up full-sized
            !isempty(buffer) && buffered_dims + k_dims > SHARD_MAX_DIMS[] && flush!()
            push!(buffer, entry)
            buffered_dims += k_dims
            n_ok += 1
            @info "$label -> read"
        catch err
            # a parameter the caller has to supply would fail every file the
            # same way: stop at the first one with the explanation instead of
            # logging hundreds of identical warnings (keeping what was read)
            if err isa MissingParameterError
                isempty(buffer) || flush!()
                n_ok > 0 && generate_catalog(; db_path)
                rethrow()
            end
            n_err += 1
            @warn "$label -> FAILED" exception=(err, catch_backtrace())
        end
    end
    isempty(buffer) || flush!()

    n_ok > 0 && generate_catalog(; db_path)

    summary = (total=length(files), inserted=n_ok, skipped=n_skip, failed=n_err)
    @info "Bulk ingest complete" summary...
    return summary
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API — legacy JLD2 ingestion (pre-TNRKit-v0.5 output)
# ─────────────────────────────────────────────────────────────────────────────

"""
    ingest_jld2_old!(filepath, model; symmetry, algorithm, db_path, allow_duplicate)
    ingest_jld2_old!(filepath, params::RunParameters; db_path, allow_duplicate)

Like [`ingest_jld2!`](@ref), but for JLD2 output from **before TNRKit
v0.5.0**, which stored CFT data in a different, less structured shape: each
iteration is a plain `Dict` (keyed `"c"` for the central charge and, for
every other key, a charge/irrep label) rather than a `(normalization,
CFTData)` tuple, and the top-level bond-dimension key is named `ndimtrunc`
instead of `chi` (`K` means the same thing in both). There is no recorded
per-iteration normalization in these files, so
`CFTResults.normalization` comes back `missing` for every iteration
ingested this way.

`symmetry`/`algorithm` are never stored in old files, so unlike
[`ingest_jld2!`](@ref) there's nothing to fall back to — **both are
required keywords here, with no default**, specifically so a call that
forgets them fails loudly instead of silently inserting entries with a
blank symmetry/algorithm (which happened twice before this was enforced —
check the sector labels' charge field name to tell them apart: `charge` is
`U1Irrep` → `"U(1)"`, `n` is `ZNIrrep` → `"Z2"`/`"Z3"`/…, `j` alone is
`SU2Irrep`, `j`+`s` is `CU1Irrep` → `"O(2)"`).

# Example
```julia
using TACOBELL
ingest_jld2_old!("Com_CFT_μ0-0.01_λ0.01_K10truncrank16_niter15.jld2", "phi4_complex";
                  symmetry="U(1)", algorithm="LoopTNR")
ingest_jld2_old!("Re_CFT_TransLine_K10truncrank16_niter15_μ0-0.01_λ0.01.jld2", "phi4_real";
                  symmetry="Z2", algorithm="LoopTNR")
```
"""
function ingest_jld2_old!(
        filepath :: String,
        model    :: AbstractString;
        symmetry :: AbstractString,
        algorithm :: AbstractString,
        db_path  :: String  = DB_PATH(),
        allow_duplicate :: Bool = false,
    )
    params = _blank_params(model; symmetry, algorithm)
    ingest_jld2_old!(filepath, params; db_path, allow_duplicate)
end

function ingest_jld2_old!(
        filepath :: String,
        params   :: RunParameters;
        db_path  :: String  = DB_PATH(),
        allow_duplicate :: Bool = false,
    )
    entry = _read_jld2_old(filepath, params)
    insert_runs!([entry]; db_path, allow_duplicate)
    @info "Inserted run $(entry.id) ($(length(entry.iterations)) iterations)"
    return entry
end

"""Read a pre-TNRKit-v0.5 JLD2 file into a `DatabaseEntry` (not inserted)."""
function _read_jld2_old(filepath::String, params::RunParameters)
    HDF5.h5open(filepath, "r") do f
        chi    = Int(HDF5.read(f["ndimtrunc"]))
        K      = Int(HDF5.read(f["K"]))
        mu0_sq = Float64(HDF5.read(f["μ0"]))
        lambda = Float64(HDF5.read(f["λ"]))
        t      = Float64(HDF5.read(f["t"]))

        symmetry  = haskey(f, "symmetry")  ? String(HDF5.read(f["symmetry"]))  : params.symmetry
        algorithm = haskey(f, "algorithm") ? String(HDF5.read(f["algorithm"])) : params.algorithm

        filled_params = RunParameters(
            model=params.model, symmetry=symmetry, algorithm=algorithm,
            chi=chi, K=K, mu0_sq=mu0_sq, lambda=lambda,
        )

        data_refs = HDF5.read(f["data"])   # Vector of object references, one per iteration
        iters = CFTResults[]

        for (i, ref) in enumerate(data_refs)
            step = i - 1
            cc, sectors = _parse_old_iteration(f, ref)
            push!(iters, CFTResults(
                iteration      = step,
                normalization  = missing,
                central_charge = cc,
                sectors        = sectors,
            ))
        end

        DatabaseEntry(params = filled_params, iterations = iters,
                      source_file = basename(filepath), runtime_s = t)
    end
end

# Old-format iteration: a plain Dict serialized the same kvvec way as the
# new format's sector `structure` (see `_parse_sectors`), but here *both*
# the key and the value of each pair are References, and there's no
# separate "data"/"structure" split — a sector's value vector IS its list
# of scaling dimensions directly. One key is special: the string "c" maps
# to the central charge instead of a sector.
function _parse_old_iteration(f, data_ref)
    dict_val = f[data_ref][]
    kvvec = f[_field(dict_val, "kvvec")][]

    cc = missing
    sectors = ScalingDimSector[]
    for pair_ref in kvvec
        pair_val = f[pair_ref][]
        key = f[_field(pair_val, "first")][]
        val = f[_field(pair_val, "second")][]

        if key isa AbstractString
            if key == "c"
                cc = Float64(_field(val, "re"))
            elseif key == "sc_dim" && val isa AbstractVector
                # runs without a symmetry-resolved spectrum (e.g. the
                # "TransLine" files) store every Δ unsectored under "sc_dim":
                # keep them as one sector with an empty (trivial) charge
                dims = Float64[_field(v, "re") for v in val]
                isempty(dims) || push!(sectors, ScalingDimSector(dims=sort(dims)))
            end
            # any other string-keyed scalar is an unrecognized diagnostic — skip it
        else
            charge = Dict{String, Float64}(
                String(name) => _charge_component(getfield(key, name))
                for name in propertynames(key)
            )
            dims = Float64[_field(v, "re") for v in val]
            isempty(dims) && continue
            push!(sectors, ScalingDimSector(charge=charge, dims=sort(dims)))
        end
    end
    return cc, sectors
end

"""
    ingest_directory_old!(dir; infer_params, db_path, extension=".jld2", allow_duplicate=false)

Like [`ingest_directory!`](@ref), but for a directory of pre-TNRKit-v0.5
output — every matching file is read with [`ingest_jld2_old!`](@ref)
instead of [`ingest_jld2!`](@ref). Same resilience behavior (a bad or
already-ingested file doesn't abort the batch) and the same `infer_params`
contract (return a model name, or a full `RunParameters`).
"""
function ingest_directory_old!(
        dir :: String;
        infer_params,
        db_path     :: String  = DB_PATH(),
        extension   :: String  = ".jld2",
        allow_duplicate :: Bool = false,
    )
    reader(fp, p::RunParameters) = _read_jld2_old(fp, p)
    reader(fp, ::AbstractString) = throw(ArgumentError(
        "old-format files never store symmetry/algorithm: infer_params must return a full " *
        "RunParameters with both set, not just a model name (see ingest_jld2_old!)"))
    _ingest_directory_generic!(reader, dir; infer_params, db_path, extension, allow_duplicate)
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API — add_data!: one entry point for any TNRKit output, file or folder
# ─────────────────────────────────────────────────────────────────────────────

"""
    MissingParameterError(msg)

Thrown by [`add_data!`](@ref) when a run parameter is neither stored in the
file nor inferable from it, so it has to be passed as a keyword.
"""
struct MissingParameterError <: Exception
    msg::String
end
Base.showerror(io::IO, e::MissingParameterError) = print(io, "MissingParameterError: ", e.msg)

"""
    add_data!(path, model; symmetry, algorithm, chi, K, mu0_sq, lambda,
              db_path, extension=".jld2", allow_duplicate=false)
        -> (; total, inserted, skipped, failed)

Add TNRKit output to the database: `path` is one `.jld2` file or a folder,
searched recursively. This is the recommended way to add data; it works for
both the current and the pre-TNRKit-v0.5 file layout (detected per file, so a
folder can mix them).

`model` is the only thing you must give, since files never store it in the
database's naming (`"phi4_real"`, `"phi4_complex"`, …). For a folder mixing
models, pass a function `filepath -> model` instead of a string.

Every other parameter is read from each file. Any of them can also be passed
as a keyword, which then **overrides** the file's value for every file:

- `chi`, `K`, `mu0_sq`, `lambda`: always stored in the file.
- `symmetry`: stored in current-format files. For old files it's inferred
  from the sector labels (`charge` → `"U(1)"`, `n` (0/1 only) → `"Z2"`,
  `j` → `"SU(2)"`, `j` + `s` → `"O(2)"`). Old files without labelled
  sectors need it passed explicitly.
- `algorithm`: stored in current-format files; old files need it passed.

If a needed parameter can't be found, a [`MissingParameterError`](@ref)
stops the batch at the first such file, saying which keyword to pass (any
runs read before it are kept). Files already in the database are skipped,
and an unreadable file is logged and skipped without stopping the batch.
`db/CATALOG.md` is regenerated at the end.

# Examples
```julia
using TACOBELL
add_data!("data/O2_runs", "phi4_complex")                    # everything from the files
add_data!("data/old_runs", "phi4_complex"; algorithm="LoopTNR")
add_data!("data/transline", "phi4_real"; symmetry="Z2", algorithm="LoopTNR")
add_data!("data/mixed", fp -> occursin("Com_", basename(fp)) ? "phi4_complex" : "phi4_real")
```
From the command line: `julia --project=. scripts/add_data.jl <path> <model> [key=value ...]`.
"""
function add_data!(
        path   :: AbstractString,
        model;
        symmetry  = nothing,
        algorithm = nothing,
        chi       = nothing,
        K         = nothing,
        mu0_sq    = nothing,
        lambda    = nothing,
        db_path   :: String = DB_PATH(),
        extension :: String = ".jld2",
        allow_duplicate :: Bool = false,
    )
    overrides = (; symmetry, algorithm, chi, K, mu0_sq, lambda)
    infer = model isa AbstractString ? (_ -> String(model)) : model
    reader(fp, m) = _read_any(fp, m, overrides)
    _ingest_directory_generic!(reader, String(path); infer_params = infer, db_path, extension, allow_duplicate)
end

# `path` may be a single file (taken as is) or a folder (searched recursively).
function _find_files(path::AbstractString, extension::AbstractString)
    isfile(path) && return [String(path)]
    isdir(path) || throw(ArgumentError("no such file or directory: $path"))
    files = String[]
    for (root, _, fnames) in walkdir(path), fn in fnames
        endswith(lowercase(fn), lowercase(extension)) && push!(files, joinpath(root, fn))
    end
    return sort!(files)
end

"""Which TNRKit output layout a file uses: `:new` (v0.5+) or `:old`."""
function _jld2_format(filepath::AbstractString)
    HDF5.h5open(filepath, "r") do f
        haskey(f, "chi") && return :new
        haskey(f, "ndimtrunc") && return :old
        error("not a TNRKit output file (neither a `chi` nor an `ndimtrunc` key)")
    end
end

# Symmetry from the charge field names of the sector labels (TensorKit
# irreps: U1Irrep has `charge`, ZNIrrep `n`, SU2Irrep `j`, CU1Irrep `j`+`s`).
# `nothing` when there are no labelled sectors or the labels are ambiguous.
function _infer_symmetry(iters::Vector{CFTResults})
    fieldsets = Set(Set(keys(s.charge)) for r in iters for s in r.sectors if !isempty(s.charge))
    length(fieldsets) == 1 || return nothing
    ks = only(fieldsets)
    ks == Set(["charge"]) && return "U(1)"
    ks == Set(["j", "s"]) && return "O(2)"
    ks == Set(["j"])      && return "SU(2)"
    if ks == Set(["n"])
        # ZN with N unknown: only 0/1 seen means Z2 (Z3+ always shows n=2 too)
        maximum(s.charge["n"] for r in iters for s in r.sectors if !isempty(s.charge)) <= 1 && return "Z2"
    end
    return nothing
end

# Read either layout, then fill in / override the parameters.
function _read_any(fp::AbstractString, model, o::NamedTuple)
    model isa AbstractString ||
        throw(ArgumentError("the model function must return a String, got $(repr(model))"))
    e = _jld2_format(fp) === :new ? _read_jld2(fp, _blank_params(model)) :
                                    _read_jld2_old(fp, _blank_params(model))
    p = e.params
    pick(override, fromfile) = override === nothing ? fromfile : override
    for (k, v) in pairs(o)
        stored = getfield(p, k)
        v === nothing || stored == "" || stored == v ||
            @warn "$(basename(fp)): $k=$(repr(v)) overrides the value stored in the file ($(repr(stored)))" maxlog=5
    end

    sym = pick(o.symmetry, isempty(p.symmetry) ? something(_infer_symmetry(e.iterations), "") : p.symmetry)
    alg = pick(o.algorithm, p.algorithm)
    isempty(sym) && throw(MissingParameterError(
        "$(basename(fp)): the file doesn't record its symmetry and it can't be inferred " *
        "from the sector labels. Pass it, e.g. add_data!(path, model; symmetry=\"Z2\")."))
    isempty(alg) && throw(MissingParameterError(
        "$(basename(fp)): the file doesn't record its algorithm (pre-TNRKit-v0.5 output). " *
        "Pass it, e.g. add_data!(path, model; algorithm=\"LoopTNR\")."))

    params = RunParameters(
        model = p.model, symmetry = String(sym), algorithm = String(alg),
        chi = Int(pick(o.chi, p.chi)), K = Int(pick(o.K, p.K)),
        mu0_sq = Float64(pick(o.mu0_sq, p.mu0_sq)), lambda = Float64(pick(o.lambda, p.lambda)),
    )
    DatabaseEntry(id = e.id, created_at = e.created_at, source_file = e.source_file,
                  runtime_s = e.runtime_s, params = params, iterations = e.iterations)
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API — reading
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_entry(id; db_path) -> DatabaseEntry

Load a single entry by its UUID string.
"""
function load_entry(id::String; db_path::String = DB_PATH())
    idx = _load_index(db_path)
    i = findfirst(r -> r["id"] == id, idx)
    i === nothing && error("No entry with id=$id")
    return only(_read_entries(idx[i]["_shard"], [idx[i]["_row"]]))
end

"""
    query_runs(; model, symmetry, algorithm, chi, K, mu0_sq, lambda,
                 mu0_sq_tol, lambda_tol, db_path) -> Vector{DatabaseEntry}

Return all entries matching the given keyword filters.
Floating-point filters use `±tol` windows; all others are exact.

Omit any keyword to match all values.

# Example
```julia
# All O(2) runs with chi=16
query_runs(symmetry="O(2)", chi=16)

# Near a specific coupling point
query_runs(mu0_sq=-2.0, lambda=1.0, mu0_sq_tol=1e-8)
```
"""
function query_runs(
        ;
        db_path    :: String                  = DB_PATH(),
        model      :: Union{String,Nothing}   = nothing,
        symmetry   :: Union{String,Nothing}   = nothing,
        algorithm  :: Union{String,Nothing}   = nothing,
        chi        :: Union{Int,Nothing}      = nothing,
        K          :: Union{Int,Nothing}      = nothing,
        mu0_sq     :: Union{Float64,Nothing}  = nothing,
        lambda     :: Union{Float64,Nothing}  = nothing,
        mu0_sq_tol :: Float64                 = 1e-10,
        lambda_tol :: Float64                 = 1e-10,
    )
    _entries_for_rows(_query_rows(; db_path, model, symmetry, algorithm, chi, K,
                                    mu0_sq, lambda, mu0_sq_tol, lambda_tol))
end

function _query_rows(; db_path = DB_PATH(), model = nothing, symmetry = nothing, algorithm = nothing,
                     chi = nothing, K = nothing, mu0_sq = nothing, lambda = nothing,
                     mu0_sq_tol = 1e-10, lambda_tol = 1e-10)
    filter(_load_index(db_path)) do row
        (isnothing(model)     || row["model"]     == model)     &&
        (isnothing(symmetry)  || row["symmetry"]  == symmetry)  &&
        (isnothing(algorithm) || row["algorithm"] == algorithm) &&
        (isnothing(chi)       || row["chi"]       == chi)       &&
        (isnothing(K)         || row["K"]         == K)         &&
        (isnothing(mu0_sq)    || abs(row["mu0_sq"] - mu0_sq) < mu0_sq_tol) &&
        (isnothing(lambda)    || abs(row["lambda"] - lambda)  < lambda_tol)
    end
end

# Convenience: pass RunParameters directly
function query_runs(p::RunParameters; db_path::String = DB_PATH())
    query_runs(; db_path,
        model=p.model, symmetry=p.symmetry, algorithm=p.algorithm,
        chi=p.chi, K=p.K, mu0_sq=p.mu0_sq, lambda=p.lambda)
end

"""
    get_iteration(entry, step) -> CFTResults

Return the `CFTResults` for a specific iteration step from an entry.
"""
function get_iteration(entry::DatabaseEntry, step::Int)
    for r in entry.iterations
        r.iteration == step && return r
    end
    error("No iteration step $step in entry $(entry.id)")
end

"""
    final_iteration(entry) -> CFTResults

The stored iteration with the highest RG step index — the usual "give me
the answer" shortcut, equivalent to
`get_iteration(entry, maximum(r.iteration for r in entry.iterations))`.
"""
final_iteration(entry::DatabaseEntry) = argmax(r -> r.iteration, entry.iterations)

"""
    find_sector(result; charge...) -> Union{ScalingDimSector, Nothing}

Look up one symmetry sector's scaling dimensions from a single iteration's
result by its charge/irrep label, given as keywords matching whatever
quantum number(s) that symmetry uses — the same names shown by `sec.charge`
or by printing a `ScalingDimSector`. All components must be given. Returns
`nothing` if no sector matches (including "wasn't populated at this
iteration" — see [`ScalingDimSector`](@ref)).

# Example
```julia
find_sector(result; charge=1)         # U(1): "charge" is U1Irrep's field name
find_sector(result; n=1)              # Z_N: "n" is ZNIrrep's field name
find_sector(result; j=1, s=2)         # O(2) / CU1Irrep
find_sector(result)                   # Trivial (no symmetry): no keywords
```
"""
function find_sector(result::CFTResults; charge...)
    query = Dict{String, Float64}(String(k) => Float64(v) for (k, v) in charge)
    for sec in result.sectors
        sec.charge == query && return sec
    end
    return nothing
end

"""
    central_charge_trajectory(entry) -> Vector{Float64}

The central charge at every stored iteration that has one, in order of
increasing RG step. Feed this to [`plateau_estimate`](@ref), or plot it
directly to see the RG flow.
"""
function central_charge_trajectory(entry::DatabaseEntry)
    sorted = sort(entry.iterations, by = r -> r.iteration)
    return Float64[r.central_charge for r in sorted if r.central_charge !== missing]
end

"""
    plateau_estimate(values; nwin=5, skipfrac=0.3, avoid_zero=false, zero_tol=1e-6,
                      recency_tol=1.5)
        -> (value, err, range, suspect)

Estimate a converged value from a noisy RG trajectory by scanning windows of
`nwin` consecutive values and keeping the flattest one (smallest standard
deviation), rather than trusting the literal last value — TNR trajectories
are often stable for a stretch and then drift or blow up in the last few
steps as finite-χ truncation error compounds. The first `skipfrac` fraction
of `values` is excluded from the search, since it's usually a lattice-scale
transient rather than the converged plateau.

Near a critical point, a trajectory can also "walk" close to a non-trivial
fixed point for several steps before eventually flowing to its true IR
fixed point (often the trivial `c = 0` one) — and that earlier walk can
occasionally be *marginally* flatter than the genuine, later plateau it
flows to. To avoid locking onto that earlier false plateau, ties are broken
by recency: among windows whose standard deviation is within a factor
`recency_tol` of the global minimum, the *latest* one is used, not just
whichever happens to be flattest. A trajectory that genuinely blows up late
(the case this function was originally written to guard against) has a
window standard deviation well above that tolerance, so it's unaffected.

Set `avoid_zero=true` to also exclude windows containing a value with
`abs(x) < zero_tol` when a clean window exists elsewhere — useful for
scaling dimensions, where an under-converged eigensolver can pad a sector
with a numerical near-zero. Leave it `false` for the central charge, where
`c → 0` can be a genuine physical answer.

Returns `value` (the window mean), `err` (the window's standard deviation,
a rough convergence error bar), `range` (which indices into `values` were
used), and `suspect` (true if no window satisfying `avoid_zero` was found,
so `range` had to fall back to the plain flattest window).
"""
function plateau_estimate(
        values::AbstractVector{<:Real};
        nwin::Int = 5, skipfrac::Float64 = 0.3,
        avoid_zero::Bool = false, zero_tol::Float64 = 1.0e-6,
        recency_tol::Float64 = 1.5,
    )
    n = length(values)
    n == 0 && error("plateau_estimate: empty trajectory")
    w = min(nwin, n)
    windows = [values[i:(i + w - 1)] for i in 1:(n - w + 1)]
    # Population (uncorrected) std: every window has the same width, so this
    # is just a uniform rescaling of the corrected estimator and doesn't
    # change any window's ranking -- but it stays finite for single-point
    # windows (w == 1), where the corrected estimator is NaN (division by
    # n - 1 == 0), which would otherwise poison every comparison below.
    stds = Statistics.std.(windows; corrected = false)
    start_min = max(1, ceil(Int, skipfrac * n))
    isbad(win) = avoid_zero && any(x -> abs(x) < zero_tol, win)
    clean = findall(i -> !isbad(windows[i]) && i >= start_min, eachindex(windows))
    candidates = isempty(clean) ? eachindex(windows) : clean
    best_std = minimum(stds[candidates])
    near_best = filter(i -> stds[i] <= recency_tol * best_std, candidates)
    i_best = maximum(near_best)
    return (
        value   = Statistics.mean(windows[i_best]),
        err     = stds[i_best],
        range   = i_best:(i_best + w - 1),
        suspect = isempty(clean),
    )
end

"""
    plateau_central_charge(entry; kwargs...)

Shorthand for `plateau_estimate(central_charge_trajectory(entry); kwargs...)`.
"""
plateau_central_charge(entry::DatabaseEntry; kwargs...) =
    plateau_estimate(central_charge_trajectory(entry); kwargs...)

"""
    plot_central_charge(entry; kwargs...)

Plot the central charge as a function of RG iteration for one run — the
quickest way to eyeball whether/where it converged. `kwargs` are passed
through to `Plots.plot` (e.g. `title`, `ylims`).

# Example
```julia
using TACOBELL
plot_central_charge(runs[1])
```
"""
function plot_central_charge(entry::DatabaseEntry; kwargs...)
    sorted = sort(entry.iterations, by = r -> r.iteration)
    xs = Int[r.iteration for r in sorted if r.central_charge !== missing]
    ys = Float64[r.central_charge for r in sorted if r.central_charge !== missing]

    p = entry.params
    Plots.plot(xs, ys;
        marker = :circle, legend = false,
        xlabel = "RG iteration", ylabel = "central charge c",
        title = "$(p.model), $(p.symmetry), χ=$(p.chi), μ₀²=$(p.mu0_sq), λ=$(p.lambda)",
        titlefontsize = 10,
        kwargs...,
    )
end

"""
    find_closest(mu0_sq, lambda; filters...) -> Union{DatabaseEntry, Nothing}

Among entries matching the discrete keyword filters, return the one whose
(μ₀², λ) is nearest in Euclidean distance to `(mu0_sq, lambda)`.
"""
function find_closest(
        mu0_sq::Float64, lambda::Float64;
        model=nothing, symmetry=nothing, algorithm=nothing,
        chi=nothing, K=nothing, db_path=DB_PATH(),
    )
    rows = _query_rows(; db_path, model, symmetry, algorithm, chi, K)
    isempty(rows) && return nothing
    dist(r) = (r["mu0_sq"] - mu0_sq)^2 + (r["lambda"] - lambda)^2
    return only(_entries_for_rows([rows[argmin(dist.(rows))]]))
end

"""
    rebuild_index!(; db_path) -> Int

The run index is derived from the shard files on the fly and cached in
memory (nothing is stored on disk, so it can never go stale or conflict in
a merge). This drops that cache and re-reads every shard — only needed if
shard files were changed in a way their modification time doesn't reflect.
Returns the number of runs.
"""
function rebuild_index!(; db_path::String = DB_PATH())
    lock(_CACHE_LOCK) do
        empty!(_META_CACHE)
        empty!(_DETAIL_CACHE)
    end
    n = length(_load_index(db_path))
    @info "Rebuilt index: $n entries"
    return n
end

"""
    compact!(; db_path, max_dims=TACOBELL.SHARD_MAX_DIMS[]) -> Int

Merge small shard files (e.g. left behind by many single [`insert_run!`](@ref)
calls) into as few full-size shards as possible. Shards that are already at
least half full are left alone. Returns the number of shard files removed.

Unlike every other write, this *replaces* existing files, so git stores the
merged data again — run it occasionally, not after every insert.
"""
function compact!(; db_path::String = DB_PATH(), max_dims::Int = SHARD_MAX_DIMS[])
    removed = 0
    groups = Dict{DataType,Vector{String}}()
    for path in _shard_files(db_path)
        n = h5open(f -> length(f["dims"]), path, "r")
        n < max_dims ÷ 2 && push!(get!(groups, _shard_dims_type(path), String[]), path)
    end
    for (T, paths) in groups
        length(paths) < 2 && continue
        entries = reduce(vcat, (_read_shard(p) for p in paths))
        for chunk in _chunk_by_dims(entries; max_dims)
            _write_shard(_new_shard_path(db_path), chunk; dims_type = T)
        end
        foreach(rm, paths)
        removed += length(paths)
    end
    @info "compact!: merged $removed small shard(s)"
    return removed
end

"""
    import_toml_runs!(dir_or_files; db_path, allow_duplicate=false) -> Int

Import runs stored in the old one-TOML-file-per-run format (`<uuid>.toml`
files, as in `db/runs/` before the switch to HDF5 shards), keeping their ids
and creation dates. Pass a directory, or a vector of file paths to control
the order runs are grouped into shards. Runs whose parameters are already
in the database are skipped. Files are streamed in bulk, so memory use
stays bounded. Returns the number of runs imported.
"""
import_toml_runs!(dir::String; kwargs...) =
    import_toml_runs!(sort!(filter(f -> endswith(f, ".toml"), readdir(dir; join = true))); kwargs...)

function import_toml_runs!(files::Vector{String}; db_path::String = DB_PATH(), allow_duplicate::Bool = false)
    seen = _existing_keys(db_path)
    buffer = DatabaseEntry[]; nd = 0; n = 0
    for (i, f) in enumerate(files)
        e = _entry_from_dict(TOML.parsefile(f))
        k = _param_key(e.params)
        if !allow_duplicate && haskey(seen, k)
            @warn "Skipping $(basename(f)): same parameters as $(seen[k])"
            continue
        end
        seen[k] = e.id
        k_dims = _ndims_total(e)
        # flush *before* overflowing, so every shard ends up full-sized
        if !isempty(buffer) && nd + k_dims > SHARD_MAX_DIMS[]
            insert_runs!(buffer; db_path, allow_duplicate = true)
            empty!(buffer); nd = 0
        end
        push!(buffer, e); nd += k_dims; n += 1
        i % 1000 == 0 && @info "import_toml_runs!: $i / $(length(files)) files read"
    end
    insert_runs!(buffer; db_path, allow_duplicate = true)
    @info "Imported $n run(s) from $(length(files)) file(s)"
    return n
end

"""
    correct_field!(field, old, new; db_path) -> Int

Fix a typo (or rename a value) in one `RunParameters` field across every
matching entry already in the database, in place — same `id`, only that
field changes. `field` is `:model`, `:symmetry`, or `:algorithm`. The shard
files holding affected runs are rewritten (keeping their scaling-dimension
precision); `db/CATALOG.md` is regenerated afterward if anything changed.
Returns the number of entries fixed.

# Example
```julia
correct_field!(:model, "phi_complex", "phi4_complex")
```
"""
function correct_field!(field::Symbol, old, new; db_path::String = DB_PATH())
    field in (:model, :symmetry, :algorithm) ||
        error("correct_field! only supports :model, :symmetry, or :algorithm, got :$field")

    n = 0
    affected = unique(row["_shard"] for row in _load_index(db_path) if row[String(field)] == old)
    for path in affected
        T = _shard_dims_type(path)
        entries = map(_read_shard(path)) do entry
            getfield(entry.params, field) == old || return entry
            n += 1
            fixed_params = RunParameters(;
                (k => (k === field ? new : getfield(entry.params, k))
                 for k in fieldnames(RunParameters))...
            )
            DatabaseEntry(;
                id=entry.id, created_at=entry.created_at, source_file=entry.source_file,
                runtime_s=entry.runtime_s, params=fixed_params, iterations=entry.iterations,
            )
        end
        _write_shard(path, entries; dims_type = T)
    end

    n > 0 && generate_catalog(; db_path)
    @info "correct_field!" field old new fixed=n
    return n
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API — inspection
# ─────────────────────────────────────────────────────────────────────────────

"""
    summarize_db(; db_path)

Print a human-readable table of all runs.
"""
function summarize_db(; db_path::String = DB_PATH())
    idx = _load_index(db_path)
    println("CFT Database — $(length(idx)) entries")
    println("─"^88)
    @printf("%-8s  %-12s  %-6s  %-10s  %4s  %4s  %8s  %8s  %5s  %6s\n",
            "ID", "model", "sym", "algorithm", "χ", "K",
            "μ₀²", "λ", "iters", "c")
    println("─"^88)
    for row in sort(idx, by=r->r["created_at"])
        c = isnan(row["central_charge"]) ? "  —   " :
            @sprintf("%6.4f", row["central_charge"])
        @printf("%-8s  %-12s  %-6s  %-10s  %4d  %4d  %8.4f  %8.4f  %5d  %s\n",
                row["id"][1:8],
                row["model"], row["symmetry"], row["algorithm"],
                row["chi"], row["K"],
                row["mu0_sq"], row["lambda"],
                row["n_iterations"], c)
    end
end

"""
    list_algorithms(; db_path) -> Vector{String}
    list_symmetries(; db_path) -> Vector{String}

Inspect which distinct values are present in the database.
"""
list_algorithms(; db_path=DB_PATH()) = unique(r["algorithm"] for r in _load_index(db_path))
list_symmetries(; db_path=DB_PATH()) = unique(r["symmetry"]  for r in _load_index(db_path))

# ─────────────────────────────────────────────────────────────────────────────
# Public API — export
# ─────────────────────────────────────────────────────────────────────────────

_csv_field(x::AbstractString) =
    (occursin(",", x) || occursin("\"", x) || occursin("\n", x)) ?
        "\"" * replace(x, "\"" => "\"\"") * "\"" : x
_csv_field(::Missing) = ""
_csv_field(x) = string(x)
_csv_row(vals) = join(_csv_field.(vals), ",")

"""
    export_csv(entries=query_runs(); out="tacobell_export.csv") -> String

Write a summary table — one row per run — to a CSV file: the simplest way
to get data out of the database and into Excel, pandas, R, or anything else
that isn't Julia. Pass a filtered result from [`query_runs`](@ref) or
[`find_closest`](@ref) to export just a selection instead of everything.
Returns the path written.

# Example
```julia
export_csv(query_runs(symmetry="O(2)"); out="o2_results.csv")
```
"""
function export_csv(entries::Vector{DatabaseEntry} = query_runs(); out::String = "tacobell_export.csv")
    io = IOBuffer()
    println(io, _csv_row((
        "id", "model", "symmetry", "algorithm", "chi", "K", "mu0_sq", "lambda",
        "n_iterations", "last_iteration", "central_charge", "source_file",
    )))
    for e in entries
        p = e.params
        cc = _best_c(e.iterations)
        last_it = isempty(e.iterations) ? missing : maximum(r -> r.iteration, e.iterations)
        println(io, _csv_row((
            e.id, p.model, p.symmetry, p.algorithm, p.chi, p.K, p.mu0_sq, p.lambda,
            length(e.iterations), last_it, isnan(cc) ? missing : cc, e.source_file,
        )))
    end
    write(out, String(take!(io)))
    @info "Wrote $(length(entries)) row(s) to $out"
    return out
end

"""
    export_csv(entry::DatabaseEntry; out=entry.id[1:8]*".csv") -> String

Write one run's full per-iteration, per-sector scaling-dimension data to a
CSV file in long format (one row per (iteration, sector, dimension)) — a
way to get everything in a single ~hundreds-of-KB run file into a flat
table for plotting or analysis outside Julia. Returns the path written.
"""
function export_csv(entry::DatabaseEntry; out::String = first(entry.id, 8) * ".csv")
    # Column names for the charge/irrep quantum number(s) depend on the
    # symmetry (e.g. "j","s" for O(2); "charge" for U(1); "n" for Z_N) —
    # collect whichever ones actually appear so this works for any symmetry.
    charge_keys = sort(unique(k for r in entry.iterations for sec in r.sectors for k in keys(sec.charge)))

    io = IOBuffer()
    println(io, _csv_row((
        "iteration", "normalization", "central_charge", charge_keys..., "dim_index", "delta",
    )))
    for r in sort(entry.iterations, by = x -> x.iteration)
        cc = r.central_charge === missing ? missing : r.central_charge
        if isempty(r.sectors)
            println(io, _csv_row((r.iteration, r.normalization, cc, fill(missing, length(charge_keys))..., missing, missing)))
        end
        for sec in r.sectors
            charge_vals = (get(sec.charge, k, missing) for k in charge_keys)
            for (i, d) in enumerate(sec.dims)
                println(io, _csv_row((r.iteration, r.normalization, cc, charge_vals..., i, d)))
            end
        end
    end
    write(out, String(take!(io)))
    @info "Wrote full trajectory for $(first(entry.id, 8))… to $out"
    return out
end

# Minimal hand-rolled JSON writer — the data here is just numbers, strings,
# dicts and arrays (plus `missing`/NaN, mapped to `null`), so this avoids
# pulling in a JSON dependency for one export function.
_to_json(x::Missing) = "null"
_to_json(x::Bool) = x ? "true" : "false"
_to_json(x::Integer) = string(x)
# Rounded to 7 significant digits: still far more precision than physically
# meaningful at any χ this package deals with, but roughly halves the size
# of exports with thousands of runs' worth of scaling dimensions in them.
_to_json(x::AbstractFloat) = (isnan(x) || isinf(x)) ? "null" : string(round(x, sigdigits=7))
_to_json(x::AbstractString) = "\"" * replace(replace(replace(replace(replace(
    x, "\\" => "\\\\"), "\"" => "\\\""), "\n" => "\\n"), "\r" => "\\r"), "\t" => "\\t") * "\""
_to_json(x::AbstractDict) = "{" * join((_to_json(string(k)) * ":" * _to_json(v) for (k, v) in x), ",") * "}"
_to_json(x::Union{AbstractVector, Tuple}) = "[" * join((_to_json(v) for v in x), ",") * "]"

"""
    export_json(entries=query_runs(); out="tacobell_export.json") -> String

Write a lightweight **index** of the database to a single JSON file: one
object per run with its parameters and (plateau-based) central charge, but
not its per-iteration trajectory or sectors — see the `DatabaseEntry`
method below for that. A portable, language-agnostic snapshot for any other
tool to consume (the Explorer webpage uses the more compact
[`export_web`](@ref) instead). Pass a filtered result from
[`query_runs`](@ref) to export just a selection. Returns the path written.
"""
function export_json(entries::Vector{DatabaseEntry} = query_runs(); out::String = "tacobell_export.json")
    runs = map(entries) do e
        p = e.params
        Dict{String, Any}(
            "id" => e.id, "model" => p.model, "symmetry" => p.symmetry, "algorithm" => p.algorithm,
            "chi" => p.chi, "K" => p.K, "mu0_sq" => p.mu0_sq, "lambda" => p.lambda,
            "n_iterations" => length(e.iterations),
            "central_charge" => (c = _best_c(e.iterations); isnan(c) ? missing : c),
        )
    end
    payload = Dict{String, Any}("generated_at" => string(now()), "runs" => runs)
    write(out, _to_json(payload))
    @info "Wrote $(length(entries)) run(s) to $out"
    return out
end

"""
    export_json(entry::DatabaseEntry; out=entry.id*".json") -> String

Write one run's full per-iteration **detail** — normalization, central
charge, and every populated sector's scaling dimensions, at every stored
iteration — to a single JSON file, at full precision. The counterpart to the
index written by the `Vector{DatabaseEntry}` method. Returns the path
written.
"""
function export_json(entry::DatabaseEntry; out::String = entry.id * ".json")
    iters = [
        Dict{String, Any}(
            "iteration" => r.iteration,
            "normalization" => r.normalization,
            "central_charge" => r.central_charge,
            "sectors" => [Dict{String, Any}("charge" => sec.charge, "dims" => sec.dims) for sec in r.sectors],
        )
        for r in sort(entry.iterations, by = x -> x.iteration)
    ]
    payload = Dict{String, Any}("id" => entry.id, "iterations" => iters)
    write(out, _to_json(payload))
    @info "Wrote full trajectory for $(first(entry.id, 8))… to $out"
    return out
end

"""
    export_json_runs(entries=query_runs(); dir="tacobell_runs_json") -> Int

Write the per-run detail JSON (see [`export_json`](@ref)'s `DatabaseEntry`
method) for every given entry into `dir`, one `<id>.json` per run. Pass a
filtered selection: the whole database as JSON is several GB. Returns the
number of files written.
"""
function export_json_runs(entries::Vector{DatabaseEntry} = query_runs(); dir::String = "tacobell_runs_json")
    mkpath(dir)
    for e in entries
        export_json(e; out = joinpath(dir, e.id * ".json"))
    end
    @info "Wrote $(length(entries)) run detail file(s) to $dir"
    return length(entries)
end

"""
    generate_catalog(; db_path, out) -> String

Write a Markdown overview of the database to `out` (default
`<db_path>/CATALOG.md`): one row per setup (model, symmetry, algorithm, χ,
K) with its number of runs, the λ values and μ₀² range covered. This is
what renders directly on GitHub for anyone browsing the repository without
Julia installed; for individual runs, use the Explorer webpage or
[`query_runs`](@ref). Regenerated automatically by
[`ingest_directory!`](@ref). Returns the path written.
"""
function generate_catalog(; db_path::String = DB_PATH(), out::String = joinpath(db_path, "CATALOG.md"))
    idx = _load_index(db_path)
    groups = Dict{Tuple,Vector{Dict{String,Any}}}()
    for r in idx
        push!(get!(groups, (r["model"], r["symmetry"], r["algorithm"], r["chi"], r["K"]), Dict{String,Any}[]), r)
    end

    fmt(x) = @sprintf("%g", x)
    function lambdas(rs)
        ls = sort(unique(r["lambda"] for r in rs))
        length(ls) <= 6 ? join(fmt.(ls), ", ") : "$(length(ls)) values in [$(fmt(first(ls))), $(fmt(last(ls)))]"
    end

    io = IOBuffer()
    println(io, "<!-- Auto-generated by `TACOBELL.generate_catalog()`. Do not edit by hand. -->")
    println(io)
    println(io, "# TACO-BELL catalog")
    println(io)
    println(io, "$(length(idx)) run(s) in $(length(groups)) setup(s), stored in ",
                "$(length(_shard_files(db_path))) shard file(s) under `shards/`.")
    println(io)
    println(io, "| model | symmetry | algorithm | χ | K | runs | λ | μ₀² range |")
    println(io, "|---|---|---|---|---|---|---|---|")
    for key in sort!(collect(keys(groups)))
        rs = groups[key]
        mus = [r["mu0_sq"] for r in rs]
        println(io, "| ", join(string.(key), " | "), " | ", length(rs), " | ", lambdas(rs),
                " | ", fmt(minimum(mus)), " … ", fmt(maximum(mus)), " |")
    end

    mkpath(dirname(out))
    write(out, String(take!(io)))
    @info "Wrote catalog ($(length(groups)) setups, $(length(idx)) runs) to $out"
    return out
end

include("webexport.jl")

end # module

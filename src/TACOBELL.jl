"""
    TACOBELL

Tensor Archive of Conformal Output — Best Ever Lattice Labour.

A lightweight, file-based database for CFT data (central charges, scaling
dimensions, ...) extracted from Tensor Network Renormalization (TNR)
calculations, designed for use alongside
[TNRKit](https://github.com/QuantumKitHub/TNRKit.jl/).

See the package README for a full walkthrough. The core entry points are
[`insert_run!`](@ref)/[`ingest_jld2!`](@ref)/[`ingest_directory!`](@ref) (add
data), [`query_runs`](@ref)/[`find_closest`](@ref) (look it up), and
[`export_csv`](@ref)/[`generate_catalog`](@ref) (get it out again, as a CSV
file or as a Markdown table for browsing on GitHub).
"""
module TACOBELL

using TOML, UUIDs, Dates, Printf, Statistics
using HDF5

export RunParameters, ScalingDimSector, CFTResults, DatabaseEntry
export insert_run!, query_runs, load_entry, get_iteration, final_iteration, summarize_db
export find_closest, find_sector, list_algorithms, list_symmetries, rebuild_index!, correct_field!
export ingest_jld2!, ingest_directory!, generate_catalog, export_csv, export_json
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
- `normalization`   : log of the per-site tensor norm at this step
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
    normalization   :: Float64
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
    populated = count(s -> !isempty(s.dims), r.sectors)
    print(io, "CFTResults(iteration=", r.iteration, ", norm=", @sprintf("%.4f", r.normalization),
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
# Serialisation helpers
# ─────────────────────────────────────────────────────────────────────────────

DB_PATH() = joinpath(@__DIR__, "..", "db")

function _sector_to_dict(sec::ScalingDimSector)
    Dict("charge" => sec.charge, "dims" => sec.dims)
end

function _sector_from_dict(d)
    if haskey(d, "charge")
        charge = Dict{String,Float64}(String(k) => Float64(v) for (k, v) in d["charge"])
    else
        # Pre-generalization files only ever stored O(2) sectors as (twice_j, s).
        charge = Dict{String,Float64}("j" => d["twice_j"] / 2, "s" => Float64(d["s"]))
    end
    ScalingDimSector(charge=charge, dims=Float64.(d["dims"]))
end

function _results_to_dict(r::CFTResults)
    d = Dict{String,Any}(
        "iteration"     => r.iteration,
        "normalization" => r.normalization,
        "notes"         => r.notes,
        "sectors"       => [_sector_to_dict(s) for s in r.sectors],
    )
    r.central_charge === missing || (d["central_charge"] = r.central_charge)
    d
end

function _results_from_dict(d)
    CFTResults(
        iteration       = d["iteration"],
        normalization   = d["normalization"],
        central_charge  = get(d, "central_charge",  missing),
        sectors         = [_sector_from_dict(s) for s in get(d, "sectors", [])],
        notes           = get(d, "notes", ""),
    )
end

function _params_to_dict(p::RunParameters)
    Dict(
        "model"     => p.model,
        "symmetry"  => p.symmetry,
        "algorithm" => p.algorithm,
        "chi"       => p.chi,
        "K"         => p.K,
        "mu0_sq"    => p.mu0_sq,
        "lambda"    => p.lambda,
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

function _entry_to_dict(e::DatabaseEntry)
    d = Dict{String,Any}(
        "id"          => e.id,
        "created_at"  => e.created_at,
        "source_file" => e.source_file,
        "params"      => _params_to_dict(e.params),
        "iterations"  => [_results_to_dict(r) for r in e.iterations],
    )
    e.runtime_s === missing || (d["runtime_s"] = e.runtime_s)
    d
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

# ─────────────────────────────────────────────────────────────────────────────
# Index helpers
# ─────────────────────────────────────────────────────────────────────────────

_index_path(db_path=DB_PATH()) = joinpath(db_path, "index.toml")

function _load_index(db_path=DB_PATH())
    path = _index_path(db_path)
    isfile(path) || return Dict{String,Any}[]
    get(TOML.parsefile(path), "entry", Dict{String,Any}[])
end

function _save_index(entries, db_path=DB_PATH())
    mkpath(db_path)
    open(_index_path(db_path), "w") do io
        TOML.print(io, Dict("entry" => entries))
    end
end

function _best_c(iters::Vector{CFTResults})
    # The literal last iteration can be numerically unstable (finite-χ
    # truncation error compounding under the RG flow), so use the flattest
    # plateau rather than trusting it blindly — see plateau_estimate.
    traj = Float64[r.central_charge for r in sort(iters, by = r -> r.iteration) if r.central_charge !== missing]
    isempty(traj) && return NaN
    return plateau_estimate(traj).value
end

function _index_row(e::DatabaseEntry)
    p = e.params
    Dict(
        "id"             => e.id,
        "created_at"     => e.created_at,
        "model"          => p.model,
        "symmetry"       => p.symmetry,
        "algorithm"      => p.algorithm,
        "chi"            => p.chi,
        "K"              => p.K,
        "mu0_sq"         => p.mu0_sq,
        "lambda"         => p.lambda,
        "n_iterations"   => length(e.iterations),
        "central_charge" => _best_c(e.iterations),
    )
end

# ─────────────────────────────────────────────────────────────────────────────
# Public API — writing
# ─────────────────────────────────────────────────────────────────────────────

"""
    insert_run!(params, iterations; db_path, allow_duplicate, source_file, runtime_s)

Insert a new entry (parameter set + all iteration results) into the database.
Returns the created `DatabaseEntry`.

Raises an error if an identical parameter set is already in the database,
unless `allow_duplicate=true`.
"""
function insert_run!(
        params      :: RunParameters,
        iterations  :: Vector{CFTResults};
        db_path     :: String             = DB_PATH(),
        allow_duplicate :: Bool           = false,
        source_file :: String             = "",
        runtime_s   :: Union{Float64,Missing} = missing,
    )
    if !allow_duplicate
        existing = query_runs(params; db_path)
        isempty(existing) || error(
            "Duplicate: entry $(existing[1].id) has identical parameters. " *
            "Pass allow_duplicate=true to insert anyway."
        )
    end

    entry = DatabaseEntry(
        params       = params,
        iterations   = iterations,
        source_file  = source_file,
        runtime_s    = runtime_s,
    )

    run_dir = joinpath(db_path, "runs")
    mkpath(run_dir)
    open(joinpath(run_dir, entry.id * ".toml"), "w") do io
        TOML.print(io, _entry_to_dict(entry))
    end

    idx = _load_index(db_path)
    push!(idx, _index_row(entry))
    _save_index(idx, db_path)

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
    params = RunParameters(
        model=String(model), symmetry=String(symmetry), algorithm=String(algorithm),
        chi=0, K=0, mu0_sq=0.0, lambda=0.0,
    )
    ingest_jld2!(filepath, params; db_path, allow_duplicate)
end

function ingest_jld2!(
        filepath :: String,
        params   :: RunParameters;
        db_path  :: String  = DB_PATH(),
        allow_duplicate :: Bool = false,
    )
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

        insert_run!(filled_params, iters;
                    db_path, allow_duplicate,
                    source_file=basename(filepath),
                    runtime_s=t)
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
        label      = _field(pair_val, "first")   # sector label is stored inline, not a Reference
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
            ingest_jld2!(filepath, params; db_path, allow_duplicate)
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

# ─────────────────────────────────────────────────────────────────────────────
# Public API — reading
# ─────────────────────────────────────────────────────────────────────────────

"""
    load_entry(id; db_path) -> DatabaseEntry

Load a single entry by its UUID string.
"""
function load_entry(id::String; db_path::String = DB_PATH())
    path = joinpath(db_path, "runs", id * ".toml")
    isfile(path) || error("No entry with id=$id")
    _entry_from_dict(TOML.parsefile(path))
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
    idx = _load_index(db_path)
    hits = filter(idx) do row
        (isnothing(model)     || row["model"]     == model)     &&
        (isnothing(symmetry)  || row["symmetry"]  == symmetry)  &&
        (isnothing(algorithm) || row["algorithm"] == algorithm) &&
        (isnothing(chi)       || row["chi"]       == chi)       &&
        (isnothing(K)         || row["K"]         == K)         &&
        (isnothing(mu0_sq)    || abs(row["mu0_sq"] - mu0_sq) < mu0_sq_tol) &&
        (isnothing(lambda)    || abs(row["lambda"] - lambda)  < lambda_tol)
    end
    [load_entry(row["id"]; db_path) for row in hits]
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
    plateau_estimate(values; nwin=5, skipfrac=0.3, avoid_zero=false, zero_tol=1e-6)
        -> (value, err, range, suspect)

Estimate a converged value from a noisy RG trajectory by scanning windows of
`nwin` consecutive values and keeping the flattest one (smallest standard
deviation), rather than trusting the literal last value — TNR trajectories
are often stable for a stretch and then drift or blow up in the last few
steps as finite-χ truncation error compounds. The first `skipfrac` fraction
of `values` is excluded from the search, since it's usually a lattice-scale
transient rather than the converged plateau.

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
    )
    n = length(values)
    n == 0 && error("plateau_estimate: empty trajectory")
    w = min(nwin, n)
    windows = [values[i:(i + w - 1)] for i in 1:(n - w + 1)]
    stds = Statistics.std.(windows)
    start_min = max(1, ceil(Int, skipfrac * n))
    isbad(win) = avoid_zero && any(x -> abs(x) < zero_tol, win)
    clean = findall(i -> !isbad(windows[i]) && i >= start_min, eachindex(windows))
    candidates = isempty(clean) ? eachindex(windows) : clean
    i_best = candidates[argmin(stds[candidates])]
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

!!! note "Requires Plots.jl"
    This method is provided by a package extension and only becomes
    available once you `using Plots`. Run `using Pkg; Pkg.add("Plots")`
    once if you don't already have it — note that, like any `Pkg.add`,
    this will add Plots as an ordinary dependency of this environment.

# Example
```julia
using TACOBELL, Plots
plot_central_charge(runs[1])
```
"""
function plot_central_charge end

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
    cands = query_runs(; db_path, model, symmetry, algorithm, chi, K)
    isempty(cands) && return nothing
    dist(e) = (e.params.mu0_sq - mu0_sq)^2 + (e.params.lambda - lambda)^2
    cands[argmin(dist.(cands))]
end

"""
    rebuild_index!(; db_path)

Rebuild `index.toml` by scanning all run files.
Run this after merging a pull request that adds run files without updating the index.
"""
function rebuild_index!(; db_path::String = DB_PATH())
    run_dir = joinpath(db_path, "runs")
    isdir(run_dir) || return
    idx = Dict{String,Any}[]
    for f in readdir(run_dir, join=true)
        endswith(f, ".toml") || continue
        try
            e = _entry_from_dict(TOML.parsefile(f))
            push!(idx, _index_row(e))
        catch err
            @warn "Skipping $f: $err"
        end
    end
    _save_index(idx, db_path)
    @info "Rebuilt index: $(length(idx)) entries"
end

"""
    correct_field!(field, old, new; db_path) -> Int

Fix a typo (or rename a value) in one `RunParameters` field across every
matching entry already in the database, in place — same `id`, same file,
only that field changes. `field` is `:model`, `:symmetry`, or `:algorithm`.
Rebuilds `index.toml` and regenerates `db/CATALOG.md` afterward if anything
changed. Returns the number of entries fixed.

# Example
```julia
correct_field!(:model, "phi_complex", "phi4_complex")
```
"""
function correct_field!(field::Symbol, old, new; db_path::String = DB_PATH())
    field in (:model, :symmetry, :algorithm) ||
        error("correct_field! only supports :model, :symmetry, or :algorithm, got :$field")

    n = 0
    for row in _load_index(db_path)
        row[String(field)] == old || continue
        entry = load_entry(row["id"]; db_path)
        fixed_params = RunParameters(;
            (k => (k === field ? new : getfield(entry.params, k))
             for k in fieldnames(RunParameters))...
        )
        fixed_entry = DatabaseEntry(;
            id=entry.id, created_at=entry.created_at, source_file=entry.source_file,
            runtime_s=entry.runtime_s, params=fixed_params, iterations=entry.iterations,
        )
        open(joinpath(db_path, "runs", entry.id * ".toml"), "w") do io
            TOML.print(io, _entry_to_dict(fixed_entry))
        end
        n += 1
    end

    if n > 0
        rebuild_index!(; db_path)
        generate_catalog(; db_path)
    end
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
_to_json(x::Real) = (isnan(x) || isinf(x)) ? "null" : string(x)
_to_json(x::AbstractString) = "\"" * replace(replace(replace(replace(replace(
    x, "\\" => "\\\\"), "\"" => "\\\""), "\n" => "\\n"), "\r" => "\\r"), "\t" => "\\t") * "\""
_to_json(x::AbstractDict) = "{" * join((_to_json(string(k)) * ":" * _to_json(v) for (k, v) in x), ",") * "}"
_to_json(x::Union{AbstractVector, Tuple}) = "[" * join((_to_json(v) for v in x), ",") * "]"

"""
    export_json(entries=query_runs(); out="tacobell_export.json") -> String

Write a summary of the database to a single JSON file: one object per run,
with its parameters and the scaling-dimension sectors of its final
iteration (not the full per-iteration trajectory — use [`export_csv`](@ref)
for that). This is what feeds the static "browse the database" webpage
(see `web/` in the repository), but is also just a portable,
language-agnostic snapshot for any other tool to consume. Pass a filtered
result from [`query_runs`](@ref) to export just a selection. Returns the
path written.
"""
function export_json(entries::Vector{DatabaseEntry} = query_runs(); out::String = "tacobell_export.json")
    runs = map(entries) do e
        p = e.params
        fi = isempty(e.iterations) ? nothing : final_iteration(e)
        Dict{String, Any}(
            "id" => e.id, "model" => p.model, "symmetry" => p.symmetry, "algorithm" => p.algorithm,
            "chi" => p.chi, "K" => p.K, "mu0_sq" => p.mu0_sq, "lambda" => p.lambda,
            "n_iterations" => length(e.iterations),
            "central_charge" => (c = _best_c(e.iterations); isnan(c) ? missing : c),
            "final_iteration" => fi === nothing ? missing : fi.iteration,
            "sectors" => fi === nothing ? [] :
                [Dict{String, Any}("charge" => sec.charge, "dims" => sec.dims) for sec in fi.sectors],
        )
    end
    payload = Dict{String, Any}("generated_at" => string(now()), "runs" => runs)
    write(out, _to_json(payload))
    @info "Wrote $(length(entries)) run(s) to $out"
    return out
end

"""
    generate_catalog(; db_path, out) -> String

Write a Markdown table of every run in the database to `out`
(default `<db_path>/CATALOG.md`). This is what renders directly on GitHub
for anyone browsing the repository without Julia installed.

Call this after [`insert_run!`](@ref) or [`ingest_jld2!`](@ref) and commit
the resulting file alongside the new `db/runs/<uuid>.toml` and
`db/index.toml`. Returns the path written.
"""
function generate_catalog(; db_path::String = DB_PATH(), out::String = joinpath(db_path, "CATALOG.md"))
    idx = _load_index(db_path)
    rows = sort(idx, by = r -> (r["model"], r["symmetry"], r["algorithm"], r["chi"], r["K"], r["mu0_sq"], r["lambda"]))

    io = IOBuffer()
    println(io, "<!-- Auto-generated by `TACOBELL.generate_catalog()`. Do not edit by hand. -->")
    println(io)
    println(io, "# TACO-BELL catalog")
    println(io)
    println(io, "$(length(rows)) run(s), last generated $(now()).")
    println(io)
    println(io, "| ID | model | symmetry | algorithm | χ | K | μ₀² | λ | iterations | c |")
    println(io, "|---|---|---|---|---|---|---|---|---|---|")
    for row in rows
        c = isnan(row["central_charge"]) ? "—" : @sprintf("%.4f", row["central_charge"])
        println(io, @sprintf("| [%s](runs/%s.toml) | %s | %s | %s | %d | %d | %.4f | %.4f | %d | %s |",
            row["id"][1:8], row["id"], row["model"], row["symmetry"], row["algorithm"],
            row["chi"], row["K"], row["mu0_sq"], row["lambda"], row["n_iterations"], c))
    end

    mkpath(dirname(out))
    write(out, String(take!(io)))
    @info "Wrote catalog with $(length(rows)) entries to $out"
    return out
end

end # module

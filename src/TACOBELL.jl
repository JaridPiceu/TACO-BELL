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
export find_closest, find_sector, list_algorithms, list_symmetries, rebuild_index!
export ingest_jld2!, ingest_directory!, generate_catalog, export_csv
export central_charge_trajectory, plateau_estimate, plateau_central_charge

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

The scaling dimensions within one symmetry sector, labelled by
the TensorKit fusion-tree quantum numbers `(j, s)`.

- `twice_j` : twice the SU(2) / O(2) spin label `j`  (so `j = twice_j/2`)
- `s`       : second quantum number as stored (sector index, Z₂ charge, …)
- `dims`    : vector of scaling dimensions Δ in this sector, sorted ascending
"""
Base.@kwdef struct ScalingDimSector
    twice_j :: Int
    s       :: Int
    dims    :: Vector{Float64}
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
- `free_energy`     : free energy per site (optional)
- `correlation_len` : correlation length in lattice units (optional)
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
    free_energy     :: Union{Float64, Missing}              = missing
    correlation_len :: Union{Float64, Missing}              = missing
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

function Base.show(io::IO, s::ScalingDimSector)
    j_str = iseven(s.twice_j) ? string(s.twice_j ÷ 2) : string(s.twice_j, "/2")
    if isempty(s.dims)
        print(io, "ScalingDimSector(j=", j_str, ", s=", s.s, ", 0 dims)")
    else
        print(io, "ScalingDimSector(j=", j_str, ", s=", s.s, ", ", length(s.dims),
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
    Dict("twice_j" => sec.twice_j, "s" => sec.s, "dims" => sec.dims)
end

function _sector_from_dict(d)
    ScalingDimSector(twice_j=d["twice_j"], s=d["s"], dims=Float64.(d["dims"]))
end

function _results_to_dict(r::CFTResults)
    d = Dict{String,Any}(
        "iteration"     => r.iteration,
        "normalization" => r.normalization,
        "notes"         => r.notes,
        "sectors"       => [_sector_to_dict(s) for s in r.sectors],
    )
    r.central_charge  === missing || (d["central_charge"]  = r.central_charge)
    r.free_energy     === missing || (d["free_energy"]     = r.free_energy)
    r.correlation_len === missing || (d["correlation_len"] = r.correlation_len)
    d
end

function _results_from_dict(d)
    CFTResults(
        iteration       = d["iteration"],
        normalization   = d["normalization"],
        central_charge  = get(d, "central_charge",  missing),
        sectors         = [_sector_from_dict(s) for s in get(d, "sectors", [])],
        free_energy     = get(d, "free_energy",     missing),
        correlation_len = get(d, "correlation_len", missing),
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
    # Return the central charge at the last iteration that has one
    for r in reverse(iters)
        r.central_charge === missing || return r.central_charge
    end
    return NaN
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
    ingest_jld2!(filepath, params; db_path, allow_duplicate)

Read a JLD2 output file produced by TNRKit and insert it into the database.

`params` supplies the one field that is never stored inside the JLD2:
`model` (`"phi4_real"` / `"phi4_complex"`). Newer TNRKit output also stores
`symmetry` and `algorithm` directly, in which case those are read from the
file and whatever you pass in `params` for them is ignored — pass
placeholders for older files that predate this. Everything else (χ, K, μ₀²,
λ, iterations, normalization, central charge, scaling dimensions) is always
read from the file.

JLD2 files are HDF5 containers, so this reads them directly with HDF5.jl
rather than depending on JLD2.jl.

# Example
```julia
using TACOBELL

params = RunParameters(
    model     = "phi4_complex",
    symmetry  = "",   # overwritten from file if present
    algorithm = "",   # overwritten from file if present
    chi = 0, K = 0, mu0_sq = 0.0, lambda = 0.0,  # always overwritten from file
)
ingest_jld2!("Com_PD_O2_mu0-2_0_lam1_0_K8_chi16_iter20.jld2", params)
```
"""
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

# Parse the sector structure out of the JLD2 scaling_dimensions.structure field
function _parse_sectors(f, struct_ref, all_dims::Vector{Float64})
    struct_val = f[struct_ref][]
    kvvec_ref  = _field(struct_val, "kvvec")
    kvvec      = f[kvvec_ref][]          # Vector of object references, one per sector

    sectors = ScalingDimSector[]
    for ref in kvvec
        pair_val = f[ref][]              # named tuple: first=sector label, second=indices ref
        label      = _field(pair_val, "first")   # sector label is stored inline, not a Reference
        second_ref = _field(pair_val, "second")  # indices vector IS a Reference

        j_val   = _field(label, "j")     # also stored inline
        twice_j = Int(_field(j_val, "twice"))
        s       = Int(_field(label, "s"))

        indices   = Int.(f[second_ref][])   # 1-based into all_dims
        isempty(indices) && continue   # nothing found in this sector — same information as
                                        # "sector absent", so skip storing it (files have
                                        # dozens of these; they're most of the file size)
        dims = sort(all_dims[indices])

        push!(sectors, ScalingDimSector(twice_j=twice_j, s=s, dims=dims))
    end
    return sectors
end

"""
    ingest_directory!(dir; infer_params, db_path, extension=".jld2", allow_duplicate=false)

Recursively ingest every matching file under `dir` (hundreds of files is the
intended use case). `infer_params` is a function `filepath -> RunParameters`
supplying the fields not stored in the file (`model`, `symmetry`,
`algorithm`) for that particular file — pass a closure that returns the same
`RunParameters` for every call if a whole directory shares one setup, or
inspect `filepath` (e.g. with a regex on its filename or parent folder) if
different files need different labels.

A file that is already in the database is skipped and counted, not treated
as an error. A file that fails to read or parse is logged with `@warn` and
skipped, without aborting the rest of the batch. `db/CATALOG.md` is
regenerated once at the end if anything new was inserted.

Returns `(; total, inserted, skipped, failed)`.

# Example
```julia
using TACOBELL

summary = ingest_directory!("data/loop_tnr_runs";
    infer_params = _ -> RunParameters(
        model="phi4_complex", symmetry="O(2)", algorithm="LoopTNR",
        chi=0, K=0, mu0_sq=0.0, lambda=0.0,   # overwritten from each file
    ),
)
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
    find_sector(result, twice_j, s) -> Union{ScalingDimSector, Nothing}

Look up one symmetry sector's scaling dimensions from a single iteration's
result by its `(twice_j, s)` label, or `nothing` if that sector wasn't
populated at this iteration.
"""
function find_sector(result::CFTResults, twice_j::Int, s::Int)
    for sec in result.sectors
        sec.twice_j == twice_j && sec.s == s && return sec
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
    io = IOBuffer()
    println(io, _csv_row((
        "iteration", "normalization", "central_charge", "twice_j", "s", "dim_index", "delta",
    )))
    for r in sort(entry.iterations, by = x -> x.iteration)
        cc = r.central_charge === missing ? missing : r.central_charge
        if isempty(r.sectors)
            println(io, _csv_row((r.iteration, r.normalization, cc, missing, missing, missing, missing)))
        end
        for sec in r.sectors, (i, d) in enumerate(sec.dims)
            println(io, _csv_row((r.iteration, r.normalization, cc, sec.twice_j, sec.s, i, d)))
        end
    end
    write(out, String(take!(io)))
    @info "Wrote full trajectory for $(first(entry.id, 8))… to $out"
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

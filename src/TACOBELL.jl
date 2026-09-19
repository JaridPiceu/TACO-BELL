"""
    TACOBELL

Tensor Archive of Conformal Output — Best Ever Lattice Labour.

A lightweight, file-based database for CFT data (central charges, scaling
dimensions, ...) extracted from Tensor Network Renormalization (TNR)
calculations, designed for use alongside
[TNRKit](https://github.com/QuantumKitHub/TNRKit.jl/).

See the package README for a full walkthrough. The core entry points are
[`insert_run!`](@ref), [`query_runs`](@ref) and [`ingest_jld2!`](@ref)
(the latter requires `using HDF5`, see its docstring).
"""
module TACOBELL

using TOML, UUIDs, Dates, Printf

export RunParameters, ScalingDimSector, CFTResults, DatabaseEntry
export insert_run!, query_runs, load_entry, get_iteration, summarize_db
export find_closest, list_algorithms, list_symmetries, rebuild_index!
export ingest_jld2!

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
# Public API — JLD2 ingestion (implemented by the HDF5 package extension)
# ─────────────────────────────────────────────────────────────────────────────

"""
    ingest_jld2!(filepath, params; db_path, allow_duplicate)

Read a JLD2 output file produced by TNRKit and insert it into the database.

The `params` argument supplies fields that are *not* stored inside the JLD2
(model name, symmetry string, algorithm name). Everything else (χ, K, μ₀²,
λ, iterations, normalization, central charge, scaling dimensions) is read
from the file.

!!! note "Requires HDF5.jl"
    This method is provided by a package extension and only becomes
    available once you `using HDF5` (JLD2 files are HDF5 containers under
    the hood, so TACOBELL reads them directly with HDF5.jl to avoid a hard
    dependency on JLD2.jl). Run `using Pkg; Pkg.add("HDF5")` once if you
    don't already have it.

# Example
```julia
using TACOBELL, HDF5

params = RunParameters(
    model     = "phi4_complex",
    symmetry  = "O(2)",
    algorithm = "LoopTNR",
    # chi, K, mu0_sq, lambda are filled from the file:
    chi = 0, K = 0, mu0_sq = 0.0, lambda = 0.0,
)
ingest_jld2!("Com_PD_O2_mu0-2_0_lam1_0_K8_chi16_iter20.jld2", params)
```
"""
function ingest_jld2! end

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

end # module

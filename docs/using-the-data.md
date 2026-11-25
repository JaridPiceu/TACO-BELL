# Using the data in Julia

```julia
using TACOBELL
```

## Look runs up

```julia
runs = query_runs(symmetry="O(2)", chi=32)        # filter by any combination of fields
run  = find_closest(-1.9, 0.5; symmetry="O(2)", chi=32, algorithm="LoopTNR")  # nearest (μ₀², λ)
summarize_db()                                     # every run as a table
list_symmetries(); list_algorithms()
```

Filters: `model`, `symmetry`, `algorithm`, `chi`, `K`, `mu0_sq`, `lambda`
(the floats match within `mu0_sq_tol`/`lambda_tol`, default 1e-10).

## What a run contains

```julia-repl
julia> run
DatabaseEntry(6036132c…, phi4_complex/O(2)/LoopTNR, χ=32, K=12, μ₀²=-1.9, λ=0.5, 31 iters, c=1.0003)

julia> fi = final_iteration(run)
CFTResults(iteration=30, norm=0.6087, c=1.0002, 19/26 sectors populated)

julia> find_sector(fi; j=1, s=2)
ScalingDimSector(j=1, s=2, 27 dims, Δ∈[0.0271, 3.6466])

julia> fi.dims[1:3]          # all sectors together, sorted
3-element Vector{Float64}:
 -0.0
  0.0271…
  0.1083…
```

| Type | Fields |
|---|---|
| `DatabaseEntry` | `id`, `params`, `iterations`, `created_at`, `source_file`, `runtime_s` |
| `RunParameters` | `model`, `symmetry`, `algorithm`, `chi`, `K`, `mu0_sq`, `lambda` |
| `CFTResults` (one RG step) | `iteration`, `normalization`, `central_charge`, `sectors`, `notes`, and `.dims` (all Δ, sorted) |
| `ScalingDimSector` | `charge` (e.g. `Dict("j"=>1.0, "s"=>2.0)`), `dims` (sorted Δ) |

Sector labels use TensorKit's own field names, so any symmetry works:

| Symmetry | Look up a sector with |
|---|---|
| none (`Trivial`) | `find_sector(r)` |
| U(1) | `find_sector(r; charge=1)` |
| Z_N | `find_sector(r; n=1)` |
| SU(2) | `find_sector(r; j=0.5)` |
| O(2) | `find_sector(r; j=1, s=2)` |

Half-integers are stored as their real value (`0.5`, not doubled). A sector
in which nothing converged is either absent or has no dims, and
`find_sector` returns `nothing` if no sector matches.

## The converged value: plateau, not last iteration

The last RG step can drift at finite χ. Use the plateau estimate instead:
the flattest window of the flow.

```julia
plateau_central_charge(run)   # (value = 1.00035, err = 0.00014, range = 27:31, suspect = false)
plot_central_charge(run)      # the whole RG flow
```

`err` is the window's standard deviation. `suspect=true` means no clean
window was found. For other trajectories (e.g. one Δ), use
`plateau_estimate(values)`. Its docstring lists the tuning knobs. The
plateau value is what `summarize_db`, the catalog, CSV exports and the
website show.

## Export

```julia
export_csv(query_runs(symmetry="O(2)"); out="o2.csv")   # one row per run
export_csv(run; out="run.csv")                          # one run, every iteration/sector/Δ
export_json(query_runs(chi=32); out="chi32.json")       # parameters + plateau c
export_json(run; out="run.json")                        # one run in full
```

Prefer to skip Julia altogether? The database is plain HDF5; see
[Storage](storage.md#reading-it-without-julia).

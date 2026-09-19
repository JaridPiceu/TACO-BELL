# TACO-BELL

Tensor Archive of Conformal Output — Best Ever Lattice Labour

A lightweight, file-based database for CFT data (central charges, scaling
dimensions, ...) extracted from Tensor Network Renormalization (TNR)
calculations of the φ⁴ model in 2D.

Designed for use with [TNRKit](https://github.com/QuantumKitHub/TNRKit.jl/)
and intended to live on GitHub alongside your calculations, so that an
expensive run only ever needs to happen once.

## Why this design?

| Choice | Reason |
|--------|--------|
| Plain TOML files | Human-readable, git-diffable, no server or binary blobs |
| One file per run | Clean history; each commit = one new result; trivial to share via PR |
| Flat `index.toml` | Fast filtering without opening every run file |
| Direct JLD2 ingest | One function call to go from TNRKit output to database entry |
| HDF5 as a package extension | No hard dependency — `TACOBELL` loads instantly; `ingest_jld2!` becomes available the moment you `using HDF5` |

## Installation

Clone the repository and activate it as a Julia environment:

```julia
using Pkg
Pkg.develop(path="path/to/TACO-BELL")   # or Pkg.add(url="https://github.com/<you>/TACO-BELL")
using TACOBELL
```

or, working directly inside a checkout:

```julia
using Pkg
Pkg.activate("path/to/TACO-BELL")
Pkg.instantiate()

using TACOBELL
```

Only the Julia standard library (`TOML`, `UUIDs`, `Dates`, `Printf`) is
required for the core database. `HDF5.jl` is an optional dependency needed
only for [`ingest_jld2!`](#ingest-a-jld2-file-recommended-workflow) — see
below.

## Quick start

### Manual insertion

```julia
using TACOBELL

params = RunParameters(
    model="phi4_real", symmetry="Z2", algorithm="BTRG",
    chi=24, K=30, mu0_sq=-0.2, lambda=1.0,
)

iters = [
    CFTResults(
        iteration=20, normalization=0.693,
        central_charge=0.4998,
        sectors=[
            ScalingDimSector(twice_j=0, s=0, dims=[0.0, 1.0, 2.0]),
            ScalingDimSector(twice_j=2, s=1, dims=[0.125, 1.125]),
        ],
    ),
]

insert_run!(params, iters)
```

See [`scripts/insert_run_example.jl`](scripts/insert_run_example.jl) for a
runnable version.

### Ingest a JLD2 file (recommended workflow)

`ingest_jld2!` reads a TNRKit `.jld2` output directly with `HDF5.jl` (JLD2
files are HDF5 containers, so no `JLD2.jl` dependency is needed). This
method only exists once HDF5 is loaded:

```julia
using TACOBELL, HDF5

# Supply only what is not stored in the JLD2: model, symmetry, algorithm.
params = RunParameters(
    model     = "phi4_complex",
    symmetry  = "O(2)",
    algorithm = "LoopTNR",
    chi = 0, K = 0, mu0_sq = 0.0, lambda = 0.0,  # overwritten from file
)

entry = ingest_jld2!("Com_PD_O2_mu0-2_0_lam1_0_K8_chi16_iter20.jld2", params)
```

If `HDF5.jl` isn't installed yet, run `using Pkg; Pkg.add("HDF5")` once.
See [`scripts/ingest_jld2_example.jl`](scripts/ingest_jld2_example.jl) for a
runnable version:

```
julia --project=. scripts/ingest_jld2_example.jl path/to/your/file.jld2
```

`ingest_jld2!` currently expects the HDF5 key layout used by the reference
example file (`chi`, `K`, `μ0`, `λ`, `t`, and a `data` array of per-iteration
records with `central_charge` and `scaling_dimensions` grouped by the
`TensorKit` fusion-tree sector `(j, s)`). If your TNRKit output uses
different key names, adjust
[`ext/TACOBELLHDF5Ext.jl`](ext/TACOBELLHDF5Ext.jl) to match — that's the only
file that needs to know about the on-disk JLD2 structure.

### Query

```julia
# All O(2) runs with chi=16
runs = query_runs(symmetry="O(2)", chi=16)

# Overview table
summarize_db()

# A specific iteration from an entry
iter20 = get_iteration(runs[1], 20)
println(iter20.central_charge)
println(iter20.sectors)   # Vector{ScalingDimSector}

# Nearest parameter point to (μ₀²=-1.9, λ=1.0) with same setup
best = find_closest(-1.9, 1.0; symmetry="O(2)", chi=16, algorithm="LoopTNR")
```

### Bulk ingestion (hundreds of files)

`ingest_directory!` walks a directory recursively, ingests every matching
file, skips ones already in the database, and logs (without aborting) any
file that fails to read — built for pointing at a folder of hundreds of
TNRKit outputs at once:

```julia
using TACOBELL, HDF5

summary = ingest_directory!("data/loop_tnr_runs";
    infer_params = _ -> RunParameters(
        model="phi4_complex", symmetry="O(2)", algorithm="LoopTNR",
        chi=0, K=0, mu0_sq=0.0, lambda=0.0,   # overwritten from each file
    ),
)
# summary == (total=.., inserted=.., skipped=.., failed=..)
```

If a directory mixes several setups, `infer_params` can inspect the
filepath instead of returning a constant (e.g. a regex on the filename or
parent folder) to pick `model`/`symmetry`/`algorithm` per file. `db/CATALOG.md`
is regenerated automatically once at the end if anything new was inserted.

For hundreds of files, run this as a **script from a terminal** rather than
pasting into the REPL — see
[`scripts/bulk_ingest.jl`](scripts/bulk_ingest.jl), edit the `infer_params`
function near the top for your setup, then:

```
julia --project=. scripts/bulk_ingest.jl path/to/data/dir
```

It's safe to re-run on the same directory later (e.g. after adding more
files) — already-ingested files are skipped, not duplicated.

### Browse without Julia

After inserting new runs, regenerate the Markdown catalog so it's visible
directly on GitHub:

```julia
generate_catalog()   # writes db/CATALOG.md
```


## Parameter reference

### `RunParameters`

| Field | Type | Description |
|-------|------|-------------|
| `model` | String | `"phi4_real"` or `"phi4_complex"` |
| `symmetry` | String | Manifest symmetry in tensor: `"O(2)"`, `"U(1)"`, `"Z2"`, `"none"`, … |
| `algorithm` | String | `"LoopTNR"`, `"BTRG"`, `"GILT-TNR"`, … |
| `chi` | Int | Bond dimension χ |
| `K` | Int | Initial truncation |
| `mu0_sq` | Float64 | Bare mass squared μ₀² |
| `lambda` | Float64 | Quartic coupling λ |

### `CFTResults` (per iteration)

| Field | Type | Description |
|-------|------|-------------|
| `iteration` | Int | RG step index (0 = initial tensor) |
| `normalization` | Float64 | Per-site tensor norm at this step |
| `central_charge` | Float64? | Extracted central charge *c* |
| `sectors` | `Vector{ScalingDimSector}` | Scaling dims grouped by symmetry sector |
| `free_energy` | Float64? | Free energy per site (optional) |
| `correlation_len` | Float64? | Correlation length in lattice units (optional) |
| `notes` | String | Free-text annotation |

### `ScalingDimSector`

Each sector is labelled by the TensorKit fusion-tree quantum numbers `(j, s)`:

| Field | Type | Description |
|-------|------|-------------|
| `twice_j` | Int | `2j` (integer), so `j = twice_j/2` |
| `s` | Int | Second quantum number (Z₂ charge or sector index) |
| `dims` | `Vector{Float64}` | Scaling dimensions Δ in this sector, sorted ascending |

For O(2) symmetry as used in this example: `j` is the O(2) angular momentum
and `s ∈ {0,1,2}` distinguishes different irrep types
(singlet/vector/adjoint or similar).

## Database layout

```
db/
├── index.toml          ← lightweight summary of all runs (auto-managed)
├── CATALOG.md           ← human-readable table for browsing on GitHub (auto-managed)
└── runs/
    ├── <uuid>.toml     ← one file per run, all iterations inside
    └── …
```

The TOML format for a run uses TOML's `[[array of tables]]` syntax so each
iteration block is self-contained and human-readable. Never edit `index.toml`
or `CATALOG.md` by hand — both are generated from the `runs/*.toml` files.

`CATALOG.md` is what makes the database browsable without Julia: it's a
plain Markdown table that GitHub renders on the repo page, with each row
linking to its full `runs/<uuid>.toml` file. This repo is currently
**private**, so "browsable by anyone" means anyone with access to it —
switch the repo to public in GitHub's settings if you want it open to
everyone.

## Sharing results / contributing

Because every run is a single TOML file, contributing a new result is a
small, self-contained pull request:

1. Ingest your run locally (`insert_run!` or `ingest_jld2!`) so a new
   `db/runs/<uuid>.toml` and an updated `db/index.toml` are created.
2. Regenerate the catalog and commit all three files:

```julia
using TACOBELL
generate_catalog()
```

```
git add db/runs/<uuid>.toml db/index.toml db/CATALOG.md
git commit -m "Add run: <short description>"
```

3. If your branch was based on an older `main` and other runs were merged in
   the meantime, resolve any conflict in `db/index.toml` / `db/CATALOG.md`
   by regenerating both instead of hand-editing the diff:

```julia
using TACOBELL
rebuild_index!()
generate_catalog()
```

then re-commit those two files. Individual `db/runs/*.toml` files never need
manual edits or conflict resolution — they're independent by construction.

## Running the tests

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

The test suite exercises insertion, duplicate detection, querying,
`find_closest`, and index rebuilding against a temporary database, so it
never touches the real `db/` directory.

## Roadmap

- [ ] Export to CSV / HDF5 for plotting pipelines
- [ ] TNRKit callback hook: `insert_run!` triggered automatically at convergence
- [ ] Pluto notebook for interactive browsing and phase-diagram plotting
- [ ] `diff_runs(id1, id2)` to compare two parameter sets side by side

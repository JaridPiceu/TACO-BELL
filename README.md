# TACO-BELL

Tensor Archive of Conformal Output - Best Ever Lattice Labour


A lightweight, file-based database for CFT data extracted from
Tensor Network Renormalization (TNR) calculations of the φ⁴ model in 2D.

Designed for use with [TNRKit](https://github.com/QuantumKitHub/TNRKit.jl/) and intended to live on GitHub alongside your calculations.

## Why this design?

| Choice | Reason |
|--------|--------|
| Plain TOML files | Human-readable, git-diffable, no server or binary blobs |
| One file per run | Clean history; each commit = one new result; trivial to share via PR |
| Flat `index.toml` | Fast filtering without opening every run file |
| Direct JLD2 ingest | One function call to go from TNRKit output to database entry |
| Pure Julia stdlib + HDF5 | Minimal dependencies; HDF5.jl is already a TNRKit dependency |

## Quick start

```julia
using Pkg
Pkg.activate("path/to/TACOBELL")
Pkg.instantiate()

using TACOBELL
```

### Ingest a JLD2 file (recommended workflow)

```julia
# Supply only what is not stored in the JLD2: model, symmetry, algorithm.
params = RunParameters(
    model     = "phi4_complex",
    symmetry  = "O(2)",
    algorithm = "LoopTNR",
    chi = 0, K = 0, mu0_sq = 0.0, lambda = 0.0,  # overwritten from file
)

entry = ingest_jld2!("Com_PD_O2_mu0-2_0_lam1_0_K8_chi16_iter20.jld2", params)
```

### Manual insertion

```julia
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
└── runs/
    ├── <uuid>.toml     ← one file per run, all iterations inside
    └── …
```

The TOML format for a run uses TOML's `[[array of tables]]` syntax so each
iteration block is self-contained and human-readable.  Never edit `index.toml`
by hand.  If it gets out of sync after a merge, call `rebuild_index!()`.

## Sharing results / contributing

Because every run is a single TOML file, contributing a new result is a
one-file pull request.  After merging, call:

```julia
CFTDatabase.rebuild_index!()
```

to regenerate `index.toml` from all run files.

## Roadmap

- [ ] Export to CSV / HDF5 for plotting pipelines
- [ ] TNRKit callback hook: `insert_run!` triggered automatically at convergence
- [ ] Pluto notebook for interactive browsing and phase-diagram plotting
- [ ] `diff_runs(id1, id2)` to compare two parameter sets side by side

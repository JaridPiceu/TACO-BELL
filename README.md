# TACO-BELL

Tensor Archive of Conformal Output — Best Ever Lattice Labour

A database for CFT data (central charges, scaling dimensions, ...) computed
with Tensor Network Renormalization (TNR) for the φ⁴ model in 2D — built so
that an expensive calculation only ever has to be run once. Add a result,
look it up later by its physical parameters, and export it for anyone
(Julia or not) to use.

There's no server and nothing to install beyond Julia packages: every run
is just a plain text file in [`db/runs/`](db/runs/), tracked in this git
repo. Works for any manifest symmetry TensorKit supports (`Trivial`, `U(1)`,
`Z_N`, `SU(2)`, `O(2)`, ...) — sector labels are read generically, not
hardcoded to one symmetry. Designed for use with
[TNRKit](https://github.com/QuantumKitHub/TNRKit.jl/).

## Install

```julia
using Pkg
Pkg.develop(path="path/to/TACO-BELL")   # or Pkg.add(url="https://github.com/<you>/TACO-BELL")
using TACOBELL
```

Or, working directly inside a checkout of this repo:

```julia
using Pkg
Pkg.activate("path/to/TACO-BELL")
Pkg.instantiate()
using TACOBELL
```

## The three things you'll do

### 1. Add a result

From a TNRKit `.jld2` output file — the normal path. `model` is the one
thing you supply; everything else (χ, K, μ₀², λ, symmetry, algorithm, all
iteration data) is read straight from the file:

```julia
using TACOBELL
ingest_jld2!("Com_PD_O2_mu0-2_0_lam1_0_K8_chi16_iter20.jld2", "phi4_complex")
```

Don't have a compatible file yet? See
[Producing a compatible JLD2 file from TNRKit](#producing-a-compatible-jld2-file-from-tnrkit)
in the Reference section below — that's the piece that goes in your TNRKit
driver script, not in TACOBELL.

Got hundreds of files in a folder? Point `ingest_directory!` at it instead —
see [Bulk ingestion](#bulk-ingestion-hundreds-of-files) below. Building the
`CFTResults`/`RunParameters` yourself instead of from a file also works —
see [`scripts/insert_run_example.jl`](scripts/insert_run_example.jl).

### 2. Look it up

```julia
using TACOBELL

runs = query_runs(symmetry="O(2)", chi=32)     # filter by any combination of fields
best = find_closest(-1.9, 0.5; symmetry="O(2)", chi=32, algorithm="LoopTNR")  # nearest (μ₀², λ)
summarize_db()                                  # print everything as a table
```

Every value you get back prints as a short, readable summary instead of a
wall of nested numbers:

```julia-repl
julia> runs
16-element Vector{DatabaseEntry}:
 DatabaseEntry(6d40ca82…, phi4_complex/O(2)/LoopTNR, χ=32, K=12, μ₀²=-0.5, λ=0.5, 31 iters, c=-0.0000)
 DatabaseEntry(8061c8f0…, phi4_complex/O(2)/LoopTNR, χ=32, K=12, μ₀²=-0.6, λ=0.5, 31 iters, c=-0.0000)
 ⋮

julia> fi = final_iteration(runs[1])
CFTResults(iteration=30, norm=1.0000, c=-0.0000, 4/4 sectors populated)

julia> fi.sectors
4-element Vector{ScalingDimSector}:
 ScalingDimSector(j=0, s=0, 35 dims, Δ∈[-0.0000, 8.3959])
 ScalingDimSector(j=1, s=2, 31 dims, Δ∈[5.4419, 8.0624])
 ScalingDimSector(j=0, s=1, 25 dims, Δ∈[6.0167, 8.3648])
 ScalingDimSector(j=2, s=2, 17 dims, Δ∈[6.1429, 8.4344])

julia> find_sector(fi; j=1, s=2)
ScalingDimSector(j=1, s=2, 31 dims, Δ∈[5.4419, 8.0624])

julia> fi.dims[1:3]   # all sectors flattened into one sorted list, like TNRKit's own indexing
3-element Vector{Float64}:
 -0.0
  5.441933747711265
  5.4514434554069
```

The raw numbers are still all there underneath (`.central_charge`, `.dims`,
`sec.charge`, etc.) — this only changes how they *print*. A run's literal
last iteration can be numerically unstable (finite-χ truncation error
compounding under the RG flow), so for a more robust "converged" estimate,
run a window-based plateau search instead of trusting the last value
blindly:

```julia
plateau_central_charge(runs[1])
# (value = ..., err = ..., range = 22:26, suspect = false)
```

`err` is that window's standard deviation (a rough convergence error bar);
`suspect=true` means no clean window was found and the result may be
unreliable. See the [`plateau_estimate`](src/TACOBELL.jl) docstring for the
tunable knobs (`nwin`, `skipfrac`, `avoid_zero`).

To see the whole RG flow at a glance rather than just the endpoint:

```julia
using TACOBELL, Plots   # requires Plots.jl — see note below
plot_central_charge(runs[1])
```

`plot_central_charge` is provided by a package extension and only exists
once `Plots` is loaded — TACOBELL itself doesn't depend on it, since Plots
is a heavy stack and this is purely a convenience. The first time you use
it, `Pkg.add("Plots")` if you don't already have it; that adds Plots as an
ordinary dependency of your environment, same as any `Pkg.add`.

### 3. Get it out again

For a table you can open in Excel/pandas/R — a summary row per run, or the
full per-iteration/per-sector detail of one run:

```julia
export_csv(query_runs(symmetry="O(2)"); out="o2_results.csv")   # one row per run
export_csv(runs[1]; out="run_detail.csv")                        # everything in one run, long format
```

For browsing on GitHub without Julia at all, regenerate the Markdown table
after adding new runs:

```julia
generate_catalog()   # writes db/CATALOG.md
```

Or browse it interactively in a webpage — pick a symmetry/algorithm, type in
μ₀²/λ, and get the matched central charge and scaling-dimension sectors back,
with a small phase-diagram chart. See
[The Explorer webpage](#the-explorer-webpage) in the Reference section.

This repo is currently **private**, so "browsable by anyone" means anyone
with access to it — switch it to public in GitHub's settings if you want
that to mean the public internet.

---

## Reference

### Bulk ingestion (hundreds of files)

`ingest_directory!` walks a directory recursively, ingests every matching
file, skips ones already in the database, and logs (without aborting) any
file that fails to read:

```julia
using TACOBELL

summary = ingest_directory!("data/loop_tnr_runs"; infer_params = _ -> "phi4_complex")
# summary == (total=.., inserted=.., skipped=.., failed=..)
```

If a directory mixes several models, `infer_params` can inspect the
filepath instead of returning a constant (e.g. a regex on the filename or
parent folder). `db/CATALOG.md` regenerates automatically at the end if
anything new was inserted, and it's safe to re-run on the same directory
later — already-ingested files are skipped, not duplicated.

For hundreds of files, run this as a script from a terminal rather than
pasting into the REPL — edit
[`scripts/bulk_ingest.jl`](scripts/bulk_ingest.jl)'s `infer_params` for your
setup, then:

```
julia --project=. scripts/bulk_ingest.jl path/to/data/dir
```

### Producing a compatible JLD2 file from TNRKit

This is the piece that goes in **your TNRKit driver script** — the code
that runs the TNR scheme and saves the result in a shape `ingest_jld2!` can
read. TNRKit's own `finalize!` only returns the tensor's normalization at
each step, not any CFT data, so pair it with `CFTData(scheme)` via a custom
`Finalizer`:

```julia
using TNRKit, JLD2

# 1. A Finalizer that returns (normalization, CFTData) at every RG step.
#    Write one method per scheme type you use (LoopTNR here; same idea for
#    BTRG, TRG, ...).
function my_finalization(scheme::LoopTNR)
    n = finalize!(scheme)
    data = CFTData(scheme)
    return n, data
end
custom_Finalizer = Finalizer(my_finalization, Tuple{Float64, Any})

# 2. Build your scheme/tensor as usual (χ, K, μ₀², λ, symmetry all baked in
#    however your setup already does that), then run with that Finalizer.
#    `data` comes back as Vector{Tuple{Float64,CFTData}} — one entry per
#    RG step, exactly what `ingest_jld2!` expects.
t = @elapsed data = run!(scheme, trscheme, criterion, custom_Finalizer)

# 3. Save with the keys ingest_jld2! looks for. `symmetry`/`algorithm` are
#    optional but recommended (see "1. Add a result" above); `model` is
#    never stored here — you supply it later, at ingest time.
jldsave("Com_PD_O2_mu0$(μ0)_lam$(λ)_K$(K)_chi$(chi)_iter$(length(data) - 1).jld2";
    chi=chi, K=K, μ0=μ0, λ=λ, t=t, data=data,
    symmetry="O(2)", algorithm="LoopTNR", niter=length(data) - 1,
)
```

If your own driver script already saves under different key names, either
match these names in your `jldsave` call, or adjust the key names
`ingest_jld2!` looks for in [`src/TACOBELL.jl`](src/TACOBELL.jl) — either
works, but keeping the file itself standard means any TNRKit script anyone
writes later can feed straight into this database without modification.

### Windows path gotcha

Write file paths as a normal string: `raw"C:\Users\you\data\file.jld2"` or
`"C:/Users/you/data/file.jld2"`. `r"..."` is a *regex* literal in Julia
(unlike Python's raw strings) — using it for a Windows path errors on the
backslashes.

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

When ingesting from a JLD2 file, every field except `model` is a
placeholder that gets overwritten — see [`ingest_jld2!`](src/TACOBELL.jl).

### `CFTResults` (per iteration)

| Field | Type | Description |
|-------|------|-------------|
| `iteration` | Int | RG step index (0 = initial tensor) |
| `normalization` | Float64 | Per-site tensor norm at this step |
| `central_charge` | Float64? | Extracted central charge *c* |
| `sectors` | `Vector{ScalingDimSector}` | Scaling dims grouped by symmetry sector |
| `notes` | String | Free-text annotation |
| `.dims` *(derived)* | `Vector{Float64}` | All sectors' Δ flattened into one sorted list |

### `ScalingDimSector` — symmetry-agnostic sector labels

Each sector is labelled by a `charge::Dict{String,Float64}` holding
whatever quantum number(s) TensorKit's own charge/irrep type uses for that
symmetry — the same field names TensorKit itself uses, so this needs no
changes to support a new symmetry:

| Symmetry | TensorKit type | `charge` keys |
|----------|-----------------|---------------|
| none | `Trivial` | `{}` (empty) |
| U(1) | `U1Irrep` | `charge` |
| Z_N | `ZNIrrep{N}` | `n` |
| SU(2) | `SU2Irrep` | `j` |
| O(2) | `CU1Irrep` | `j`, `s` |

Half-integer values (spins, U(1)/O(2) charges) are stored as their real
value (e.g. `0.5`), not doubled. Look a sector up with
[`find_sector`](src/TACOBELL.jl), passing the same keyword names:

```julia
find_sector(fi; j=1, s=2)   # O(2)
find_sector(fi; charge=1)   # U(1)
find_sector(fi; n=1)        # Z_N
find_sector(fi)             # Trivial — no keywords
```

`dims::Vector{Float64}` holds the scaling dimensions Δ in that sector,
sorted ascending. Sectors where nothing converged aren't stored at all (an
absent sector means the same thing as an empty one, and TNR output can have
dozens of these per iteration).

### Database layout

```
db/
├── index.toml       ← lightweight summary of all runs (auto-managed)
├── CATALOG.md        ← human-readable table for browsing on GitHub (auto-managed)
└── runs/
    ├── <uuid>.toml  ← one file per run, all iterations inside
    └── …
```

Never edit `index.toml` or `CATALOG.md` by hand — both are generated from
the `runs/*.toml` files (`rebuild_index!()` / `generate_catalog()`).

### The Explorer webpage

[`web/index.html`](web/index.html) is a single static page (no server, no
build step) that lets you pick a model/symmetry/algorithm, type in μ₀²/λ,
and see the matched run's central charge and scaling-dimension sectors,
plus a small phase-diagram chart of *c* vs μ₀² for that setup. It reads
[`web/data/tacobell.json`](web/data/tacobell.json), a snapshot written by
`export_json`:

```julia
export_json(; out="web/data/tacobell.json")
```

Regenerate that snapshot (and commit it) the same way you'd regenerate
`CATALOG.md` — after adding new runs, before sharing the page. To view it
locally, **serve** the `web/` folder rather than double-clicking
`index.html` — browsers block a page's own `fetch()` of a local file opened
via `file://`. In VS Code, the simplest way is the "Live Server" extension:
right-click `web/index.html` → *Open with Live Server*. Without VS Code, any
static file server works, e.g. `python -m http.server 8000` from inside
`web/`, then open `http://localhost:8000`. To put this online for real, enable
GitHub Pages for this repo pointed at the `web/` folder — note that on
GitHub's free plan, Pages only serves **public** repos, so this waits for
the same "make it public" step as the rest of the database.

### Sharing results / contributing

Because every run is a single TOML file, contributing a new result is a
small, self-contained pull request:

1. Ingest your run locally so a new `db/runs/<uuid>.toml` and an updated
   `db/index.toml` are created.
2. Regenerate the catalog and commit all three files:

```julia
generate_catalog()
```
```
git add db/runs/<uuid>.toml db/index.toml db/CATALOG.md
git commit -m "Add run: <short description>"
```

3. If your branch was based on an older `main` and other runs were merged
   in the meantime, resolve any conflict by regenerating instead of
   hand-editing the diff, then re-commit:

```julia
rebuild_index!()
generate_catalog()
```

Individual `db/runs/*.toml` files never need manual edits or conflict
resolution — they're independent by construction.

### JLD2 key layout

`ingest_jld2!` reads TNRKit's `.jld2` output directly with `HDF5.jl` (JLD2
files are HDF5 containers, so no extra `JLD2.jl` dependency is needed). It
currently expects the key layout of the reference example files (`chi`,
`K`, `μ0`, `λ`, `t`, optionally `symmetry`/`algorithm`, and a `data` array
of per-iteration records with `central_charge` and `scaling_dimensions`).
Sector labels themselves are read generically via reflection (see
[ScalingDimSector](#scalingdimsector--symmetry-agnostic-sector-labels)
above), so a new *symmetry* needs no code change — but if your TNRKit
output uses different top-level key *names* entirely, adjust `ingest_jld2!`
in [`src/TACOBELL.jl`](src/TACOBELL.jl) to match; that's the only place
that needs to know about the on-disk structure.

### Running the tests

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

The test suite runs entirely against temporary databases, so it never
touches the real `db/` directory.

## Roadmap

- [ ] TNRKit callback hook: `insert_run!` triggered automatically at convergence
- [ ] GitHub Pages deployment for the Explorer webpage, once the repo goes public
- [ ] `diff_runs(id1, id2)` to compare two parameter sets side by side

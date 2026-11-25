# TACO-BELL

[![Test](https://github.com/JaridPiceu/TACO-BELL/actions/workflows/test.yml/badge.svg)](https://github.com/JaridPiceu/TACO-BELL/actions/workflows/test.yml)

**Tensor Archive of Conformal Output — Best Ever Lattice Labour**

A database of CFT data (central charges and scaling dimensions per symmetry
sector) computed with Tensor Network Renormalization for the 2D φ⁴ model, so
that every expensive calculation only has to be run once.

## 🌮 Browse the data: [jaridpiceu.github.io/TACO-BELL](https://jaridpiceu.github.io/TACO-BELL/)

The quickest way in, with nothing to install. Pick a model, symmetry and
algorithm, type in μ₀² and λ, and get the central charge, the scaling
dimensions per sector, the full RG flow and the phase diagram around that
point. → [How to use the website](docs/website.md)

For an overview of what has been computed, see
[`db/CATALOG.md`](db/CATALOG.md): every setup with its number of runs and
the λ/μ₀² it covers.

## Use it from Julia

```julia
using Pkg
Pkg.add(url="https://github.com/JaridPiceu/TACO-BELL.git")
# or, inside a clone:  Pkg.activate("path/to/TACO-BELL"); Pkg.instantiate()
```

```julia
using TACOBELL

run = find_closest(-1.9, 0.5; symmetry="O(2)", chi=32)  # nearest run in (μ₀², λ)
plateau_central_charge(run)                              # converged c ± error
find_sector(final_iteration(run); j=1, s=2)              # Δ's in one sector
runs = query_runs(symmetry="O(2)", chi=32)               # filter on any parameter
export_csv(runs; out="o2.csv")                           # for Excel / pandas / R
```

→ [Using the data](docs/using-the-data.md)

## Add your own runs

```
julia --project=. scripts/bulk_ingest.jl path/to/jld2/files
```

Then commit the new files in `db/shards/` and open a pull request. The
website updates itself once it's merged.
→ [Adding data](docs/adding-data.md), including how to save compatible
files from TNRKit

## Documentation

| Guide | For |
|---|---|
| [Using the website](docs/website.md) | browsing results without Julia |
| [Using the data](docs/using-the-data.md) | queries, sectors, plateau estimates, exports |
| [Adding data](docs/adding-data.md) | TNRKit output → database → pull request |
| [How the data is stored](docs/storage.md) | file format, precision, reading from Python, maintenance |

## Development

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

The tests run on every pull request (Julia 1.9 and latest) and must pass,
as the required `test` check, before merging into `main`. All workflow
actions are pinned to commit SHAs and kept up to date by Dependabot.

## License

- **Code**: [MIT](LICENSE)
- **Data** (`db/`): [CC BY 4.0](db/LICENSE). Reuse freely, crediting this repository.

# How the data is stored

```
db/
├── CATALOG.md                        overview per setup (generated)
├── LICENSE                           CC BY 4.0
└── shards/
    ├── 20261125-124711-54fd4e4d.h5   a batch of runs (≤ ~35 MB each)
    └── …
```

All runs live in compressed **HDF5 shard files**. Adding data always writes
*new* shards and never modifies existing ones. As a result, git stores each
run only once, data pull requests never conflict, and bulk ingest is fast.
The run index that `query_runs` and friends search is built from the shards
on first use and kept in memory. Nothing derived is stored on disk.

## Inside a shard

Each shard stores its data in columns:

| Group | Columns |
|---|---|
| `/runs` | `id`, `created_at`, `source_file`, `model`, `symmetry`, `algorithm`, `chi`, `K`, `mu0_sq`, `lambda`, `runtime_s`, `n_iterations` |
| `/iterations` | `iteration`, `normalization`, `central_charge`, `n_sectors` (+ `notes_index`/`notes_text` for non-empty notes) |
| `/sectors` | `charge` (index into `charge_labels`, e.g. `"j=1.0;s=2.0"`; `""` = trivial), `n_dims` |
| `/dims` | every scaling dimension, flat |

The counts link the levels together: a run owns the next `n_iterations`
iteration rows, an iteration owns the next `n_sectors` sector rows, and a
sector owns the next `n_dims` values of `/dims`. `missing` is stored as
`NaN`. The full specification is in [`src/storage.jl`](../src/storage.jl).

### Reading it without Julia

```python
import h5py, numpy as np
with h5py.File("db/shards/<file>.h5") as f:
    c    = f["iterations/central_charge"][:]
    n_it = f["runs/n_iterations"][:]
    first = np.concatenate([[0], np.cumsum(n_it)])   # run r → iteration rows first[r]:first[r+1]
```

## Precision

Scaling dimensions are stored as `Float32`. The error is at most 2·10⁻⁶
(measured over the whole database), far below what TNR resolves. This makes
the database ~2.5× smaller than `Float64`, whose last digits are noise that
doesn't compress. Central charges, normalizations and parameters are exact
`Float64`. For lossless storage, set `TACOBELL.DIMS_TYPE[] = Float64` before
ingesting. Each shard records its own type.

## Maintenance

| Task | Command |
|---|---|
| Merge many small shards (e.g. after `insert_run!` in a loop) | `compact!()` |
| Fix a name across runs | `correct_field!(:model, "old", "new")` |
| Regenerate `db/CATALOG.md` | `generate_catalog()` |
| Rebuild the website data locally | `julia --project=. scripts/update_web.jl` |
| Import runs in the pre-2026-11 one-TOML-file-per-run format | `import_toml_runs!("path/to/runs")` |

`compact!` and `correct_field!` rewrite existing shards, so git stores
those bytes again. Use them occasionally, not routinely.

## The website data

The website does not read the shards directly. On every push to `main` that
changes the data, the
[deploy workflow](../.github/workflows/deploy-pages.yml) runs
`scripts/update_web.jl`. That script writes `web/data/index.json` (all runs,
~1 MB) and one small `web/data/runs/<id>.bin` per run (~10 KB, with Δ
rounded to 4 decimals). Neither is committed.

To preview the site locally, run `scripts/update_web.jl`, then serve the
`web/` folder (e.g. `cd web && python -m http.server 8000`). Opening
`index.html` directly doesn't work, because browsers block `fetch` from
`file://` pages.

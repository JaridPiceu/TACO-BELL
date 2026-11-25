# Adding data

Runs come from [TNRKit](https://github.com/QuantumKitHub/TNRKit.jl/) `.jld2`
output files. Add them to your local copy of the database, then open a pull
request.

## 1. Save a compatible JLD2 file from TNRKit

Put this in **your TNRKit driver script**. TNRKit's own `finalize!` only
returns the normalization, so pair it with `CFTData` at every RG step:

```julia
using TNRKit, JLD2

function my_finalization(scheme::LoopTNR)   # one method per scheme type you use
    n = finalize!(scheme)
    return n, CFTData(scheme)
end
custom_Finalizer = Finalizer(my_finalization, Tuple{Float64, Any})

t = @elapsed data = run!(scheme, trscheme, criterion, custom_Finalizer)

jldsave("Com_PD_O2_mu0$(μ0)_lam$(λ)_K$(K)_chi$(chi)_iter$(length(data) - 1).jld2";
    chi=chi, K=K, μ0=μ0, λ=λ, t=t, data=data,
    symmetry="O(2)", algorithm="LoopTNR", niter=length(data) - 1)
```

These key names are what TACOBELL reads. `symmetry` and `algorithm` are
optional but recommended. The model name is never stored in the file; you
give it when ingesting.

## 2. Ingest

A whole folder (the usual case). Files already in the database are skipped,
and a broken file is logged without stopping the batch:

```
julia --project=. scripts/bulk_ingest.jl path/to/data/dir
```

Edit `infer_params` in [`scripts/bulk_ingest.jl`](../scripts/bulk_ingest.jl)
first if the folder isn't all `phi4_complex`. From Julia, the same thing is:

```julia
using TACOBELL
ingest_directory!("path/to/data/dir"; infer_params = _ -> "phi4_complex")
```

A single file:

```julia
ingest_jld2!("Com_PD_O2_mu0-2_0_lam1_0_K8_chi16_iter20.jld2", "phi4_complex")
```

Every ingest adds **new** files under `db/shards/` and never changes existing
ones. `ingest_directory!` also updates `db/CATALOG.md`; after other inserts,
run `generate_catalog()`.

**Results you computed some other way** can be added by building the
structs yourself; see
[`scripts/insert_run_example.jl`](../scripts/insert_run_example.jl). For
many runs, collect them and call `insert_runs!(entries)` once: each
`insert_run!` call writes its own small file. `compact!()` merges small
files afterwards.

### Old files (before TNRKit v0.5)

These store the data differently and never record the symmetry or
algorithm, so you must supply both:

```julia
ingest_jld2_old!("Com_CFT_μ0-0.01_λ0.01_K10truncrank16_niter15.jld2", "phi4_complex";
                 symmetry="U(1)", algorithm="LoopTNR")
```

For a folder, use [`scripts/bulk_ingest_old.jl`](../scripts/bulk_ingest_old.jl).
Not sure which symmetry? The sector labels give it away: `charge` → `U(1)`,
`n` → `Z2`/`ZN`, `j` alone → `SU(2)`, `j` + `s` → `O(2)`. These files have no
normalization, so it is `missing`.

## 3. Share it

```
git switch -c add-<description>
git add db/shards db/CATALOG.md
git commit -m "Add <n> runs: <description>"
git push -u origin add-<description>
```

Then open a pull request. `main` is protected: the tests must pass before
merging. Once it's merged, the [website](website.md) updates automatically.

Two pull requests that both add data never conflict on the data, because
each adds its own new files. If `db/CATALOG.md` conflicts, run
`generate_catalog()` after merging and commit the result.

## Fixing mistakes

A typo in a model, symmetry or algorithm name across many runs:

```julia
correct_field!(:model, "phi_complex", "phi4_complex")
```

## Tips

- On Windows, write paths as `raw"C:\data\file.jld2"` or `"C:/data/file.jld2"`.
  (`r"..."` is a regex in Julia, not a raw string.)
- If your JLD2 files use different key names, adjust `_read_jld2` in
  [`src/TACOBELL.jl`](../src/TACOBELL.jl). That's the only place that knows
  the file layout.

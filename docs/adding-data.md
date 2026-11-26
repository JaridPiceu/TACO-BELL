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

## 2. Add it to the database

Give a folder (searched recursively) or a single file, and the model name.
Everything else (χ, K, μ₀², λ, symmetry, algorithm) is read from each file:

```
julia --project=. scripts/add_data.jl path/to/data/dir phi4_complex
```

or, from Julia:

```julia
using TACOBELL
add_data!("path/to/data/dir", "phi4_complex")
```

Old and new TNRKit files are both recognized, even mixed in one folder.
Files whose name contains `_PARTIAL` (runs that stopped early) are added
too, tagged so they can be filtered out with `is_partial`.
Files already in the database are skipped, so it's safe to re-run on a
folder you've added to. A broken file is logged and skipped without stopping
the rest.

**Supplying or overriding a parameter.** Pass it as `key=value` to the
script, or as a keyword in Julia. It then applies to every file, replacing
the value stored in it (with a warning if the two differ):

```
julia --project=. scripts/add_data.jl path/to/old/runs phi4_complex algorithm=LoopTNR
julia --project=. scripts/add_data.jl path/to/runs phi4_real symmetry=Z2 algorithm=LoopTNR
```

Keys: `symmetry`, `algorithm`, `chi`, `K`, `mu0_sq`, `lambda`. If a file
doesn't record something it needs, `add_data!` stops at the first such file
and tells you which key to pass.

**A folder with several models:** in Julia, pass a function of the file path
instead of a model name:

```julia
add_data!("path/to/data", fp -> startswith(basename(fp), "Com_") ? "phi4_complex" : "phi4_real")
```

Every addition creates **new** files under `db/shards/` and never changes
existing ones. `db/CATALOG.md` is updated automatically. Before opening a
pull request you can check the result the same way CI does:

```
julia --project=. scripts/check_db.jl
```

**Results you computed some other way** can be added by building the
structs yourself; see
[`scripts/insert_run_example.jl`](../scripts/insert_run_example.jl). For
many runs, collect them and call `insert_runs!(entries)` once: each
`insert_run!` call writes its own small file. `compact!()` merges small
files afterwards.

### Old files (before TNRKit v0.5)

`add_data!` handles these too, but they never record the algorithm, so pass
`algorithm=...`. The symmetry is worked out from the sector labels
(`charge` → `U(1)`, `n` → `Z2`, `j` alone → `SU(2)`, `j` + `s` → `O(2)`).
Files without symmetry-labelled sectors (e.g. the "TransLine" runs) need
`symmetry=...` as well. These files have no normalization, so it is
`missing`.

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

`add_data!` warns when you use a model name the database doesn't have yet,
since a different spelling of an existing model splits its data in two. If
it happened anyway and the same runs now exist under both names,
`correct_field!` refuses (it would create duplicates). Delete one copy
instead:

```julia
delete_runs!(query_runs(model="Complex_phi4"))
```

`scripts/check_db.jl` reports such duplicates.

## Tips

- On Windows, write paths as `raw"C:\data\file.jld2"` or `"C:/data/file.jld2"`.
  (`r"..."` is a regex in Julia, not a raw string.)
- If your JLD2 files use different key names, adjust `_read_jld2` in
  [`src/TACOBELL.jl`](../src/TACOBELL.jl). That's the only place that knows
  the file layout.

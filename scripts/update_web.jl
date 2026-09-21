"""
update_web.jl

Regenerate everything the Explorer webpage (web/index.html) reads from the
database: the lightweight index and every run's per-iteration detail file.
Run this after adding new runs, before committing/redeploying the page.

Run from the repository root:
    julia --project=. scripts/update_web.jl
"""

using TACOBELL

index_path = export_json(; out=joinpath("web", "data", "tacobell.json"))
n = export_json_runs(; dir=joinpath("web", "data", "runs"))

println()
println("Updated $index_path and $n run detail file(s) in web/data/runs/.")
println("Commit and push (or redeploy) to publish the changes.")

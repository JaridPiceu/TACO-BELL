"""
update_web.jl

Regenerate everything the Explorer webpage (web/index.html) reads from the
database: web/data/index.json and one compact web/data/runs/<id>.bin per
run. The GitHub Pages workflow runs this on every deploy; run it yourself
only to preview the page locally (web/data/ is not committed).

Run from the repository root:
    julia --project=. scripts/update_web.jl
"""

using TACOBELL

dir = export_web(; dir = joinpath("web", "data"))
runs = filter(endswith(".bin"), readdir(joinpath(dir, "runs"); join = true))

println()
println("Wrote $(joinpath(dir, "index.json")) ($(round(filesize(joinpath(dir, "index.json")) / 1e6; digits = 2)) MB) and ",
        "$(length(runs)) run file(s) ($(round(sum(filesize, runs; init = 0) / 1e6; digits = 1)) MB) in $(joinpath(dir, "runs")).")

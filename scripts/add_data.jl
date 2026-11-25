"""
add_data.jl

Add TNRKit .jld2 output (one file, or a whole folder searched recursively)
to the database. Old (pre-TNRKit-v0.5) and current files are both handled,
even mixed in one folder. Everything except the model name is read from the
files; pass `key=value` to supply or override a parameter for every file.

Run from the repository root:

    julia --project=. scripts/add_data.jl <file-or-folder> <model> [key=value ...]

Examples:

    julia --project=. scripts/add_data.jl data/O2_runs phi4_complex
    julia --project=. scripts/add_data.jl data/old_runs phi4_complex algorithm=LoopTNR
    julia --project=. scripts/add_data.jl data/transline phi4_real symmetry=Z2 algorithm=LoopTNR

Keys: symmetry, algorithm, chi, K, mu0_sq, lambda (overrides), and
allow_duplicate=true, db_path=<folder> (default: db/). See `?add_data!` for details.
"""

using TACOBELL

length(ARGS) >= 2 || error("Usage: julia --project=. scripts/add_data.jl <file-or-folder> <model> [key=value ...]")
path, model = ARGS[1], ARGS[2]

const PARSERS = Dict(
    "symmetry" => identity, "algorithm" => identity,
    "chi" => s -> parse(Int, s), "K" => s -> parse(Int, s),
    "mu0_sq" => s -> parse(Float64, s), "lambda" => s -> parse(Float64, s),
    "allow_duplicate" => s -> parse(Bool, s), "db_path" => identity,
)

kwargs = Dict{Symbol,Any}()
for arg in ARGS[3:end]
    key, val = occursin('=', arg) ? split(arg, '='; limit = 2) : error("Expected key=value, got: $arg")
    haskey(PARSERS, key) || error("Unknown key '$key'. Valid keys: $(join(sort(collect(keys(PARSERS))), ", "))")
    kwargs[Symbol(key)] = PARSERS[key](String(val))
end

summary = try
    add_data!(path, model; kwargs...)
catch err
    err isa MissingParameterError || rethrow()
    println(stderr, "\n", err.msg, "\nOn the command line: add `key=value`, e.g. symmetry=Z2 or algorithm=LoopTNR.")
    exit(1)
end

println()
println("Done: $(summary.inserted) inserted, $(summary.skipped) already in the database, ",
        "$(summary.failed) failed (of $(summary.total) files).")
summary.inserted > 0 && println("Next: commit db/shards and db/CATALOG.md and open a pull request (see docs/adding-data.md).")

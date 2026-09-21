"""
bulk_ingest_old.jl

Ingest a directory of pre-TNRKit-v0.5 JLD2 output — the older, less
structured format (see `ingest_jld2_old!`'s docstring for what differs).
Same resilience contract as `bulk_ingest.jl`: safe to re-run, a bad file is
logged and skipped rather than aborting the batch.

Run from the repository root:
    julia --project=. scripts/bulk_ingest_old.jl path/to/old/data/dir
"""

using TACOBELL

isempty(ARGS) && error("Usage: julia --project=. scripts/bulk_ingest_old.jl <directory>")
data_dir = ARGS[1]
isdir(data_dir) || error("Not a directory: $data_dir")

# Old files never stored `symmetry`/`algorithm` at all, so — unlike the
# current-format script — you need to set those for real here too, not
# just `model`. `ingest_jld2_old!`/`ingest_directory_old!` require both
# explicitly (no default) for exactly this reason: it's easy to otherwise
# forget and silently insert blank-symmetry entries, which happened twice.
#
# Update this to match the batch you're ingesting — don't assume it's the
# same as last time. Check a sector's charge field name if unsure:
# "charge" -> U1Irrep -> "U(1)" | "n" -> ZNIrrep -> "Z2" | "j" -> SU2Irrep
# | "j"+"s" -> CU1Irrep -> "O(2)".
infer_params(filepath) = RunParameters(
    model="phi4_real", symmetry="Z2", algorithm="LoopTNR",
    chi=0, K=0, mu0_sq=0.0, lambda=0.0,  # overwritten from each file
)

summary = ingest_directory_old!(data_dir; infer_params)

println()
println("Done: $(summary.inserted) inserted, $(summary.skipped) skipped " *
        "(already in db), $(summary.failed) failed, out of $(summary.total) files.")
println()
summarize_db()

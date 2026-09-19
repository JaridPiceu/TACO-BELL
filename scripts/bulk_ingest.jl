"""
bulk_ingest.jl

Ingest every JLD2 file under a directory (recursively) into the database in
one go — built for hundreds of files at a time.

Safe to re-run: files already in the database are skipped and counted, and
a file that fails to read is logged with a warning and skipped rather than
aborting the whole batch.

Run from the repository root:
    julia --project=. scripts/bulk_ingest.jl path/to/data/dir
"""

using TACOBELL

isempty(ARGS) && error("Usage: julia --project=. scripts/bulk_ingest.jl <directory>")
data_dir = ARGS[1]
isdir(data_dir) || error("Not a directory: $data_dir")

# Every file under `data_dir` is assumed here to share the same
# model/symmetry/algorithm. If your directory mixes several setups, inspect
# `filepath` instead (e.g. a regex on `basename(filepath)` or on the
# subfolder it lives in) and return a different RunParameters per file.
infer_params(filepath) = RunParameters(
    model     = "phi4_complex",
    symmetry  = "O(2)",
    algorithm = "LoopTNR",
    chi = 0, K = 0, mu0_sq = 0.0, lambda = 0.0,  # overwritten from each file
)

summary = ingest_directory!(data_dir; infer_params)

println()
println("Done: $(summary.inserted) inserted, $(summary.skipped) skipped " *
        "(already in db), $(summary.failed) failed, out of $(summary.total) files.")
println()
summarize_db()

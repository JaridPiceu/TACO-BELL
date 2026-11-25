"""
migrate_toml_to_shards.jl

One-off migration from the old storage (one TOML file per run in db/runs/,
plus db/index.toml) to compressed HDF5 shards in db/shards/. Every run is
kept — ids, creation dates and all — including any that share a parameter
set. Afterwards every run is read back from the shards and compared against
its TOML source; the script errors out if anything differs beyond the
Float32 rounding of the scaling dimensions.

Run from the repository root (takes several minutes for ~18k runs):
    julia --project=. scripts/migrate_toml_to_shards.jl

It does not delete db/runs/ or db/index.toml; remove those yourself once
you're happy (`git rm -r --cached db/runs db/index.toml` + delete).
"""

using TACOBELL, TOML

db      = TACOBELL.DB_PATH()
run_dir = joinpath(db, "runs")
isdir(run_dir) || error("No $run_dir to migrate")
isempty(TACOBELL._shard_files(db)) || error("db/shards/ already has shard files; refusing to migrate twice")

files = filter(f -> endswith(f, ".toml"), readdir(run_dir; join = true))

# Group runs of the same setup (and λ slice) into the same shards, so e.g. a
# phase-diagram slice is read from one file. Uses the old index for the
# parameters, so the 2.9 GB of run files only have to be parsed once.
idx_path = joinpath(db, "index.toml")
if isfile(idx_path)
    rows = Dict(r["id"] => r for r in TOML.parsefile(idx_path)["entry"])
    sortkey(f) = (r = get(rows, first(splitext(basename(f))), nothing);
                  r === nothing ? ("~", "", "", 0, 0, 0.0, 0.0) :
                  (r["model"], r["symmetry"], r["algorithm"], r["chi"], r["K"], r["lambda"], r["mu0_sq"]))
    sort!(files; by = sortkey)
end

t = @elapsed n = import_toml_runs!(files; db_path = db, allow_duplicate = true)
@info "Imported $n of $(length(files)) runs in $(round(t; digits = 1)) s"
n == length(files) || error("imported $n runs but found $(length(files)) files")

# ── verification ─────────────────────────────────────────────────────────────
maxerr = 0.0; nchecked = 0
for shard in TACOBELL._shard_files(db)
    for e in TACOBELL._read_shard(shard)
        ref = TACOBELL._entry_from_dict(TOML.parsefile(joinpath(run_dir, e.id * ".toml")))
        ok = ref.id == e.id && ref.created_at == e.created_at && ref.params == e.params &&
             ref.source_file == e.source_file && isequal(ref.runtime_s, e.runtime_s) &&
             length(ref.iterations) == length(e.iterations)
        ok || error("run $(e.id): metadata differs after migration")
        for (a, b) in zip(ref.iterations, e.iterations)
            (a.iteration == b.iteration && isequal(a.normalization, b.normalization) &&
             isequal(a.central_charge, b.central_charge) && a.notes == b.notes &&
             length(a.sectors) == length(b.sectors)) || error("run $(e.id): iteration $(a.iteration) differs")
            for (s, r) in zip(a.sectors, b.sectors)
                (s.charge == r.charge && length(s.dims) == length(r.dims)) ||
                    error("run $(e.id): sector layout differs")
                isempty(s.dims) || (global maxerr = max(maxerr, maximum(abs.(s.dims .- r.dims))))
            end
        end
        global nchecked += 1
    end
end
nchecked == length(files) || error("verified $nchecked runs, expected $(length(files))")
maxerr < 1e-5 || error("scaling dimensions differ by up to $maxerr")

shards = TACOBELL._shard_files(db)
println()
println("Verified all $nchecked runs: metadata, central charges and normalizations identical,")
println("scaling dimensions within $(maxerr) (Float32 rounding).")
println("$(length(shards)) shard file(s), $(round(sum(filesize, shards) / 1e6; digits = 1)) MB total.")
generate_catalog(; db_path = db)

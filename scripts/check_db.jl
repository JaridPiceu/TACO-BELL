"""
check_db.jl

Integrity check of the committed database (db/): reads every run in every
shard and checks that the data is consistent and that db/CATALOG.md is up
to date. CI runs this on every pull request; run it yourself before opening
one that adds data.

Run from the repository root:
    julia --project=. scripts/check_db.jl [db_path]

Exits with status 1 if any check fails. Things that are allowed but worth
knowing (e.g. runs without scaling dimensions) are printed as notes.
"""

using TACOBELL, HDF5

const db = isempty(ARGS) ? TACOBELL.DB_PATH() : ARGS[1]
const problems = String[]
const notes = String[]
problem(msg) = (push!(problems, msg); println("  ✗ ", msg))
note(msg) = (push!(notes, msg); println("  · ", msg))

shards = TACOBELL._shard_files(db)
println("Checking $(length(shards)) shard(s) in $db")
isempty(shards) && problem("no shard files found")

leftovers = filter(endswith(".tmp"), readdir(joinpath(db, "shards")))
isempty(leftovers) || problem("half-written shard file(s) left behind: $(join(leftovers, ", "))")

seen_ids = Dict{String,String}()
seen_keys = Dict{Any,String}()
n_runs = 0; n_dims = 0
no_dims = Dict{Tuple{String,String,String},Int}()

for path in shards
    name = basename(path)
    size_mb = filesize(path) / 1e6
    size_mb < 95 || problem("$name is $(round(size_mb; digits=1)) MB, too close to GitHub's 100 MB file limit")
    fmt = h5open(f -> haskey(attributes(f), "tacobell_format") ? read(attributes(f)["tacobell_format"]) : nothing, path, "r")
    fmt == TACOBELL.SHARD_FORMAT || problem("$name: format $(repr(fmt)), expected $(TACOBELL.SHARD_FORMAT)")

    entries = try
        TACOBELL._read_shard(path)
    catch err
        problem("$name: unreadable ($(sprint(showerror, err)))"); continue
    end
    println("  $name: $(length(entries)) runs, $(round(size_mb; digits=1)) MB")

    for e in entries
        global n_runs += 1
        p = e.params
        what = "$name: run $(e.id) ($p)"
        haskey(seen_ids, e.id) && problem("$what: id also used in $(seen_ids[e.id])")
        seen_ids[e.id] = name
        k = TACOBELL._param_key(p)
        haskey(seen_keys, k) && problem("$what: same parameters as run $(seen_keys[k])")
        seen_keys[k] = e.id

        any(isempty, (p.model, p.symmetry, p.algorithm)) && problem("$what: empty model/symmetry/algorithm")
        (p.chi > 0 && p.K > 0) || problem("$what: χ and K must be positive")
        (isfinite(p.mu0_sq) && isfinite(p.lambda)) || problem("$what: non-finite μ₀² or λ")
        isempty(e.iterations) && problem("$what: no iterations")
        issorted([r.iteration for r in e.iterations]; lt = <=) || problem("$what: iterations not strictly increasing")

        nd = 0
        for r in e.iterations
            c = r.central_charge
            c === missing || isfinite(c) || problem("$what: non-finite central charge at iteration $(r.iteration)")
            for s in r.sectors
                all(isfinite, s.dims) || problem("$what: non-finite scaling dimension at iteration $(r.iteration)")
                nd += length(s.dims)
            end
        end
        global n_dims += nd
        nd == 0 && (no_dims[(p.model, p.symmetry, p.algorithm)] = get(no_dims, (p.model, p.symmetry, p.algorithm), 0) + 1)
    end
end

for ((m, s, a), n) in sort!(collect(no_dims))
    note("$n run(s) of $m / $s / $a have a central charge but no scaling dimensions")
end

# The committed catalog must match what the shards say now.
catalog = joinpath(db, "CATALOG.md")
if !isfile(catalog)
    problem("db/CATALOG.md is missing: run generate_catalog() and commit it")
else
    fresh = read(generate_catalog(; db_path = db, out = tempname() * ".md"), String)
    normalize(s) = replace(s, "\r\n" => "\n")
    normalize(fresh) == normalize(read(catalog, String)) ||
        problem("db/CATALOG.md is out of date: run generate_catalog() and commit the result")
end

println()
println("$n_runs runs, $n_dims scaling dimensions, $(length(problems)) problem(s), $(length(notes)) note(s).")
isempty(problems) || exit(1)

"""
ingest_jld2_example.jl

Shows how to ingest a TNRKit JLD2 output file into the database.
`model` is the only thing you need to supply manually — it's never stored
in the JLD2. `symmetry`/`algorithm` are read from the file when present
(recent TNRKit output has them); everything else (χ, K, μ₀², λ, all
iteration data, scaling dimensions for whatever symmetry was used) is
always read from the file.

Run from the repository root:
    julia --project=. scripts/ingest_jld2_example.jl path/to/your/file.jld2
"""

using TACOBELL

filepath = length(ARGS) >= 1 ? ARGS[1] :
    "Com_PD_O2_mu0-2_0_lam1_0_K8_chi16_iter20.jld2"

if !isfile(filepath)
    println("No JLD2 file found at: $filepath")
    println("Pass a path to a real TNRKit output file, e.g.:")
    println("    julia --project=. scripts/ingest_jld2_example.jl path/to/your/file.jld2")
    exit(1)
end

entry = ingest_jld2!(filepath, "phi4_complex")

println("Ingested entry: ", entry)
println()
summarize_db()

# ── Inspect the last stored iteration ─────────────────────────────────────────

last_iter = final_iteration(entry)
println("\nLast iteration: ", last_iter)

println("\n  Populated sectors (up to 5 dims shown each):")
for sec in last_iter.sectors
    println("    ", sec)
end

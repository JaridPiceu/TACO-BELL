"""
ingest_jld2_example.jl

Shows how to ingest a TNRKit JLD2 output file into the database.
`model` is the only thing you need to supply manually — it's never stored
in the JLD2. `symmetry`/`algorithm` are read from the file when present
(recent TNRKit output has them); everything else (χ, K, μ₀², λ, all
iteration data) is always read from the file.

Run from the repository root:
    julia --project=. scripts/ingest_jld2_example.jl path/to/your/file.jld2
"""

using Printf
using TACOBELL

filepath = length(ARGS) >= 1 ? ARGS[1] :
    "Com_PD_O2_mu0-2_0_lam1_0_K8_chi16_iter20.jld2"

if !isfile(filepath)
    println("No JLD2 file found at: $filepath")
    println("Pass a path to a real TNRKit output file, e.g.:")
    println("    julia --project=. scripts/ingest_jld2_example.jl path/to/your/file.jld2")
    exit(1)
end

# Only `model` truly needs a real value here; symmetry/algorithm/chi/K/
# mu0_sq/lambda are all placeholders — overwritten from the file.
params_template = RunParameters(
    model     = "phi4_complex",
    symmetry  = "",
    algorithm = "",
    chi = 0, K = 0, mu0_sq = 0.0, lambda = 0.0,
)

entry = ingest_jld2!(filepath, params_template)

println("Ingested $(length(entry.iterations)) iterations as entry $(entry.id[1:8])…")
println()
summarize_db()

# ── Inspect the last stored iteration ─────────────────────────────────────────

last_step = entry.iterations[end].iteration
last_iter = get_iteration(entry, last_step)
println("\nIteration $last_step:")
println("  normalization  = ", last_iter.normalization)
println("  central charge = ", last_iter.central_charge)
println("  number of sectors = ", length(last_iter.sectors))

# Print the lowest 5 scaling dims per sector (j, s)
println("\n  Lowest Δ per sector:")
for sec in sort(last_iter.sectors, by=s->(s.s, s.twice_j))
    isempty(sec.dims) && continue
    j_str = iseven(sec.twice_j) ? "$(sec.twice_j÷2)" : "$(sec.twice_j)/2"
    @printf("    j=%-4s  s=%d  →  Δ = %s\n",
            j_str, sec.s,
            join([@sprintf("%.4f", d) for d in sec.dims[1:min(5,end)]], "  "))
end

"""
insert_run_example.jl

Shows how to manually build a `RunParameters` / `CFTResults` pair and insert
it into the database, without going through JLD2 ingestion. Useful when you
computed results some other way, or want to hand-add a literature value.

Run from the repository root:
    julia --project=. scripts/insert_run_example.jl
"""

using TACOBELL

params = RunParameters(
    model="phi4_real", symmetry="Z2", algorithm="BTRG",
    chi=24, K=30, mu0_sq=-0.2, lambda=1.0,
)

iters = [
    CFTResults(
        iteration=20, normalization=0.693,
        central_charge=0.4998,
        sectors=[
            ScalingDimSector(twice_j=0, s=0, dims=[0.0, 1.0, 2.0]),
            ScalingDimSector(twice_j=2, s=1, dims=[0.125, 1.125]),
        ],
    ),
]

entry = insert_run!(params, iters)
println("Inserted entry $(entry.id[1:8])…")
println()
summarize_db()

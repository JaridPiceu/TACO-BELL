# Test fixtures

Small TNRKit output files (all real except the one marked synthetic), used by `test/runtests.jl` to check that
the JLD2 readers get the right numbers out of each file layout. Don't edit
them: the tests compare against values read from these exact files.

| file | layout | content |
|---|---|---|
| `new_format_Trivial.jld2` | TNRKit ≥ v0.5 | real φ⁴, no symmetry (`Trivial`), LoopTNR, χ=40, K=30, μ₀²=-0.72128, λ=0.5, 26 iterations |
| `new_format_O2.jld2` | TNRKit ≥ v0.5 | complex φ⁴, O(2), LoopTNR, χ=30, K=14, μ₀²=-0.02, λ=0.03, 26 iterations |
| `old_format_Z2.jld2` | before v0.5 | real φ⁴ "TransLine" run, χ=32, K=40, μ₀²=-0.9518, λ=0.7, 19 iterations, unsectored scaling dimensions (`sc_dim`) |
| `old_format_U1_synthetic.jld2` | before v0.5 | **synthetic**, made by `make_old_format_U1.jl`: 2 iterations with U(1)-labelled sectors (incl. a half-integer and an empty one), χ=16, K=10 |

To add one, pick the smallest file of the kind you need (a few hundred KB at
most) and add a test that pins a few of its values.

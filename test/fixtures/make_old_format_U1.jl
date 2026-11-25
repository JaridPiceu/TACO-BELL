# Generates old_format_U1_synthetic.jld2: a tiny file in the pre-TNRKit-v0.5
# layout with symmetry-labelled sectors (no such real file is small enough to
# commit). Each iteration is a Dict mapping "c" to the central charge and a
# U(1) irrep (a struct whose `charge` is a HalfInt, stored as `twice`) to
# that sector's scaling dimensions, all complex, as TNRKit wrote them.
#
# Not run by the tests; only needed to regenerate the file. Needs JLD2:
#     julia --project=@jld2tmp -e 'using Pkg; Pkg.add("JLD2")'
#     julia --project=@jld2tmp test/fixtures/make_old_format_U1.jl

using JLD2

struct HalfInt
    twice::Int
end
struct U1Irrep
    charge::HalfInt
end

iteration(c, sectors) = Dict{Any,Any}("c" => complex(c), (U1Irrep(HalfInt(2q)) => complex.(d) for (q, d) in sectors)...)

data = [
    iteration(1.02, [(0, [0.0, 1.31, 2.0]), (1, [0.27, 1.25]), (-1, [0.27, 1.25]), (0.5, Float64[])]),
    iteration(1.001, [(0, [0.0, 1.0, 2.0]), (1, [0.25, 1.25]), (-1, [0.25, 1.25]), (0.5, [0.0625])]),
]

jldsave(joinpath(@__DIR__, "old_format_U1_synthetic.jld2");
    ndimtrunc = 16, K = 10, μ0 = -0.01, λ = 0.01, t = 12.5, data = data)

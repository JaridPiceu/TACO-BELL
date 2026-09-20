"""
    TACOBELLPlotsExt

Package extension that implements [`TACOBELL.plot_central_charge`](@ref)
once the user loads Plots.jl (`using Plots`). Kept out of the main module
so that TACOBELL itself has no hard dependency on the (heavy) Plots.jl
stack for something that's purely a visualization convenience.
"""
module TACOBELLPlotsExt

using TACOBELL
using Plots

using TACOBELL: DatabaseEntry

function TACOBELL.plot_central_charge(entry::DatabaseEntry; kwargs...)
    sorted = sort(entry.iterations, by = r -> r.iteration)
    xs = Int[r.iteration for r in sorted if r.central_charge !== missing]
    ys = Float64[r.central_charge for r in sorted if r.central_charge !== missing]

    p = entry.params
    Plots.plot(xs, ys;
        marker = :circle, legend = false,
        xlabel = "RG iteration", ylabel = "central charge c",
        title = "$(p.model), $(p.symmetry), χ=$(p.chi), μ₀²=$(p.mu0_sq), λ=$(p.lambda)",
        titlefontsize = 10,
        kwargs...,
    )
end

end # module

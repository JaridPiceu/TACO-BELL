# Using the website

**<https://jaridpiceu.github.io/TACO-BELL/>** — browse every run in the
database from your browser. No Julia, no account, nothing to install.

## Find a result

1. **Pick a setup** in the top row:
   - **model**: `phi4_real` or `phi4_complex`, or whichever model names the
     data uses
   - **symmetry**: the symmetry enforced on the tensor (`Z2`, `U(1)`, `O(2)`, …)
   - **algorithm**: the TNR scheme (e.g. `LoopTNR`)
   - **χ** (bond dimension) and **K** (initial truncation). Leave these on
     *any* to use the most accurate run available (highest χ, then highest K).
2. **Type μ₀² and λ**, then press **Find** (or Enter).

The page shows the run at those parameters. If no run has exactly those
parameters, it shows the closest one and tells you how far off it is:

- **✓ exact match**: a run with exactly these μ₀² and λ exists.
- **≈ closest available point**: the nearest run in (μ₀², λ), with the
  difference shown next to it.

## Read the result

| Part of the page | What it shows |
|---|---|
| **Central charge** (big number) | *c* at the **plateau**: the average over the flattest stretch of the RG flow. This is the converged estimate, not the last iteration, which can drift at finite χ. |
| **Parameter chips** | The matched run's exact parameters. Hover a chip for an explanation. |
| **Iteration slider** | Step through the RG flow. It starts at the plateau; *jump to plateau* brings you back. Away from the plateau, *c* is that single iteration's value. |
| **Scaling dimensions by sector** | Δ at the selected iteration, grouped by symmetry sector (e.g. `j=1, s=2` for O(2)). Click *+N more* to see all of them. |
| **c vs iteration** | The central charge along the RG flow. The dashed line marks the selected iteration. |
| **Lowest scaling dimensions vs iteration** | The three lowest non-zero Δ along the RG flow. |
| **Phase diagram** | Plateau *c* against μ₀² for all runs with the same setup and λ, with the matched run in orange. Use *show as table* for the numbers. |

Hover over any point in a chart to read its exact value.

## Good to know

- Scaling dimensions on the website are rounded to 4 decimals. For full
  precision, use the Julia package (see [Using the data](using-the-data.md)).
- The site updates automatically whenever new data is merged into `main`
  (it takes a few minutes). The date at the bottom of the page shows when
  the data was last exported.
- Want the raw data? `data/index.json` lists every run, and
  `data/runs/<id>.bin` holds one run's trajectory. The format is described
  in [`src/webexport.jl`](../src/webexport.jl).

# Re-randomization-probability weighting (RPW): reproducible synthetic demonstration

Companion R code for:

> Rand K. *Re-randomization-probability weighting for population-adjusted comparisons of maintenance outcomes across trial designs.* Manuscript, 2026.

This archive contains the stand-alone R script `appendix_A_simulation.R` described in Appendix A of the manuscript. Run from top to bottom, it reproduces every number, table and figure in the manuscript's Synthetic demonstration section from fixed seeds. All data are synthetic. They represent no actual drug, trial or patient, and no real patient data are used or included.

## Contents

| File | Description |
| --- | --- |
| `appendix_A_simulation.R` | The complete, self-contained analysis script (MIT licence in the header) |
| `LICENSE` | MIT licence text |
| `README.md` | This file |

## Requirements

- R ≥ 4.0. The published results were produced with R 4.5.3 on Debian GNU/Linux 12 (x86_64).
- The analysis uses base R only (`stats`).
- The figures use [`ggplot2`](https://cran.r-project.org/package=ggplot2), which was version 4.0.3 for the published run.

Install the one dependency with:

```r
install.packages("ggplot2")
```

## Running

Place the script in an empty directory and run it from that directory:

```sh
Rscript appendix_A_simulation.R
```

Alternatively, `setwd()` to the script's directory and `source("appendix_A_simulation.R")` in an R session. The script writes all output relative to the working directory. A full run takes about 5 minutes on a single CPU core (measured with R 4.5.3 on x86_64 Linux).

## What the script does

The script is organized in numbered sections that match the comments in the code:

1. **Data-generating model.** Covariate distributions for three trial populations, an induction-response model in which severity is an effect modifier, and maintenance models.
2. **Three synthetic trials:**
   - **Trial A** (individual patient data): a randomized-withdrawal trial in which week-16 responders are re-randomized 2:2:1 to high dose, low dose and withdrawal. n = 800.
   - **Trial B** (aggregate data only): a randomized-withdrawal trial in which responders are re-randomized 1:2 to withdrawal and high dose. n = 700.
   - **Trial C** (aggregate data only): a treat-through trial in which everyone stays on high dose. n = 650.
3. **Baseline reconstruction.** RPW weights reconstruct the induction population's baseline covariate distribution within each maintenance arm (Table 2, Figure 2).
4. **MAIC weights.** Signorovitch method of moments, hand-coded.
5. **Estimators.** Horvitz–Thompson and Hájek forms of the induction-population estimands *q* (durable response), *m* (conditional maintenance) and *r* (induction response), plus the withdrawal-referenced contrast.
6. **True estimand values** by numerical integration over the known synthetic world. Section 6b adds the standardized anchoring-decomposition quantities quoted in the text.
7. **Primary analysis.** Unanchored and withdrawal-referenced comparisons of Trial A against Trials B and C, with a two-trial bootstrap (Table 3, Figure 3).
8. **Monte Carlo** over 2,000 replications of Trial A, checking that the adjusted point estimators are centred on the truth (Figure 4).
9. **Figures 1–4.**
10. **Unit checks** of the design-weight identities. These use `stopifnot()`, so the script stops if any check fails.
11. **Operating characteristics and stress scenarios** (Tables 4 and 5):
    - stratified assignment probabilities;
    - interval coverage of the paired two-trial bootstrap;
    - an omitted effect modifier;
    - correlated covariates with marginal-only matching.

Each stochastic stage has its own `set.seed()` call: 20260807 for data generation, then 1–10 for the later stages. Changing one stage therefore does not perturb the random-number stream of the others.

## Outputs

| Output | Content |
| --- | --- |
| `simulation_summary.csv` | Primary analysis estimates with bootstrap SEs and true values (Table 3) |
| `operating_characteristics_B.csv` | Bias, empirical SE, RMSE, mean bootstrap SE and coverage (Table 4) |
| `operating_characteristics_stress.csv` | Stress-scenario results (Table 5) |
| `figures/fig1_designs.{png,pdf}` | Figure 1, the trial designs |
| `figures/fig2_reconstruction.{png,pdf}` | Figure 2, baseline reconstruction |
| `figures/fig3_results.{png,pdf}` | Figure 3, primary results |
| `figures/fig4_montecarlo.{png,pdf}` | Figure 4, Monte Carlo distribution of the estimators |
| `sessionInfo.txt` | R session information for the run |

PNG figures are saved at 600 dpi on a white background. The remaining numbers quoted in the text, such as the anchoring-decomposition quantities and the replicate SD of the withdrawal-referenced estimator, are printed to the console.

## Reproducibility

The released script has been re-run in a clean directory. Its three CSV outputs were identical to those of the run reported in the manuscript. Results may differ in the last digits on other platforms or BLAS libraries, or with R versions that change the default random-number generator.

## Citation

If you use this code, please cite both the article and this archive:

> Rand K. *Re-randomization-probability weighting for population-adjusted comparisons of maintenance outcomes across trial designs.* 2026. [journal reference to be added]
>
> Rand K. *R code accompanying paper on re-randomization probability weighting for use in maintenance PAIC where trials include responder re-randomization* [R code]. Zenodo; 2026. https://doi.org/10.5281/zenodo.23240397

## Licence

MIT. See `LICENSE` and the script header. Copyright © 2026 Kim Rand.

## Contact

Kim Rand, Health Services Research Centre, Akershus University Hospital, Lørenskog, Norway and Maths in Health B.V. Klimmen, the Netherlands

# Multilevel network meta-regression for disconnected networks

R code to reproduce the analysis in the paper *Extending multilevel network
meta-regression to disconnected networks and single-arm studies*.

The example uses the plaque psoriasis data that come with the `multinma`
R package, so no extra data files are needed.

## Setup

The analysis needs the development version of `multinma` (it uses
`baseline_synthesis()` and `con()`, which are not yet on CRAN). Install it once:

```r
install.packages("remotes")
remotes::install_github("dmphillippo/multinma", ref = "develop")
```

Other packages used:

```r
install.packages(c("tidyverse", "ggtext", "patchwork", "ggraph", "loo"))
```

`multinma` uses Stan, so a working C++ toolchain is needed. See the
[rstan installation guide](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

## Running

Open `Paper Code Final.R` and run it from top to bottom (or `source()` it).
It fits all the models in the paper, so it takes a while.

## Outputs

- Network plots of the full, disconnected and external networks (shown, not saved)
- Bayesian R² and LOO predictive metrics for each model (printed)
- Plots saved to a `Plots/` folder:
  - `Probability_plot_PASI75.pdf`, `_PASI90.pdf`, `_PASI100.pdf`: predicted probabilities by bridge
  - `Probit_plot_PASI75.pdf`, `_PASI90.pdf`, `_PASI100.pdf`: marginal probit differences by bridge

PASI 75 is the main-text figure; PASI 90 and PASI 100 are in the appendix.

## Licence

MIT. See `LICENSE`.

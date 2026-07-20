# autoEnrich

Sweeps a numeric column (e.g. log fold-change) across a range of
thresholds, runs EnrichR pathway enrichment on the resulting gene set at
each threshold, and plots which pathways recur across multiple
thresholds.

## Installation

Once this package is on GitHub (see "Publishing" below), anyone can
install it with:

```r
# install.packages("devtools")
devtools::install_github("yourusername/autoEnrich")
```

For local development / before it's on GitHub, from inside this folder:

```r
# install.packages(c("devtools", "roxygen2"))
devtools::document()   # (re)generates NAMESPACE and man/ pages from roxygen comments
devtools::install()    # builds and installs the package into your R library
```

## Usage

```r
library(autoEnrich)

plot <- enrichR_summary(
  matrix             = my_matrix,
  gene_column        = "gene",
  threshold_column   = "avg_logFC",
  threshold_sequence = c(0, 5, 0.1),
  database           = "WikiPathways_2024_Mouse"
)

plot
```

See `?enrichR_summary` for full argument documentation once installed.

## Publishing so others can install it

1. Create a new GitHub repo (e.g. `autoEnrich`) and push this folder to it:

   ```bash
   cd autoEnrich
   git init
   git add .
   git commit -m "Initial commit"
   git branch -M main
   git remote add origin https://github.com/yourusername/autoEnrich.git
   git push -u origin main
   ```

2. Anyone with R can then run:

   ```r
   devtools::install_github("yourusername/autoEnrich")
   ```

   No CRAN submission needed for this to work — `install_github` builds
   straight from the repo.

3. (Optional) If you want it discoverable/installable via
   `install.packages("autoEnrich")` for anyone, that requires a CRAN
   submission, which has its own review process and stricter
   requirements (examples that run in under 5s, full test coverage,
   etc.). Most lab/internal tools stop at GitHub distribution.

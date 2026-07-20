# autoEnrich

Sweeps a numeric column (e.g. log fold-change) across a range of
thresholds, runs EnrichR pathway enrichment on the resulting gene set at
each threshold, and plots which pathways recur across multiple
thresholds. Useful for checking how stable an enrichment result is to
the exact cutoff used, rather than relying on a single arbitrary
threshold.

## Installation

```r
# install.packages("devtools")
devtools::install_github("gagraf/autoEnrich")
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

## Arguments

| Argument | Description |
|---|---|
| `matrix` | data.frame/tibble with a gene column and a numeric column to threshold on |
| `gene_column` | name of the column holding gene/feature IDs |
| `threshold_column` | name of the numeric column to threshold on |
| `threshold_sequence` | `c(start, end, by)` passed to `seq()` |
| `database` | single EnrichR database name, e.g. `"WikiPathways_2024_Mouse"` |
| `directions` | `"Positive"`, `"Negative"`, and/or `"None"` (magnitude only, ignores sign) |
| `threshold_type` | `"adj_pvalue_0.05"` (default) or `"pvalue_0.01"` |
| `min_genes` | minimum genes required at a threshold/direction to run enrichment (default 5) |
| `clean_gene_ids` | strip everything after `"_"` in gene IDs before enrichment (default TRUE) |

## Requirements

- Internet access (calls the EnrichR API)
- R packages: `enrichR`, `dplyr`, `tidyr`, `purrr`, `ggplot2` (installed automatically as dependencies)

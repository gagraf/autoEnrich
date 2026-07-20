#' Recurrent pathway enrichment across a threshold sweep
#'
#' Sweeps a numeric column across a range of thresholds, runs EnrichR on the
#' resulting gene set at each threshold (both up and down, or magnitude
#' only), and returns a single plot showing which pathways recur across
#' multiple thresholds.
#'
#' @param matrix data.frame/tibble containing a gene column and a numeric
#'   column to threshold on.
#' @param gene_column Name of the column holding gene/feature IDs.
#' @param threshold_column Name of the numeric column to threshold on.
#' @param threshold_sequence c(start, end, by) passed to seq() to build the
#'   thresholds swept over, e.g. c(0, 5, 0.1).
#' @param database Single EnrichR database name, e.g. "WikiPathways_2024_Mouse".
#' @param directions Which direction(s) to sweep. One or more of "Positive",
#'   "Negative", "None". "Positive"/"Negative" split on sign of value.
#'   "None" ignores sign entirely: thresholds are applied to abs(value),
#'   sweeping low to high, and genes on either side of zero are pooled into
#'   a single set per threshold. Useful for magnitude-only columns, or as a
#'   quick option for something like a p-value column (though for p-values,
#'   note that "None" still keeps the *largest* values as the threshold
#'   increases -- if smaller means more significant, transform to
#'   -log10(p) first and use "Positive" instead).
#' @param threshold_type Which EnrichR significance filter to apply before
#'   ranking terms: "adj_pvalue_0.05" (default) or "pvalue_0.01".
#' @param min_genes Minimum number of unique genes required at a given
#'   threshold/direction to bother running enrichment. Default 5.
#' @param clean_gene_ids If TRUE (default), strips everything from the
#'   first "_" onward in gene IDs before running enrichment.
#'
#' @details Per threshold/direction, only the top 20 EnrichR terms (by
#'   overlap fraction) are kept. Of the terms that recur across more than
#'   one threshold, only the top 25 per direction (by best overlap fraction
#'   seen) are plotted, using a fixed 25-color palette so lines stay
#'   distinguishable -- these caps are hardcoded rather than exposed as
#'   arguments, since more than 25 lines on one plot isn't readable anyway.
#'
#' @return A single ggplot object: overlap fraction vs threshold, one line
#'   per pathway that showed up at more than one threshold, faceted by
#'   direction.
#'
#' @examples
#' \dontrun{
#' plot <- enrichR_summary(
#'   matrix             = my_matrix,
#'   gene_column        = "gene",
#'   threshold_column   = "avg_logFC",
#'   threshold_sequence = c(0, 5, 0.1),
#'   database           = "WikiPathways_2024_Mouse"
#' )
#' plot
#' }
#'
#' @importFrom enrichR enrichr
#' @importFrom dplyr rename filter mutate group_by ungroup arrange
#'   summarise slice_head semi_join select bind_rows rowwise n
#' @importFrom tidyr expand_grid
#' @importFrom purrr pmap
#' @importFrom ggplot2 ggplot aes geom_line geom_point scale_color_manual
#'   scale_y_continuous facet_wrap theme_bw labs theme element_text
#' @importFrom magrittr %>%
#' @importFrom stats setNames
#' @export
enrichR_summary <- function(matrix,
                             gene_column,
                             threshold_column,
                             threshold_sequence = c(0, 2, 0.1),
                             database,
                             directions = c("Positive", "Negative"),
                             threshold_type = "adj_pvalue_0.05",
                             min_genes = 5,
                             clean_gene_ids = TRUE) {

  df <- as.data.frame(matrix) %>%
    dplyr::rename(feature = !!gene_column, value = !!threshold_column)

  if (clean_gene_ids) df$feature <- sub("_.*", "", df$feature)

  if (length(threshold_sequence) != 3) {
    stop("threshold_sequence must be c(start, end, by), e.g. c(0, 2, 0.1).")
  }
  thresholds <- seq(threshold_sequence[1], threshold_sequence[2], by = threshold_sequence[3])

  if (!all(directions %in% c("Positive", "Negative", "None"))) {
    stop("directions must be one or more of 'Positive', 'Negative', 'None'.")
  }
  if ("None" %in% directions && length(directions) > 1) {
    warning("'None' pools all genes by abs(value) regardless of sign; mixing it with ",
            "'Positive'/'Negative' in the same call will produce overlapping, ",
            "redundant gene sets. Consider using directions = 'None' on its own.")
  }

  # -------------------------------
  # 1. Build gene sets for every threshold x direction combo
  # -------------------------------

  grid <- tidyr::expand_grid(threshold = thresholds, direction = directions)

  gene_sets <- grid %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      genes = list(
        if (direction == "Positive")    unique(df$feature[df$value >= threshold])
        else if (direction == "Negative") unique(df$feature[df$value <= -threshold])
        else if (direction == "None")   unique(df$feature[abs(df$value) >= threshold])
        else stop("Unknown direction: ", direction, ". Use 'Positive', 'Negative', or 'None'.")
      ),
      n_genes = length(genes)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n_genes >= min_genes)

  if (nrow(gene_sets) == 0) {
    stop("No threshold/direction combination produced at least min_genes genes.")
  }

  # -------------------------------
  # 2. Run EnrichR at every threshold x direction combo
  # -------------------------------

  pval_col <- switch(
    threshold_type,
    "adj_pvalue_0.05" = "Adjusted.P.value",
    "pvalue_0.01"     = "P.value",
    stop("Unknown threshold_type: ", threshold_type)
  )
  pval_cutoff <- switch(threshold_type, "adj_pvalue_0.05" = 0.05, "pvalue_0.01" = 0.01)

  enrichment_results <- purrr::pmap(
    list(gene_sets$genes, gene_sets$direction, gene_sets$threshold),
    function(genes, direction, threshold) {

      enr <- enrichR::enrichr(genes, databases = database)
      if (is.null(enr[[database]]) || nrow(enr[[database]]) == 0) return(NULL)

      edf <- enr[[database]]
      edf$OverlapDecimal <- sapply(as.character(edf$Overlap), function(x) {
        nums <- as.numeric(strsplit(x, "/")[[1]])
        nums[1] / nums[2]
      })

      edf <- edf[edf[[pval_col]] < pval_cutoff, ]
      if (nrow(edf) == 0) return(NULL)

      edf <- edf[order(edf$OverlapDecimal, decreasing = TRUE), ]
      edf <- head(edf, 20)  # top terms per threshold/direction, by overlap fraction

      edf %>%
        dplyr::select(Term, OverlapDecimal) %>%
        dplyr::mutate(direction = direction, threshold = threshold)
    }
  )

  res_df <- dplyr::bind_rows(enrichment_results)

  if (nrow(res_df) == 0) {
    stop("No enrichment results passed the significance filter at any threshold.")
  }

  # -------------------------------
  # 3. Keep only pathways that recur across >1 threshold
  # -------------------------------

  recurrent_pathways <- res_df %>%
    dplyr::group_by(direction, Term) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::ungroup()

  if (nrow(recurrent_pathways) == 0) {
    stop("No pathways recurred across more than one threshold.")
  }

  top_terms <- recurrent_pathways %>%
    dplyr::group_by(direction, Term) %>%
    dplyr::summarise(max_overlap = max(OverlapDecimal, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(direction) %>%
    dplyr::arrange(dplyr::desc(max_overlap), .by_group = TRUE) %>%
    dplyr::slice_head(n = 25) %>%  # top recurring terms per direction, by best overlap seen
    dplyr::ungroup()

  recurrent_pathways <- recurrent_pathways %>%
    dplyr::semi_join(top_terms, by = c("direction", "Term"))

  # -------------------------------
  # 4. Plot
  # -------------------------------

  # Fixed qualitative palette so pathway lines stay visually distinct
  # (default ggplot hue scale gets muddy past ~8 categories).
  c25 <- c(
    "dodgerblue2", "#E31A1C", "green4", "#6A3D9A", "#FF7F00",
    "black", "gold1", "skyblue2", "#FB9A99", "palegreen2",
    "#CAB2D6", "#FDBF6F", "gray70", "khaki2", "maroon",
    "orchid1", "deeppink1", "blue1", "steelblue4", "darkturquoise",
    "green1", "yellow4", "yellow3", "darkorange4", "brown"
  )

  terms <- unique(recurrent_pathways$Term)
  term_colors <- stats::setNames(rep_len(c25, length(terms)), terms)

  ggplot2::ggplot(recurrent_pathways,
                   ggplot2::aes(x = threshold, y = OverlapDecimal, color = Term, group = Term)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = term_colors) +
    ggplot2::scale_y_continuous(limits = c(0, 1)) +
    ggplot2::facet_wrap(~ direction) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      x = "Threshold", y = "Overlap",
      title = "Recurrent Pathways Across Thresholds"
    ) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 6),
                    legend.title = ggplot2::element_text(size = 8))
}

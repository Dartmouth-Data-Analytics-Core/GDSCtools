#' QC Scatter Plot with Marginal Densities
#'
#' Generates a two-panel QC scatter plot for any pair of metadata variables.
#' The top panel shows all samples overlaid with marginal density curves
#' (via [ggExtra::ggMarginal()]). The bottom panel facets by the grouping
#' variable. Axes can optionally be log10-scaled. This function assumes
#' you have extracted the metadata from your Seurat object and saved as a
#' dataframe.
#'
#' @param metadata A data frame of per-cell or per-sample QC metrics.
#' @param Xvar Character string. Column name in `metadata` for the x-axis.
#' @param Yvar Character string. Column name in `metadata` for the y-axis.
#' @param groupVar Character string. Column name in `metadata` used for color
#'   and faceting (default `"orig.ident"`).
#' @param logTransformX Logical. Apply [ggplot2::scale_x_log10()] (default
#'   `FALSE`).
#' @param logTransformY Logical. Apply [ggplot2::scale_y_log10()] (default
#'   `FALSE`).
#' @param colors Optional named or positional color vector for groups (default
#'    `NULL`, in which case, 8 default colors are used.)
#' @param outDir Optional directory path for saving the plot. If `NULL`
#'   (default), the plot is not saved.
#' @param width Numeric. Plot width in inches (default `12`).
#' @param height Numeric. Plot height in inches (default `10`).
#'
#' @return A combined [cowplot::plot_grid()] object (returned invisibly when
#'   saved to disk).
#'
#' @examples
#' \dontrun{
#' qcScatter(metadata,
#'   Xvar = "log_atac_peak_region_fragments",
#'   Yvar = "pct_reads_in_peaks",
#'   groupVar = "orig.ident",
#'   colors = my_colors
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_classic labs scale_color_manual
#'   theme element_text scale_x_log10 scale_y_log10 facet_grid ggsave
#' @importFrom ggExtra ggMarginal
#' @importFrom cowplot plot_grid
#' @import rlang
#' @export

qcScatter <- function(metadata, Xvar, Yvar, groupVar = "orig.ident",
                      logTransformX = FALSE, logTransformY = FALSE,
                      colors = NULL, outDir = NULL, width, height) {

  stopifnot(
    Xvar %in% colnames(metadata),
    Yvar %in% colnames(metadata),
    "orig.ident" %in% colnames(metadata)
  )

  n_groups <- length(unique(metadata[["orig.ident"]]))

  if (is.null(colors)) {
    default_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728",
                         "#9467BD", "#8C564B", "#E377C2", "#7F7F7F")
    if (n_groups > length(default_palette)) {
      stop("More than ", length(default_palette), " groups detected. ",
           "Please provide a color vector via the `colors` argument.")
    }
    colors <- default_palette[seq_len(n_groups)]
  }

  # Convert `variables` to a symbol
  variableX <- ensym(Xvar)
  variableY <- ensym(Yvar)
  facetter <- ensym(groupVar)

  # Handle dynamic axis labels
  x_label <- if (logTransformX) paste0("Log(", gsub("_", " ", as_string(variableX)), ")") else gsub("_", " ", as_string(variableX))
  y_label <- if (logTransformY) paste0("Log(", gsub("_", " ", as_string(variableY)), ")") else gsub("_", " ", as_string(variableY))
  title <- paste0(x_label, " vs ", y_label)

  # If Y var is logged
  p1 <- ggplot2::ggplot(metadata, aes(x = .data[[as_string(variableX)]], y = .data[[as_string(variableY)]], color = orig.ident)) +
    geom_point(size = 2, alpha = 0.3) +
    theme_classic(base_size = 16) +
    labs(x = x_label,
         y = y_label,
         color = "Sample",
         title = title) +
    scale_color_manual(values = colors) +
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"),
          legend.position = "none")

  # Conditionally apply log scales
  if (logTransformX) p1 <- p1 + scale_x_log10()
  if (logTransformY) p1 <- p1 + scale_y_log10()

  #----- Add margin annotations of distributions
  p2 <- ggExtra::ggMarginal(p1, groupColour = TRUE, groupFill = FALSE)

  p3 <- ggplot2::ggplot(metadata, aes(x = .data[[as_string(variableX)]], y = .data[[as_string(variableY)]], color = orig.ident)) +
    geom_point(size = 2, alpha = 0.3) +
    theme_classic(base_size = 16) +
    labs(x = x_label,
         y = y_label,
         color = "Sample") +
    scale_color_manual(values = colors) +
    facet_grid(~.data[[as_string(facetter)]]) +
    theme(axis.title = element_text(face = "bold"),
          axis.text = element_text(face = "bold"),
          strip.text = element_text(face = "bold"),
          legend.position = "none")

  # Apply log scales to faceted plot if needed
  if (logTransformX) p3 <- p3 + scale_x_log10()
  if (logTransformY) p3 <- p3 + scale_y_log10()

  final <- cowplot::plot_grid(p2, p3, ncol = 1)
  if (!is.null(outDir)) {
    dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(
      file.path(outDir, paste0(Xvar, "_vs_", Yvar, ".png")),
      final, width = width, height = height
    )
    return(invisible(final))
  }

  final
}

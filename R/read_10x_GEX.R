#' Read a 10x Genomics Gene Expression Sample into a Seurat Object
#'
#' Reads a single sample's filtered feature barcode matrix H5 file from a 10x
#' Genomics Cell Ranger output directory and creates a
#' [Seurat][Seurat::CreateSeuratObject] object. A log10 genes-per-UMI
#' complexity metric is added to the object metadata as `log10GenesPerUMI`.
#'
#' The expected input path is
#' `<sampleDir>/<sample>/outs/filtered_feature_bc_matrix.h5`. To process
#' multiple samples, iterate over sample names with [purrr::map()] or
#' [lapply()].
#'
#' @param sampleDir Character string. Path to the parent directory that contains
#'   one subdirectory per sample.
#' @param sample Character string. Name of the sample subdirectory within
#'   `sampleDir`.
#' @param ident Character string. Identity label to assign to every cell in the
#'   resulting Seurat object (stored in the `ident` metadata column).
#' @param min.cells Integer. Include features detected in at least this many
#'   cells. Passed to [Seurat::CreateSeuratObject()].
#' @param min.features Integer. Include cells where at least this many features
#'   are detected. Passed to [Seurat::CreateSeuratObject()].
#'
#' @return A [Seurat][Seurat::CreateSeuratObject] object for the specified
#'   sample with the `RNA` assay, an `ident` metadata column, and a
#'   `log10GenesPerUMI` complexity metric.
#'
#' @examples
#' \dontrun{
#' library(purrr)
#'
#' sample_dir <- "/path/to/data/"
#' samples    <- c("S1", "S2", "S3")
#' idents     <- c("CondA", "CondB", "CondC")
#'
#' seurat_list <- map(samples, idents, ~ read_10x_GEX(
#'   sampleDir    = sample_dir,
#'   sample       = .x,
#'   ident        = .y,
#'   min.cells    = 10,
#'   min.features = 100
#' ))
#' }
#'
#' @importFrom Seurat Read10X_h5 CreateSeuratObject
#' @importFrom purrr map
#' @export

read_10x_GEX <- function(sampleDir, sample, ident, min.cells, min.features) {
  print(sample)
  counts <- Read10X_h5(filename = paste0(sampleDir, sample, "/outs/filtered_feature_bc_matrix.h5"))
  x <- CreateSeuratObject(counts = counts$`Gene Expression`,
                          min.cells = min.cells,
                          min.features = min.features,
                          assay = "RNA")
  x[["ident"]] <- paste0(ident)

  #----- Add complexity measure
  x$log10GenesPerUMI <- log10(x$nFeature_RNA) / log10(x$nCount_RNA)
  return(x)
}

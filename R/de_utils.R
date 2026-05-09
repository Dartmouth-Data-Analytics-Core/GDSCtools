#' Run Cluster-Specific Differential Expression Analysis
#'
#' Performs differential expression between two sample identities within each
#' specified cell-type cluster using \code{Seurat::FindMarkers}. Results are
#' filtered by adjusted p-value and absolute log2 fold-change thresholds, and
#' annotated with a regulation-direction label.
#'
#' @param obj A \code{Seurat} object containing the RNA assay and metadata
#'   columns specified by \code{cluster_col} and \code{ident_col}.
#' @param clusters Character vector of cluster names (values in the
#'   \code{cluster_col} metadata column) to iterate over.
#' @param cluster_col Character string naming the metadata column that holds
#'   cell-type cluster assignments. Default \code{"final_celltype"}.
#' @param ident_col Character string naming the metadata column to set as the
#'   active identity for the DE comparison (e.g., \code{"sample"},
#'   \code{"condition"}). Default \code{"sample"}.
#' @param ident1 Character string. The first identity level (numerator in the
#'   log fold-change calculation). Positive fold-changes indicate genes
#'   up-regulated in \code{ident1}.
#' @param ident2 Character string. The second identity level (denominator).
#' @param min.pct Numeric. Minimum fraction of cells in either group that must
#'   express a gene for it to be tested. Passed to
#'   \code{Seurat::FindMarkers}. Default \code{0.25}.
#' @param padj_thresh Numeric. Maximum adjusted p-value for a gene to be
#'   retained in the filtered output. Default \code{0.05}.
#' @param logfc_thresh Numeric. Minimum absolute log2 fold-change for a gene
#'   to be retained. Default \code{0.5}.
#'
#' @return A named list of data frames (one per cluster). Each data frame
#'   contains the standard \code{FindMarkers} columns plus:
#'   \describe{
#'     \item{gene}{Gene symbol (from row names).}
#'     \item{cluster}{The cluster label.}
#'     \item{comparison}{String \code{"<ident1>_vs_<ident2>"}.}
#'     \item{res}{Direction label: \code{"Upregulated in <ident1>"} or
#'       \code{"Upregulated in <ident2>"}.}
#'   }
#'
#' @details The function subsets the Seurat object to each cluster, sets the
#'   active identity to \code{ident_col}, and runs
#'   \code{Seurat::FindMarkers(assay = "RNA")}. A frequency table of
#'   regulation direction is printed for each cluster.
#'
#' @export
run_cluster_de <- function(
    obj,
    clusters,
    cluster_col = "final_celltype",
    ident_col = "sample",
    ident1,
    ident2,
    min.pct = 0.25,
    padj_thresh = 0.05,
    logfc_thresh = 0.5
) {

  results_list <- lapply(clusters, function(clust) {

    message(paste0("Running DE for ", clust, " (", ident1, " vs ", ident2, ")"))

    sub <- obj[, obj@meta.data[[cluster_col]] == clust]

    Seurat::Idents(sub) <- ident_col

    markers <- Seurat::FindMarkers(
      sub,
      assay = "RNA",
      ident.1 = ident1,
      ident.2 = ident2,
      min.pct = min.pct
    )

    markers$gene <- rownames(markers)
    markers$cluster <- clust
    markers$comparison <- paste0(ident1, "_vs_", ident2)

    markers.filt <- subset(
      markers,
      p_val_adj <= padj_thresh & abs(avg_log2FC) > logfc_thresh
    )

    markers.filt$res <- ifelse(
      markers.filt$avg_log2FC >= logfc_thresh,
      paste0("Upregulated in ", ident1),
      paste0("Upregulated in ", ident2)
    )

    print(table(markers.filt$res))

    return(markers.filt)
  })

  names(results_list) <- clusters
  return(results_list)
}

#' Write Differential Expression or GO-ORA Results to CSV Files
#'
#' Saves a named list of DE or GO enrichment results into cluster-specific CSV
#' files. For DE results, each list element is written directly. For GO results,
#' the \code{GeneRatio} column is parsed into numeric numerator/denominator
#' columns, and results can optionally be split by regulation direction.
#'
#' @param de_list Named list of data frames (for \code{type = "DE"}) or
#'   \code{enrichResult} objects / direction-split lists (for
#'   \code{type = "GO"}).
#' @param outSubDir Character path to the output directory. Created
#'   recursively if it does not exist.
#' @param type Character, one of \code{"DE"} or \code{"GO"}. Controls output
#'   filename suffixes and parsing logic.
#' @param split_direction Logical or \code{NULL}. Required when
#'   \code{type = "GO"}. If \code{TRUE}, each element of \code{de_list} is
#'   expected to be a list with sub-elements \code{$up} and \code{$down}
#'   (each an \code{enrichResult}). If \code{FALSE}, elements are single
#'   \code{enrichResult} objects.
#'
#' @return \code{invisible(NULL)}. Called for its side effect of writing CSV
#'   files.
#'
#' @section Output files:
#' \itemize{
#'   \item DE mode: \code{<outSubDir>/<cluster>_DE_Results.csv}
#'   \item GO mode (no split): \code{<outSubDir>/<cluster>_GO_ORA_Results.csv}
#'   \item GO mode (split): \code{<outSubDir>/<cluster>_GO_ORA_up_Results.csv}
#'     and \code{..._down_Results.csv}
#' }
#'
#' @export
write_de_results <- function(de_list, outSubDir, type = c("DE", "GO"),
                             split_direction = NULL) {

  type <- match.arg(type)

  if (type == "GO" && is.null(split_direction)) {
    stop("'split_direction' must be specified when type = 'GO'")
  }

  if (!dir.exists(outSubDir)) {
    dir.create(outSubDir, recursive = TRUE)
  }

  .write_go <- function(ego, celltypeName, direction_label = NULL) {
    if (is.null(ego)) return(invisible(NULL))
    res <- as.data.frame(ego)
    if (nrow(res) == 0) {
      message("Skipping ", celltypeName, " ", direction_label %||% "", " (no results)")
      return(invisible(NULL))
    }
    tmp <- strsplit(as.character(res$GeneRatio), "/")
    res$geneRatio_num <- as.numeric(sapply(tmp, `[`, 1))
    res$geneRatio_den <- as.numeric(sapply(tmp, `[`, 2))

    suffix <- if (!is.null(direction_label)) {
      paste0("_GO_ORA_", direction_label, "_Results.csv")
    } else {
      "_GO_ORA_Results.csv"
    }

    outfile <- paste0(outSubDir, celltypeName, suffix)
    utils::write.csv(res, file = outfile, row.names = FALSE)
    message("Wrote: ", outfile)
  }

  for (i in names(de_list)) {

    celltypeName <- gsub("/", "_", gsub(" ", "_", i))

    if (type == "DE") {
      res <- as.data.frame(de_list[[i]])
      if (is.null(res) || nrow(res) == 0) {
        message("Skipping ", i, " (no results)")
        next
      }
      outfile <- paste0(outSubDir, celltypeName, "_DE_Results.csv")
      utils::write.csv(res, file = outfile, row.names = FALSE)
      message("Wrote: ", outfile)
      next
    }

    if (split_direction) {
      dir_list <- de_list[[i]]
      if (is.null(dir_list)) {
        message("Skipping ", i, " (no results)")
        next
      }
      .write_go(dir_list$up,   celltypeName, "up")
      .write_go(dir_list$down, celltypeName, "down")
    } else {
      .write_go(de_list[[i]], celltypeName)
    }
  }

  invisible(NULL)
}

#' Run Gene Ontology Over-Representation Analysis on Multiple DEG Sets
#'
#' Performs GO enrichment analysis for each cluster-specific DEG list using
#' \code{clusterProfiler::enrichGO}. Produces a QC summary table of tested
#' versus retained gene sets per cluster.
#'
#' @param de_list A named list of data frames, each containing at least a
#'   column \code{gene} with gene symbols. Typically the output of
#'   \code{\link{run_cluster_de}}.
#' @param org_db An \code{OrgDb} annotation object (e.g.,
#'   \code{org.Mm.eg.db::org.Mm.eg.db} for mouse,
#'   \code{org.Hs.eg.db::org.Hs.eg.db} for human).
#' @param universe Character vector of background gene symbols representing all
#'   genes tested in the differential expression analysis. Defines the
#'   enrichment universe.
#' @param keyType Character string specifying the gene identifier type.
#'   Default \code{"SYMBOL"}.
#' @param qvalueCutoff Numeric. q-value (BH-adjusted p-value) threshold for
#'   retaining enriched terms. Default \code{0.1}.
#' @param ont Character string specifying the GO sub-ontology to test:
#'   \code{"BP"} (Biological Process), \code{"MF"} (Molecular Function),
#'   \code{"CC"} (Cellular Component), or \code{"all"}. Default \code{"bp"}.
#' @param outputDir Character path to the directory where the QC summary CSV
#'   is written. Created if it does not exist.
#' @param outname Character filename for the QC summary CSV. Default
#'   \code{"GO_ORA_QC_summary.csv"}.
#'
#' @return A named list of \code{enrichResult} objects (one per cluster).
#'   Elements are \code{NULL} for clusters with no DEGs. A QC summary CSV is
#'   also written to \code{<outputDir>/<outname>}.
#'
#' @details For each cluster, the function:
#' \enumerate{
#'   \item Extracts unique gene symbols from the DE data frame.
#'   \item Calls \code{clusterProfiler::enrichGO} with \code{minGSSize = 5},
#'     \code{maxGSSize = 300}, BH adjustment, and \code{readable = TRUE}.
#'   \item Annotates the result with the cluster name.
#'   \item Prints a per-cluster summary of tested/retained gene sets.
#' }
#'
#' @importFrom clusterProfiler enrichGO
#' @importFrom utils write.csv
#'
#' @export
run_go_ora <- function(de_list,
                       org_db,
                       universe = universe,
                       keyType = "SYMBOL",
                       qvalueCutoff = 0.1,
                       ont = "bp",
                       outputDir = NULL,
                       outname = "GO_ORA_QC_summary.csv") {

  ora_results <- lapply(names(de_list), function(clust) {

    degs <- de_list[[clust]]

    if (is.null(degs) || nrow(degs) == 0) {
      message("Skipping ", clust, " (no DEGs)")
      return(NULL)
    }

    genes <- unique(degs$gene)

    set.seed(1234)
    ego <- clusterProfiler::enrichGO(
      gene          = genes,
      OrgDb         = org_db,
      keyType       = keyType,
      universe      = universe,
      ont           = ont,
      minGSSize     = 5,
      maxGSSize     = 300,
      pAdjustMethod = "BH",
      qvalueCutoff  = qvalueCutoff,
      readable      = TRUE
    )

    if (!is.null(ego)) {
      ego@result$cluster <- clust
    }

    return(ego)
  })

  names(ora_results) <- names(de_list)

  qc_df <- data.frame(
    set = character(),
    tested = integer(),
    retained = integer(),
    prop = numeric(),
    stringsAsFactors = FALSE
  )

  for (i in names(ora_results)) {
    message("\n==============================")
    message("GO set: ", i)
    message("==============================")
    x <- ora_results[[i]]
    if (is.null(x) || nrow(x@result) == 0) {
      message("No results for: ", i)
      qc_df <- rbind(qc_df, data.frame(
        set = i, tested = 0, retained = 0, prop = NA
      ))
      next
    }

    total <- nrow(x@result)
    y <- as.data.frame(x)
    filt <- nrow(y)
    prop <- ifelse(total > 0, round(filt / total, 3), NA)
    message("Genesets tested      : ", total)
    message("After filtering      : ", filt)
    message("Retention proportion : ", prop)

    qc_df <- rbind(qc_df, data.frame(
      set = i, tested = total, retained = filt, prop = prop
    ))
  }

  if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE)
  utils::write.csv(qc_df, file = paste0(outputDir, outname))
  return(ora_results)
}

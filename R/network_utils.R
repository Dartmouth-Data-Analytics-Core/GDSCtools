#' Run GO-ORA Overlap Networks for All Clusters in a Comparison
#'
#' Iterates over every cluster in a named list of GO enrichment results,
#' builds a Jaccard-overlap network of enriched GO terms, applies Leiden
#' community detection, and saves network plots and summary tables. Clusters
#' that lack enrichment results or have fewer than two significant pathways
#' are skipped with a message.
#'
#' @param go_list A named list of \code{enrichResult} objects (from
#'   \code{clusterProfiler::enrichGO}), one per cell-type cluster. Names must
#'   match the names in \code{de_list}.
#' @param de_list A named list of differential-expression data frames, one per
#'   cluster. Each data frame must contain at least the columns \code{gene}
#'   (gene symbol) and \code{avg_log2FC}.
#' @param comp_label Character string identifying the comparison (e.g.,
#'   \code{"CD206_Car_vs_Bleo"}). Used to build output sub-directory paths and
#'   as a prefix in plot titles.
#' @param base_outdir Character path to the root output directory. A
#'   sub-directory \code{<comp_label>/networks/} is created automatically.
#' @param overlap_threshold Numeric \code{[0, 1]}. Minimum Jaccard overlap
#'   between two GO term gene sets required to draw an edge. Default
#'   \code{0.1}.
#' @param leiden_resolution Numeric. Resolution parameter passed to
#'   \code{igraph::cluster_leiden}. Lower values produce fewer, larger
#'   communities. Default \code{0.03}.
#'
#' @return A named list of network result lists (one per cluster that produced
#'   a valid network). Each element is the list returned by
#'   \code{\link{generate_overlap_network}}.
#'
#' @seealso \code{\link{build_network_inputs}},
#'   \code{\link{generate_overlap_network}}
#'
#' @export
run_comparison_networks <- function(
    go_list,
    de_list,
    comp_label,
    base_outdir,
    overlap_threshold = 0.1,
    leiden_resolution = 0.03
) {

  net_outdir <- paste0(base_outdir, comp_label, "/networks/")
  if (!dir.exists(net_outdir)) dir.create(net_outdir, recursive = TRUE)

  results <- lapply(names(go_list), function(clust) {

    clust_raw  <- clust
    clust_safe <- gsub("/", "_", clust_raw)

    ego   <- go_list[[clust_raw]]
    de_df <- de_list[[clust_raw]]

    if (is.null(ego) || is.null(de_df) || nrow(de_df) == 0) {
      message("Skipping network for ", clust_raw)
      return(NULL)
    }

    inputs <- build_network_inputs(ego, de_df)

    if (is.null(inputs) || length(inputs$collection) < 2) {
      message("Not enough pathways for: ", clust_raw)
      return(NULL)
    }

    clust_label <- paste0(comp_label, "_", gsub(" ", "_", clust_safe))

    generate_overlap_network(
      results_set       = inputs$results_set,
      collection        = inputs$collection,
      de_genes          = inputs$de_genes,
      collection_name   = clust_label,
      file_stem         = clust_safe,
      overlap_threshold = overlap_threshold,
      leiden_resolution = leiden_resolution,
      outdir            = net_outdir
    )
  })

  names(results) <- gsub("/", "_", names(go_list))
  results <- results[!sapply(results, is.null)]

  return(results)
}

#' Build Network Input Data from GO Enrichment and DE Results
#'
#' Extracts significant GO terms (adjusted p-value < 0.05) from an
#' \code{enrichResult} object and pairs them with fold-change information from
#' a differential-expression table. The output is a list containing:
#' (1) a named list of gene-set members (\code{collection}), (2) a data frame
#' of mean log2 fold-changes per pathway (\code{results_set}), and (3) the
#' unique set of DE genes (\code{de_genes}). These three components are the
#' required inputs for \code{\link{generate_overlap_network}}.
#'
#' @param ego An \code{enrichResult} object (or coercible via
#'   \code{as.data.frame}) produced by \code{clusterProfiler::enrichGO} or a
#'   similar ORA function. Must contain columns \code{Description},
#'   \code{p.adjust}, and \code{geneID} (slash-delimited gene symbols).
#' @param de_df A data frame of differentially expressed genes with at least
#'   columns \code{gene} and \code{avg_log2FC}.
#'
#' @return A list with three elements:
#' \describe{
#'   \item{collection}{Named list of character vectors; each element contains
#'     the gene symbols associated with one GO term.}
#'   \item{results_set}{Data frame with one row per GO term and a single
#'     column \code{avg_log2FC} giving the mean fold-change of DE genes in
#'     that term's gene set.}
#'   \item{de_genes}{Character vector of unique DE gene symbols.}
#' }
#' Returns \code{NULL} if no significant terms pass filtering.
#'
#' @seealso \code{\link{generate_overlap_network}}
#'
#' @importFrom stats setNames
#'
#' @export
build_network_inputs <- function(ego, de_df) {

  .extract <- function(ego_obj) {
    if (is.null(ego_obj) || nrow(as.data.frame(ego_obj)) == 0) return(NULL)
    res <- as.data.frame(ego_obj)
    res[!is.na(res$Description) & res$p.adjust < 0.05, ]
  }

  res <- .extract(ego)

  if (is.null(res) || nrow(res) == 0) return(NULL)

  res <- res[!duplicated(res$Description), ]

  collection <- stats::setNames(
    lapply(res$geneID, function(gids) unlist(strsplit(gids, "/"))),
    res$Description
  )

  de_genes <- unique(de_df$gene)
  fc_lookup <- stats::setNames(de_df$avg_log2FC, de_df$gene)

  pathway_fc <- sapply(collection, function(genes) {
    in_de <- intersect(genes, names(fc_lookup))
    if (length(in_de) == 0) return(0)
    mean(fc_lookup[in_de])
  })

  results_set <- data.frame(
    avg_log2FC = pathway_fc,
    row.names  = names(collection)
  )

  list(collection = collection, results_set = results_set, de_genes = de_genes)
}

#' Generate a Jaccard-Overlap Network of GO Terms with Leiden Clustering
#'
#' Constructs an undirected, weighted network where nodes are enriched GO terms
#' and edge weights are Jaccard overlap coefficients between their DE gene
#' sets. Edges below \code{overlap_threshold} are removed. Leiden community
#' detection is applied to identify clusters of functionally related terms.
#' The function produces two PNG plots (a plain up/down-regulation view and a
#' Leiden-clustered view) and exports CSV summaries of node-level and
#' cluster-level statistics.
#'
#' @param results_set Data frame with row names equal to GO term descriptions
#'   and at least one column \code{avg_log2FC} (mean fold-change across member
#'   genes).
#' @param collection Named list of character vectors. Each element is a set of
#'   gene symbols belonging to one GO term. Names must match
#'   \code{rownames(results_set)}.
#' @param de_genes Character vector of all differentially expressed gene
#'   symbols. Gene sets in \code{collection} are intersected with this vector
#'   before computing overlaps.
#' @param collection_name Character label used in plot titles and console
#'   output. Typically formatted as
#'   \code{"<comparison>_<celltype>"} (e.g.,
#'   \code{"CD206_Car_vs_Bleo_Fibroblast_1"}).
#' @param overlap_threshold Numeric \code{[0, 1]}. Jaccard coefficient below
#'   which edges are set to zero. Default \code{0.1}.
#' @param leiden_resolution Numeric. Resolution parameter for
#'   \code{igraph::cluster_leiden}. Default \code{0.03}.
#' @param outdir Character path to the output directory. If \code{NULL},
#'   no files are written. The directory is created if it does not exist.
#' @param plot_to_screen Logical. If \code{TRUE} (default), the Leiden-
#'   clustered network is also drawn to the active graphics device (e.g.,
#'   the RStudio Plots pane).
#' @param file_stem Character string used as the filename prefix for all
#'   exported files. If \code{NULL} or empty, \code{collection_name} is used.
#'
#' @return A list with six elements:
#' \describe{
#'   \item{network}{An \code{igraph} graph object with vertex attributes
#'     \code{fc}, \code{regulation}, \code{size}, \code{cluster},
#'     \code{color}, and \code{shape}.}
#'   \item{clustering}{The \code{communities} object returned by
#'     \code{igraph::cluster_leiden}.}
#'   \item{overlap_matrix}{Square numeric matrix of pairwise Jaccard
#'     coefficients (before thresholding).}
#'   \item{go_term_df}{Data frame with columns \code{go_term},
#'     \code{cluster}, \code{degree}, \code{regulation}, \code{avg_log2FC},
#'     and \code{leiden_resolution}.}
#'   \item{cluster_stats}{Data frame with columns \code{cluster},
#'     \code{avg_degree}, and \code{n_terms}.}
#'   \item{leiden_resolution}{The resolution value used.}
#' }
#'
#' @section Side effects:
#' When \code{outdir} is non-\code{NULL}, the following files are written:
#' \itemize{
#'   \item \code{<file_stem>_network_summary.csv} -- node/edge/cluster counts.
#'   \item \code{<file_stem>-cluster-stats.csv} -- per-cluster average degree.
#'   \item \code{<file_stem>-plain.png} -- network coloured by regulation
#'     direction (red = up, blue = down).
#'   \item \code{<file_stem>-clustered-leiden.png} -- network coloured by
#'     Leiden cluster membership.
#' }
#'
#' @seealso \code{\link{build_network_inputs}},
#'   \code{\link{write_cluster_results}}
#'
#' @importFrom igraph graph_from_adjacency_matrix simplify V E cluster_leiden
#'   membership degree delete_vertices vcount ecount layout_with_fr
#' @importFrom stats aggregate
#' @importFrom utils write.csv
#' @importFrom grDevices png dev.off
#' @importFrom graphics plot legend
#'
#' @export
generate_overlap_network <- function(
    results_set,
    collection,
    de_genes,
    collection_name = "network",
    overlap_threshold = 0.1,
    leiden_resolution = 0.03,
    outdir = NULL,
    plot_to_screen = TRUE,
    file_stem = NULL
) {

  rownames(results_set) <- gsub("-", "_", rownames(results_set))
  de_genes <- unique(gsub("-", "_", de_genes))
  names(collection) <- gsub("-", "_", names(collection))

  collection_sub <- collection[names(collection) %in% rownames(results_set)]

  collection_sub <- lapply(collection_sub, function(g) {
    intersect(gsub("-", "_", g), de_genes)
  })
  collection_sub <- collection_sub[sapply(collection_sub, length) > 0]

  if (length(collection_sub) < 2) {
    warning(paste("Not enough pathways for:", collection_name))
  }

  nms <- names(collection_sub)
  mat <- matrix(0, nrow = length(nms), ncol = length(nms),
                dimnames = list(nms, nms))

  for (i in seq_along(nms)) {
    for (j in seq_along(nms)) {
      gi <- unique(collection_sub[[i]])
      gj <- unique(collection_sub[[j]])
      overlap    <- length(intersect(gi, gj))
      union_size <- length(unique(c(gi, gj)))
      mat[i, j]  <- ifelse(union_size > 0, overlap / union_size, 0)
    }
  }

  mat_thres <- mat
  mat_thres[mat_thres < overlap_threshold] <- 0

  network <- igraph::graph_from_adjacency_matrix(
    mat_thres, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  network <- igraph::simplify(network, remove.loops = TRUE, remove.multiple = TRUE)

  fc <- results_set$avg_log2FC[match(igraph::V(network)$name, rownames(results_set))]
  igraph::V(network)$fc <- fc
  igraph::V(network)$regulation <- ifelse(fc > 0, "up", "down")
  igraph::V(network)$size <- log(sapply(collection_sub[igraph::V(network)$name], length)) * 2

  isolates <- igraph::degree(network) == 0
  network <- igraph::delete_vertices(network, isolates)

  if (igraph::vcount(network) < 2) {
    warning(paste("Network too sparse after filtering:", collection_name))
  }

  igraph::E(network)$width <- igraph::E(network)$weight

  ceb <- igraph::cluster_leiden(network, weights = igraph::E(network)$weight,
                                resolution_parameter = leiden_resolution)
  mem <- igraph::membership(ceb)
  igraph::V(network)$cluster <- as.integer(mem)

  distinct_colors <- c(
    "#E31A1C", "#1F78B4", "#33A02C", "#FF7F00",
    "#6A3D9A", "#E6AB02", "#66C2A5", "#E7298A",
    "#1B9E77", "#A6761D", "#666666", "#D95F02"
  )
  n_clusters <- length(unique(mem))
  cluster_palette <- rep_len(distinct_colors, n_clusters)
  igraph::V(network)$color <- cluster_palette[mem]
  igraph::V(network)$shape <- ifelse(igraph::V(network)$regulation == "up", "circle", "square")

  node_degrees <- igraph::degree(network)
  go_term_df <- data.frame(
    go_term    = igraph::V(network)$name,
    cluster    = igraph::V(network)$cluster,
    degree     = node_degrees,
    regulation = igraph::V(network)$regulation,
    avg_log2FC = igraph::V(network)$fc,
    leiden_resolution = leiden_resolution,
    stringsAsFactors = FALSE
  )

  cluster_stats <- stats::aggregate(degree ~ cluster, data = go_term_df, FUN = mean)
  colnames(cluster_stats) <- c("cluster", "avg_degree")
  cluster_stats$n_terms <- as.integer(table(go_term_df$cluster))

  message("\n===== Network Summary: ", collection_name, " =====")
  message("  Nodes:              ", igraph::vcount(network))
  message("  Edges:              ", igraph::ecount(network))
  message("  Clusters:           ", n_clusters)
  message("  Leiden resolution:  ", leiden_resolution)
  message("  Overlap threshold:  ", overlap_threshold)
  message("\n  Average degree per cluster:")
  for (i in seq_len(nrow(cluster_stats))) {
    message("    Cluster ", cluster_stats$cluster[i],
            ": avg_degree = ", round(cluster_stats$avg_degree[i], 2),
            " (", cluster_stats$n_terms[i], " terms)")
  }
  message("==========================================\n")

  if (!is.null(outdir)) {
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

    if (is.null(file_stem) || !nzchar(file_stem)) {
      file_stem <- collection_name
    }

    summary_df <- data.frame(
      metric = c("nodes", "edges", "clusters"),
      count  = c(igraph::vcount(network), igraph::ecount(network), n_clusters),
      stringsAsFactors = FALSE
    )
    summary_file <- paste0(outdir, file_stem, "_network_summary.csv")
    utils::write.csv(summary_df, file = summary_file, row.names = FALSE)
    message("Saved network summary CSV: ",
            normalizePath(summary_file, winslash = "/", mustWork = FALSE),
            " (exists=", file.exists(summary_file), ")")

    cluster_stats_out <- cluster_stats[, c("cluster", "avg_degree", "n_terms")]
    cluster_stats_file <- paste0(outdir, file_stem, "-cluster-stats.csv")
    utils::write.csv(cluster_stats_out, file = cluster_stats_file, row.names = FALSE)
    message("Saved cluster stats CSV: ",
            normalizePath(cluster_stats_file, winslash = "/", mustWork = FALSE),
            " (exists=", file.exists(cluster_stats_file), ")")

    set.seed(1234)
    lo <- igraph::layout_with_fr(network)
    ppi <- 300

    parts  <- strsplit(collection_name, "_vs_")[[1]]
    ident1 <- if (length(parts) >= 1) gsub("_", " ", parts[1]) else "ident1"
    ident2 <- if (length(parts) >= 2) gsub("_", " ", parts[2]) else "ident2"

    up_down_cols <- ifelse(igraph::V(network)$regulation == "up", "#E67A6E", "#278EEA")

    out_png_plain <- paste0(outdir, file_stem, "-plain.png")
    grDevices::png(filename = out_png_plain,
                   width = 10 * ppi, height = 10 * ppi, res = ppi)
    graphics::plot(
      network,
      layout             = lo,
      vertex.color       = up_down_cols,
      vertex.label       = NA,
      vertex.frame.color = "#555555",
      main               = paste0(collection_name, " (pathway overlap)")
    )
    graphics::legend("bottomleft",
                     legend = c(paste0("Up in ", ident1), paste0("Down in ", ident2)),
                     pch = c(16, 15), col = "black", pt.cex = 2, bty = "n", cex = 0.7)
    grDevices::dev.off()

    if (identical(Sys.info()[["sysname"]], "Darwin")) {
      try(system2("chflags", c("nohidden", out_png_plain)), silent = TRUE)
    }
    message("Saved network PNG: ", normalizePath(out_png_plain, winslash = "/", mustWork = TRUE),
            " (exists=", file.exists(out_png_plain), ")")

    tryCatch({
      cluster_ids <- sort(unique(mem))
      out_png <- paste0(outdir, file_stem, "-clustered-leiden.png")

      grDevices::png(filename = out_png,
                     width = 12.5 * ppi, height = 10 * ppi, res = ppi)
      graphics::plot(
        network,
        layout             = lo,
        vertex.label       = NA,
        vertex.frame.color = "#333333",
        main               = paste0(collection_name, " (Leiden)")
      )
      graphics::legend("bottomleft",
                       legend = c(paste0("Up (", ident1, ")"), paste0("Down (", ident2, ")")),
                       pch = c(16, 15), col = "black", pt.cex = 2, bty = "n", cex = 0.7)
      graphics::legend("bottomright",
                       legend = paste("Cluster", cluster_ids),
                       pch = 21,
                       pt.bg = cluster_palette[match(cluster_ids, sort(unique(cluster_ids)))],
                       col = "black", pt.cex = 2, bty = "n", cex = 0.7)
      grDevices::dev.off()

      if (identical(Sys.info()[["sysname"]], "Darwin")) {
        try(system2("chflags", c("nohidden", out_png)), silent = TRUE)
      }
      message("Saved network PNG: ", normalizePath(out_png, winslash = "/", mustWork = FALSE),
              " (exists=", file.exists(out_png), ")")
    }, error = function(e) {
      try(grDevices::dev.off(), silent = TRUE)
      message("FAILED to save Leiden PNG for ", collection_name, ": ", conditionMessage(e))
    })
  }

  if (isTRUE(plot_to_screen)) {
    set.seed(1234)
    lo <- igraph::layout_with_fr(network)
    parts  <- strsplit(collection_name, "_vs_")[[1]]
    ident1 <- if (length(parts) >= 1) gsub("_", " ", parts[1]) else "ident1"
    ident2 <- if (length(parts) >= 2) gsub("_", " ", parts[2]) else "ident2"
    cluster_ids <- sort(unique(mem))
    graphics::plot(
      network,
      layout             = lo,
      vertex.label       = NA,
      vertex.frame.color = "#333333",
      main               = paste0(collection_name, " (Leiden)")
    )
    graphics::legend("bottomleft",
                     legend = c(paste0("Up (", ident1, ")"), paste0("Down (", ident2, ")")),
                     pch = c(16, 15), col = "black", pt.cex = 1.5, bty = "n", cex = 0.7)
    graphics::legend("bottomright",
                     legend = paste("Cluster", cluster_ids),
                     pch = 21, pt.bg = cluster_palette[cluster_ids],
                     col = "black", pt.cex = 1.5, bty = "n", cex = 0.7)
  }

  list(
    network          = network,
    clustering       = ceb,
    overlap_matrix   = mat,
    go_term_df       = go_term_df,
    cluster_stats    = cluster_stats,
    leiden_resolution = leiden_resolution
  )
}

#' Extract and Export GO Terms from a Specific Leiden Cluster
#'
#' Given a network result object and a single GO term of interest, identifies
#' which Leiden cluster that term belongs to and returns (and optionally writes)
#' a data frame of all GO terms sharing that cluster. This is useful for
#' inspecting the biological theme captured by a particular network community.
#'
#' @param net_result A list returned by \code{\link{generate_overlap_network}},
#'   which must contain at least the element \code{go_term_df}.
#' @param go_term Character string. The GO term description to look up (e.g.,
#'   \code{"extracellular matrix organization"}). Hyphens are internally
#'   replaced with underscores to match network node naming conventions.
#' @param outdir Character path to the output directory. If non-\code{NULL},
#'   a CSV file is written. The directory is created if it does not exist.
#' @param label Optional character label used in the output filename. If
#'   \code{NULL}, defaults to \code{"cluster-<id>"}.
#'
#' @return A data frame containing all GO terms in the same Leiden cluster as
#'   \code{go_term}, with columns \code{go_term}, \code{cluster},
#'   \code{degree}, \code{regulation}, \code{avg_log2FC}, and
#'   \code{leiden_resolution}. Returns \code{NULL} with a warning if the
#'   requested term is not found in the network.
#'
#' @seealso \code{\link{generate_overlap_network}}
#'
#' @importFrom utils write.csv
#'
#' @export
write_cluster_results <- function(net_result, go_term, outdir, label = NULL) {

  go_term_df <- net_result$go_term_df
  target_cluster <- go_term_df$cluster[go_term_df$go_term == gsub("-", "_", go_term)]

  if (length(target_cluster) == 0) {
    warning("GO term '", go_term, "' not found in network")
    return(NULL)
  }

  cluster_terms <- go_term_df[go_term_df$cluster == target_cluster[1], ]

  if (!is.null(outdir)) {
    if (is.null(label)) label <- paste0("cluster-", target_cluster[1])
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
    utils::write.csv(cluster_terms, paste0(outdir, "network-cluster-", label, ".csv"),
                     row.names = FALSE)
  }

  cluster_terms
}

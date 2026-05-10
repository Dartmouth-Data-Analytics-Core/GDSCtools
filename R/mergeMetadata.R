#' @title mergeMetadata
#'
#' @description
#' Efficiently merge metadata across multiple Seurat objects into one master dataframe.
#'
#' @param x A named-list of Seurat objects
#'
#' @returns A dataframe of merged metadata
#'
#' @examples # list of seurat objects with names
#' objects <- c("obj1", "obj2", "obj3", "obj4")
#' # Vector of sample IDs
#' sampleIDs <- c("S1", "S2", "S3", "S4")
#' names(objects) <- sampleIDs
#' metadata <- mergeMetadata(objects)
#'
#' @import Seurat
#'
#' @export

mergeMetadata <- function(x) {

  #----- Check that sample names match length of list
  if (!is.list(x)) {
    stop("Input must be a list of Seurat objects.")
  }

  #----- Assign sample names if not provided
  sample_names <- names(x)
  if (is.null(sample_names) || any(sample_names == "")) {
    sample_names <- paste0("sample", seq_along(x))
  }

  #----- Merge metadata
  merged_metadata <- do.call(rbind, lapply(seq_along(x), function(i) {
    meta <- x[[i]]@meta.data
    meta$orig.ident <- sample_names[i]
    meta$cell_id <- rownames(meta)
    return(meta)
  }))

  #----- Return metadata
  return(merged_metadata)
}

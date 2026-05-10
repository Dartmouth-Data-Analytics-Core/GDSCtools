test_that("write_de_results writes DE CSV files", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_de_out/")

  de_list <- list(
    ClusterA = data.frame(
      gene = c("G1", "G2"),
      avg_log2FC = c(1.0, -0.5),
      p_val_adj = c(0.01, 0.03),
      stringsAsFactors = FALSE
    ),
    ClusterB = data.frame(
      gene = c("G3"),
      avg_log2FC = c(2.0),
      p_val_adj = c(0.001),
      stringsAsFactors = FALSE
    )
  )

  write_de_results(de_list, outdir, type = "DE")

  expect_true(file.exists(paste0(outdir, "ClusterA_DE_Results.csv")))
  expect_true(file.exists(paste0(outdir, "ClusterB_DE_Results.csv")))

  unlink(outdir, recursive = TRUE)
})

test_that("write_de_results skips empty DE clusters", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_de_empty/")

  de_list <- list(
    ClusterA = data.frame(
      gene = character(0),
      avg_log2FC = numeric(0),
      stringsAsFactors = FALSE
    )
  )

  expect_message(
    write_de_results(de_list, outdir, type = "DE"),
    "Skipping"
  )

  expect_false(file.exists(paste0(outdir, "ClusterA_DE_Results.csv")))

  unlink(outdir, recursive = TRUE)
})

test_that("write_de_results requires split_direction for GO type", {
  de_list <- list(ClusterA = data.frame(gene = "G1"))

  expect_error(
    write_de_results(de_list, tempdir(), type = "GO"),
    "split_direction"
  )
})

test_that("write_de_results handles GO without split", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_go_nosplit/")

  mock_ego <- data.frame(
    Description = "immune response",
    GeneRatio = "5/100",
    p.adjust = 0.01,
    geneID = "G1/G2/G3/G4/G5",
    stringsAsFactors = FALSE
  )

  de_list <- list(ClusterA = mock_ego)

  write_de_results(de_list, outdir, type = "GO", split_direction = FALSE)

  expect_true(file.exists(paste0(outdir, "ClusterA_GO_ORA_Results.csv")))

  written <- read.csv(paste0(outdir, "ClusterA_GO_ORA_Results.csv"))
  expect_true("geneRatio_num" %in% colnames(written))
  expect_true("geneRatio_den" %in% colnames(written))
  expect_equal(written$geneRatio_num, 5)
  expect_equal(written$geneRatio_den, 100)

  unlink(outdir, recursive = TRUE)
})

test_that("write_de_results handles GO with split direction", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_go_split/")

  mock_up <- data.frame(
    Description = "pathway_up",
    GeneRatio = "3/50",
    p.adjust = 0.02,
    geneID = "G1/G2/G3",
    stringsAsFactors = FALSE
  )

  mock_down <- data.frame(
    Description = "pathway_down",
    GeneRatio = "2/50",
    p.adjust = 0.04,
    geneID = "G4/G5",
    stringsAsFactors = FALSE
  )

  de_list <- list(ClusterA = list(up = mock_up, down = mock_down))

  write_de_results(de_list, outdir, type = "GO", split_direction = TRUE)

  expect_true(file.exists(paste0(outdir, "ClusterA_GO_ORA_up_Results.csv")))
  expect_true(file.exists(paste0(outdir, "ClusterA_GO_ORA_down_Results.csv")))

  unlink(outdir, recursive = TRUE)
})

test_that("write_de_results sanitizes cluster names with slashes and spaces", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_de_sanitize/")

  de_list <- list(
    `Cell Type/Sub 1` = data.frame(
      gene = "G1",
      avg_log2FC = 1.0,
      stringsAsFactors = FALSE
    )
  )

  write_de_results(de_list, outdir, type = "DE")

  expect_true(file.exists(paste0(outdir, "Cell_Type_Sub_1_DE_Results.csv")))

  unlink(outdir, recursive = TRUE)
})

test_that("write_de_results returns invisible NULL", {
  tmp <- tempdir()
  outdir <- paste0(tmp, "/test_de_return/")

  de_list <- list(
    A = data.frame(gene = "G1", avg_log2FC = 1.0, stringsAsFactors = FALSE)
  )

  result <- write_de_results(de_list, outdir, type = "DE")

  expect_null(result)

  unlink(outdir, recursive = TRUE)
})

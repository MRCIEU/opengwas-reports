suppressPackageStartupMessages({
  library("testthat")
  library("tidyverse")
  library("here")
  library("glue")
})

source(here("funcs/utils.R"))

context("unlog")

test_that("unlog, no unlog", {
  val <- 1000
  neg_log_10_val <- -log10(val)
  unlog_val <- unlog(val)
  expect_equal(val, unlog_val)
})

test_that("unlog, -log10", {
  val <- 1000
  neg_log10_val <- -log10(val)
  unlog_val <- unlog(neg_log10_val, TRUE)
  expect_equal(val, unlog_val)
})

test_that("unlog, -log2", {
  val <- 1000
  neg_log2_val <- -log2(val)
  unlog_val <- unlog(neg_log2_val, TRUE, function(x) 2^(-x))
  expect_equal(val, unlog_val)
})

context("split by chunk")

test_that("split_by_chunk", {
  vec <- 1:10
  n_chunks <- 3
  chunk_idx <- 3
  expected_results <- list(1:4, 5:8, 9:10)
  results <- list(
    split_by_chunk(vec, n_chunks, 1),
    split_by_chunk(vec, n_chunks, 2),
    split_by_chunk(vec, n_chunks, 3)
  )
  expect_equal(expected_results, results)
})

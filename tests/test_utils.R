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

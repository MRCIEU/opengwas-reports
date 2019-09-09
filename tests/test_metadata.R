suppressPackageStartupMessages({
  library("testthat")
  library("tidyverse")
  library("here")
  library("glue")
})

source(here("funcs/metadata.R"))

context("Parse metadata")

test_that("Parse harmonised variants", {
  field <- "<ID=IEU-a:1,TotalVariants=2610590,VariantsNotRead=0,HarmonisedVariants=2610590,VariantsNotHarmonised=0,SwitchedAlleles=2610584,StudyType=Continuous>"
  expected_res <- 2610590
  res <- parse_harmonised_variants(field)
  expect_equal(expected_res, res)
})

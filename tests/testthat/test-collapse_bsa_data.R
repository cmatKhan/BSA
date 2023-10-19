test_that("collapse_bsa_data works correctly", {

  df_with_meta <- data.frame(
    A = rep(c("group1", "group2"), each = 5),
    CHR = rep("CHR1", 10),
    POS = rep(c(1,5,10,20,30), times=2),
    REF_Allele = rep('A', 10),
    ALT1 = rep('G', 10),
    RealDepth = runif(10, 1, 10),
    Depth = runif(10, 1, 10),
    Reference = runif(10, 1, 10),
    Alternative1 = runif(10, 1, 10)
  )
  grouping_conditions <- c('CHR','POS','REF_Allele','ALT1')
  sample_columns <- c("RealDepth", "Depth", "Reference", "Alternative1")

  result <- collapse_bsa_data(df_with_meta, grouping_conditions, sample_columns)

  # Add your assertions here. For instance, check if the result has the correct number of rows:
  expect_equal(nrow(result), 5)
  expected_cols <- c("CHR", "POS", "REF_Allele", "ALT1", 'sample', 'RealDepth', 'Depth', 'Reference', 'Alternative1')
  expect_equal(colnames(result), expected_cols)
})

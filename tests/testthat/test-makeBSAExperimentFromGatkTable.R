# setup fixtures/temporary files, etc ----
setup({

  # Create the GATK table in code
  gatk_table <- data.frame(
    CHROM = rep("CP022321.1", 5),
    POS = c(222, 558, 1443, 10528, 12097),
    REF = c("G", "C", "C", "A", "T"),
    ALT = c("A", "G", "T", "T", "A"),
    `MULTI-ALLELIC` = rep(FALSE, 5),
    TYPE = rep("SNP", 5),
    C8.AD = c("1,149", "3,159", "1,147", "10,95", "0,23"),
    C8.DP = c(150, 162, 148, 105, 23),
    C8.PL = c("4205,0", "4977,0", "4533,0", "1778,0", "490,0"),
    C8.GQ = c(99, 99, 99, 99, 99),
    C8.GT = c("A", "G", "T", "T", "A"),
    KN99a.AD = c("240,3", "237,0", "223,0", "42,0", "28,0"),
    KN99a.DP = c(243, 237, 223, 42, 28),
    KN99a.PL = c("0,7042", "0,7202", "0,6743", "0,1176", "0,700"),
    KN99a.GQ = c(99, 10, 6, 99, 99),
    KN99a.GT = c("G", "C", "C", "A", "T"),
    SLB0021.AD = c("225,65", "243,70", "230,61", "40,19", "31,8"),
    SLB0021.DP = c(290, 313, 292, 59, 39),
    SLB0021.PL = c("0,6352", "0,6879", "0,6604", "0,874", "0,620"),
    SLB0021.GQ = c(99, 10, 6, 99, 99),
    SLB0021.GT = c("G", "C", "C", "A", "T"),
    SLB0025.AD = c("84,17", "153,35", "143,24", "27,10", "54,4"),
    SLB0025.DP = c(101, 188, 167, 37, 58),
    SLB0025.PL = c("0,2666", "0,4700", "0,4721", "0,677", "0,1396"),
    SLB0025.GQ = c(99, 10, 6, 99, 99),
    SLB0025.GT = c("G", "C", "C", "A", "T"),
    check.names = FALSE
  )
  # Save the GATK table to a temporary file
  gatk_table_path <- tempfile()
  write.table(gatk_table, gatk_table_path, sep = "\t", row.names = FALSE)

  meta_df = data.frame(
    sample = c("SLB0021", "SLB0025"),
    group = c("ir7all", "ir7all"),
    pool = c(1, 1),
    replicate = c(1, 1),
    condition = c("inoculum", "Lung"),
    sac_day = c(0, 9),
    culture_time = c(0, 0)
  )

  meta_path <- tempfile(fileext='.csv')
  write.table(meta_df, meta_path, sep = ",", row.names = FALSE)

})

test_that("test .validate_makeBSAExperimentFromGatkTable valid input", {
  stub(.validate_makeBSAExperimentFromGatkTableInput, "file.exists", TRUE)

  expect_silent(.validate_makeBSAExperimentFromGatkTableInput(
    gatk_table_path = "path/to/nonexistent_gatk_table.txt",
    col_data_path = "path/to/nonexistent_col_data.csv",
    drop_samples = NULL,
    high_confidence_depth = 10,
    high_confidence_alt_percentage = 0.9,
    keep_multiallelic = FALSE,
    high_confidence_pl = NULL,
    high_confidence_gq = NULL
  ))
})

test_that("test .validate_makeBSAExperimentFromGatkTable error conditions", {
  stub(.validate_makeBSAExperimentFromGatkTableInput, "file.exists", FALSE)

  expect_error(
    .validate_makeBSAExperimentFromGatkTableInput(
      gatk_table_path = "nonexistent_gatk_table.txt", # Invalid because it does not exist
      col_data_path = "nonexistent_col_data.txt", # Invalid because it does not exist and wrong extension
      drop_samples = 123, # Invalid because it's not a character vector or NULL
      high_confidence_depth = -10, # Invalid because it's negative
      high_confidence_alt_percentage = 1.5, # Invalid because it's outside the 0 to 1 range
      keep_multiallelic = "yes", # Invalid because it's not a boolean
      high_confidence_pl = -5, # Invalid because it's negative
      high_confidence_gq = "high" # Invalid because it's not numeric
    ),
    # Check for a comprehensive error message containing all validation issues
    paste(
      "`keep_multiallelic` must be a boolean",
      "`drop_samples` must be a character vector or NULL",
      "`keep_multiallelic` is not yet implemented. Please set to FALSE.",
      "File nonexistent_gatk_table.txt does not exist",
      "`col_data_path` nonexistent_col_data.txt does not exist",
      "`col_data_path` must be a csv file with `.csv` as an extension. Verify that nonexistent_col_data.txt is a csv and change the extension.",
      "`high_confidence_depth` must be a positive integer",
      "`high_confidence_alt_percentage` must be a number between 0 and 1",
      "`high_confidence_pl` must be a positive number or NULL",
      "`high_confidence_gq` must be a positive number or NULL",
      sep = "\n"
    )
  )
})


test_that("test .read_in_gatk_table function", {
  # Call the function
  result <- .read_in_gatk_table(gatk_table_path)
  # Check that the result is a list with two elements: 'table' and 'samples'
  expect_type(result, "list")
  expect_equal(names(result), c("table", "samples"))

  # Check that 'table' is a data frame with the correct columns
  expect_s3_class(result$table, "data.frame")
  expect_equal(colnames(result$table), c("CHROM", "POS", "REF", "ALT", "MULTI-ALLELIC", "TYPE", "SLB0021.AD", "SLB0021.DP", "SLB0021.PL", "SLB0021.GQ", "SLB0021.GT", "KN99a.AD", "KN99a.DP", "KN99a.PL", "KN99a.GQ", "KN99a.GT"))

  # Check that 'samples' is a character vector with the correct values
  expect_type(result$samples, "character")
  expect_equal(result$samples, c("C8", "KN99a"))

  # Test the function with the 'drop_samples' parameter
  result_with_drop <- .read_in_gatk_table(gatk_table_path, drop_samples = c("C8"))
  expect_equal(result_with_drop$samples, c("KN99a", "SLB0021", "SLB0025"))
  expect_false("C8.AD" %in% colnames(result_with_drop$table))
})

test_that("parse_gt_table returns correct output", {
  # Define the input data frame and reference vector
  gt_data <- data.frame(sample1 = c("A", "G", "C", "T", "A"),
                        sample2 = c("G", "A", "C", "T", "G"))
  ref <- c("A", "G", "C", "T", "A")

  # Call the function
  result <- .parse_gt_table(gt_data, ref)

  # Check that the result is a matrix with the correct dimensions
  expect_is(result, "matrix")
  expect_equal(dim(result), c(length(ref), ncol(gt_data)))

  # Check that the result has the correct values
  expected_result <- matrix(c(0, 0, 0, 0, 0, 1, 1, 0, 0, 1), nrow = length(ref))
  colnames(expected_result) = colnames(gt_data)
  expect_equal(result, expected_result)
})

test_that("parse_ad_table returns correct output", {

  ad_data = matrix(c("1,149","3,159","1,147","10,95","0,23",
                     "240,3", "237,0","223,0","42,0","28,0"),
                   nrow=5)
  colnames(ad_data) = c('sample1', 'sample2')

  result = .parse_ad_table(ad_data)
  expect_is(result, 'matrix')
  expect_equal(dim(result), c(5, 2))

  expected_result =   ad_data = matrix(c(149,159,147,95,23,
                                         3,0,0,0,0),
                                       nrow=5)
  colnames(expected_result) = c('sample1', 'sample2')

  expect_equal(result, expected_result)
})

test_that("test makeBSAExperimentFromGatkTable constructs an appropriate BSAExperiment object", {
  # Mock or specify paths to your representative GATK table and col_data CSV files
  gatk_table_path <- system.file("extdata", "gatk_table_example.txt", package = "yourPackage")
  col_data_path <- system.file("extdata", "col_data_example.csv", package = "yourPackage")

  # Expected setup
  expected_samples <- c("SLB0021", "SLB0025") # Adjust based on your mock/real data
  expected_chromosomes <- c("chr1", "chr1") # Example values
  expected_positions <- c(100, 200) # Example positions

  # Execute function
  result <- makeBSAExperimentFromGatkTable(
    gatk_table_path = gatk_table_path,
    col_data_path = col_data_path,
    drop_samples = c('KN99a', 'C8'),
    high_confidence_depth = 10,
    high_confidence_alt_percentage = 0.9,
    keep_multiallelic = FALSE,
    high_confidence_pl = NULL,
    high_confidence_gq = NULL
  )

  # Test if the result is a BSAExperiment object
  expect_true(inherits(result, "BSAExperiment"))

  # Validate rowRanges
  row_ranges <- rowRanges(result)
  expect_equal(length(row_ranges), length(expected_samples)) # Adjust as needed
  expect_equal(seqnames(row_ranges), expected_chromosomes)
  expect_equal(start(row_ranges), expected_positions)

  # Validate samples in colData
  col_data <- colData(result)
  expect_equal(col_data$sample, expected_samples)

})

# tear down any fixtures/temporary files ----
teardown({
  unlink(gatk_table_path)
  unlink(meta_path)
})

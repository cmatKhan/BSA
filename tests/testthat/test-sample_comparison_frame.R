testthat::test_that("pool constructions frame works", {

  # NOTE! "P2.1I" is a singleton
  sample_example = c("P1.1I","P1.1L","P1.1Y",
                     "P1.2B","P1.2I","P1.2L","P1.2Y","P2.1I")
  pool_construction = tibble::tibble(sample = sample_example) %>%
    dplyr::mutate(batch = stringr::str_remove(sample, "\\w$")) %>%
    dplyr::mutate(cond = ifelse(stringr::str_detect(sample, "P[[:alnum:]].{1,3}I"),'inoculum', NA)) %>%
    dplyr::mutate(cond = ifelse(stringr::str_detect(sample, "P[[:alnum:]].{1,3}Y"),'ypd', cond)) %>%
    dplyr::mutate(cond = ifelse(stringr::str_detect(sample, "P[[:alnum:]].{1,3}L"),'lung', cond)) %>%
    dplyr::mutate(cond = ifelse(stringr::str_detect(sample, "P[[:alnum:]].{1,3}B"),'brain', cond)) %>%
    dplyr::mutate(bulk = ifelse(cond == "inoculum", 'low', 'high'))

  # NOTE! "P2.1I" is excluded b/c it does not have any comparisons
  expected_frame = dplyr::tibble(
    lowBulk = c(rep("P1.1I",2), rep("P1.2I", 3)),
    highBulk = c("P1.1L","P1.1Y","P1.2B","P1.2L","P1.2Y")
  )

  all.equal(BSA::sample_comparison_frame(pool_construction), expected_frame)

})

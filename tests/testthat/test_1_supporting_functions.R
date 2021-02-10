context("tests of supporting functions")

test_that("load the demo data", {
  result <- data(PolyHaplotyper_demo)
  expect_equal(dim(snpdos), c(30, 663))
  expect_equal(dim(ped), c(661, 4))
  expect_equal(names(ped), c("genotype", "mother", "father", "sample_nr"))
})

test_that("test several functions", {
  result <- allHaplotypes(mrknames=letters[1:4])
  expect_equal(dim(result), c(16, 4))
  expect_equal(colnames(result), letters[1:4])
  expect_true(all(result[nrow(result),] == 1))
  expect_equal(result[1:2, ncol(result)], c(0, 1))
  result <- mrkdos2mrkdid(
    mrkDosage=matrix(c(1,2,3,4,4,3,2,1), ncol=2,
                     dimnames=list(letters[1:4], c("ind1", "ind2"))),
    ploidy=4)
  expect_equal(result, setNames(c(195, 587), c("ind1", "ind2")))
  parhac <- matrix(c(1,1,2,2,5,6,1,1,1,1,1,3), ncol=2,
                   dimnames=list(NULL, c("P1", "P2")))
  result <- getFSfreqs(parhac=parhac, DRrate=0.1)
  expect_equal(names(result), c("FShac", "freq"))
  expect_equal(dim(result$FShac), c(6, 54))
  expect_equal(sum(result$freq), 1, tolerance=1e-6)
})


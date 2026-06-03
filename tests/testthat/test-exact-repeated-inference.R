test_that("source random stream follows Delphi LCG random", {
  random_draw <- screen_j_source_random_stream(9L)

  expect_equal(
    as.numeric(replicate(5L, random_draw())),
    c(
      0.282419453840702772,
      0.498396688839420676,
      0.934840948320925236,
      0.835643683793023229,
      0.861519629601389170
    ),
    tolerance = 1e-15
  )
})

test_that("the report gives the age-class shares of each group, NA without area", {
  # grid cell 2 has no area, grid cell 3 only rounding noise (below 1e-9 Mha)
  ages <- new.magpie(c("c1", "c2", "c3"), "y2000", c("ac0", "ac5"), fill = c(1, 0, 1e-22, 3, 0, 2e-22))
  localMagpieSources(list(), outputs = list("AgeClassesHighRes.input=magpie" =
                                              list(secdforest = ages, diagnostics = data.frame(year = "y2000"))))
  out <- calcAgeClassesReport("downscaledmagpie", "magpie")[["x"]]
  expect_equal(as.vector(out[["secdforest"]]["c1", , ]), c(0.25, 0.75))
  expect_true(all(is.na(out[["secdforest"]][c("c2", "c3"), , ])))
  expect_equal(out[["diagnostics"]], data.frame(year = "y2000"))
  expect_error(calcAgeClassesReport("ESM", "magpie"), "downscaledmagpie")
})

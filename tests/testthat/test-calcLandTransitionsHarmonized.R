test_that("without harmonisation the input passes unchanged, other methods stop", {
  input <- list(statesLr = array(1, c(1, 1, 1)), rules = list())
  localMagpieSources(list(), outputs = list("LandTransitionsInput.input=magpie" = input))
  expect_identical(calcLandTransitionsHarmonized("magpie", "landuseinit", NULL, "none")[["x"]], input)
  expect_error(calcLandTransitionsHarmonized("magpie", "landuseinit", c(2015, 2050), "fade"), "only \"none\"")
})

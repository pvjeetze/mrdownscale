test_that("fullDOWNSCALEDMAGPIE checks the products before any calculation", {
  expect_error(fullDOWNSCALEDMAGPIE(product = c("land", "peatland")), "Unsupported product \"peatland\"")
  expect_error(fullDOWNSCALEDMAGPIE(downscaling = "proportional", harmonization = "fade",
                                    product = c("land", "BII")),
               "gives only the product \"land\"")
})

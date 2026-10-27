# helpers ------------------------------------------------------------------------------------

# magclass object without year of the grid cells c1, c2, ... (values grid cell by grid cell, item after item)
cellValues <- function(values, nCell, names) {
  return(new.magpie(paste0("c", seq_len(nCell)), NULL, names, fill = values))
}

# magclass object without year of units (values unit by unit, item after item)
unitValues <- function(values, units, names) {
  return(new.magpie(units, NULL, names, fill = values))
}

# transitions between three land types in two clusters with source and destination caps
makeCappedFlows <- function(seed = 2) {
  set.seed(seed)
  types <- c("a", "b", "c")
  flows <- expand.grid(from = types, to = types, stringsAsFactors = FALSE)
  flows <- flows[flows[["from"]] != flows[["to"]], ]
  flowNames <- paste(flows[["from"]], flows[["to"]], sep = ">")
  nCell <- 12
  cluster <- rep(c("J1", "J2"), each = 6)
  state <- cellValues(runif(nCell * 3), nCell, types)
  # unit totals up to 15 % of the unit area of the source type
  stateLr <- toolAggregate(state, cellMapping(cluster), from = "cell", to = "cluster")
  total <- unitValues(0, c("J1", "J2"), flowNames)
  for (j in c("J1", "J2")) {
    total[j, , ] <- runif(length(flowNames)) * 0.15 * as.vector(stateLr[j, , ])[match(flows[["from"]], types)]
  }
  prior <- cellValues(runif(nCell * length(flowNames)), nCell, flowNames)
  dstCap <- cellValues(runif(nCell * 3, 0.05, 0.4), nCell, types)
  list(total = total, prior = prior, cluster = cluster, flows = flows, types = types,
       srcBlock = list(group = flows[["from"]], cap = state),
       dstBlock = list(group = flows[["to"]], cap = dstCap))
}

# valid small problem for the input checks
makeSmall <- function() {
  flows <- c("f1", "f2")
  list(total = unitValues(c(1, 2, 3, 4), c("J1", "J2"), flows), prior = cellValues(1, 6, flows),
       cluster = rep(c("J1", "J2"), each = 3), caps = list(list(group = c("g", "g"), cap = cellValues(5, 6, "g"))))
}

# tests --------------------------------------------------------------------------------------

test_that("uncapped case equals the closed form", {
  set.seed(1)
  cluster <- c("J1", "J1", "J1", "J2", "J2")
  prior <- cellValues(runif(15), 5, c("f1", "f2", "f3"))
  total <- unitValues(runif(6), c("J1", "J2"), getNames(prior))
  res <- toolAllocateByCrossEntropy(total, prior, cellMapping(cluster))
  # each flow of a unit in proportion to the prior
  .onCells <- function(xLr) toolAggregate(xLr, cellMapping(cluster), from = "cluster", to = "cell")
  expected <- prior / .onCells(unitSums(prior, cluster)) * .onCells(total)
  expect_equal(as.vector(res[["flowsHr"]]), as.vector(expected), tolerance = 1e-12)
})

test_that("prior scale within a cluster and flow does not matter", {
  tr <- makeCappedFlows()
  caps <- list(tr[["srcBlock"]])
  res1 <- toolAllocateByCrossEntropy(tr[["total"]], tr[["prior"]], cellMapping(tr[["cluster"]]), caps, tol = 1e-12)
  prior2 <- tr[["prior"]]
  cellsJ1 <- paste0("c", which(tr[["cluster"]] == "J1"))
  prior2[cellsJ1, , 2] <- prior2[cellsJ1, , 2] * 7.3
  res2 <- toolAllocateByCrossEntropy(tr[["total"]], prior2, cellMapping(tr[["cluster"]]), caps, tol = 1e-12)
  # the source block binds, so the test is not trivial
  expect_gt(sum(res1[["capFactor"]][[1]] < 1 - 1e-6), 0)
  expect_equal(res2[["flowsHr"]], res1[["flowsHr"]], tolerance = 1e-10)
})

test_that("single capped flow equals water-filling", {
  set.seed(1)
  nCell <- 8
  prior <- runif(nCell)
  cap <- runif(nCell, 0, 0.6)
  total <- 0.8 * sum(cap)
  tLevel <- uniroot(function(t) sum(pmin(cap, t * prior)) - total, c(0, 1e6), tol = 1e-14)[["root"]]
  reference <- pmin(cap, tLevel * prior)
  expect_gte(sum(tLevel * prior > cap), 2)
  res <- toolAllocateByCrossEntropy(unitValues(total, "A", "f"), cellValues(prior, nCell, "f"),
                                    cellMapping(rep("A", nCell)),
                                    list(list(group = "g", cap = cellValues(cap, nCell, "g"))), tol = 1e-12)
  expect_equal(as.vector(res[["flowsHr"]]), reference, tolerance = 1e-10)
})

test_that("two binding cap blocks satisfy the optimality conditions", {
  tr <- makeCappedFlows()
  tol <- 1e-12
  res <- toolAllocateByCrossEntropy(tr[["total"]], tr[["prior"]], cellMapping(tr[["cluster"]]),
                                    list(tr[["srcBlock"]], tr[["dstBlock"]]), tol = tol)
  x <- res[["flowsHr"]]
  beta1 <- res[["capFactor"]][[1]]
  beta2 <- res[["capFactor"]][[2]]
  expect_lt(max(res[["residual"]]), tol)
  expect_gte(sum(beta1 < 1 - 1e-6), 2)
  expect_gte(sum(beta2 < 1 - 1e-6), 2)

  # (a) cluster totals
  expect_equal(as.vector(unitSums(x, tr[["cluster"]])[getCells(tr[["total"]]), , ]), as.vector(tr[["total"]]),
               tolerance = 1e-12)
  # (b) caps: the flows of each group summed per grid cell
  .groupSum <- function(group, cap) {
    out <- cap * 0
    for (k in tr[["types"]]) {
      out[, , k] <- dimSums(x[, , group == k], dim = 3)
    }
    return(out)
  }
  srcSum <- .groupSum(tr[["flows"]][["from"]], tr[["srcBlock"]][["cap"]])
  dstSum <- .groupSum(tr[["flows"]][["to"]], tr[["dstBlock"]][["cap"]])
  expect_lte(max(srcSum - tr[["srcBlock"]][["cap"]]), 1e-9)
  expect_lte(max(dstSum - tr[["dstBlock"]][["cap"]]), 1e-9)
  # (c) beta in [0, 1]
  expect_true(all(beta1 >= 0 & beta1 <= 1 + 1e-12))
  expect_true(all(beta2 >= 0 & beta2 <= 1 + 1e-12))
  # (d) complementary slackness
  expect_lte(max(abs((tr[["srcBlock"]][["cap"]] - srcSum)[beta1 < 1 - 1e-9])), 1e-9)
  expect_lte(max(abs((tr[["dstBlock"]][["cap"]] - dstSum)[beta2 < 1 - 1e-9])), 1e-9)
  # (e) implied alpha is constant within cluster and flow
  .perFlow <- function(beta, group) {
    out <- x * 0
    for (n in seq_along(group)) {
      out[, , n] <- beta[, , group[n]]
    }
    return(out)
  }
  implied <- x / (tr[["prior"]] * .perFlow(beta1, tr[["flows"]][["from"]]) * .perFlow(beta2, tr[["flows"]][["to"]]))
  spread <- NULL
  for (j in getCells(tr[["total"]])) {
    for (flow in getNames(x)) {
      v <- as.vector(implied[paste0("c", which(tr[["cluster"]] == j)), , flow])
      v <- v[is.finite(v)]
      spread <- c(spread, diff(range(v)) / mean(v))
    }
  }
  expect_lt(max(spread), 1e-8)
})

test_that("result does not depend on the order of blocks and cells", {
  tr <- makeCappedFlows()
  res1 <- toolAllocateByCrossEntropy(tr[["total"]], tr[["prior"]], cellMapping(tr[["cluster"]]),
                                     list(tr[["srcBlock"]], tr[["dstBlock"]]), tol = 1e-12)
  set.seed(2)
  perm <- sample(ncells(tr[["prior"]]))
  capsPermuted <- list(list(group = tr[["dstBlock"]][["group"]], cap = tr[["dstBlock"]][["cap"]][perm, , ]),
                       list(group = tr[["srcBlock"]][["group"]], cap = tr[["srcBlock"]][["cap"]][perm, , ]))
  res2 <- toolAllocateByCrossEntropy(tr[["total"]], tr[["prior"]][perm, , ], cellMapping(tr[["cluster"]]),
                                     capsPermuted, tol = 1e-12)
  expect_equal(res2[["flowsHr"]][getCells(res1[["flowsHr"]]), , ], res1[["flowsHr"]], tolerance = 1e-10)
})

test_that("result does not depend on the number of processes", {
  # parallel::mclapply forks, which Windows cannot
  skip_on_os("windows")
  tr <- makeCappedFlows()
  # grid cells of the clusters interleaved, so that the sets of units are not blocks of grid cells
  set.seed(3)
  perm <- sample(ncells(tr[["prior"]]))
  caps <- list(list(group = tr[["srcBlock"]][["group"]], cap = tr[["srcBlock"]][["cap"]][perm, , ]),
               list(group = tr[["dstBlock"]][["group"]], cap = tr[["dstBlock"]][["cap"]][perm, , ]))
  .allocate <- function() {
    return(toolAllocateByCrossEntropy(tr[["total"]], tr[["prior"]][perm, , ], cellMapping(tr[["cluster"]]), caps))
  }
  res1 <- .allocate()
  withr::local_options(mrdownscale.cores = 2)
  expect_identical(.allocate(), res1)
})

test_that("the option mrdownscale.cores is checked", {
  small <- makeSmall()
  withr::local_options(mrdownscale.cores = 1.5)
  expect_error(toolAllocateByCrossEntropy(small[["total"]], small[["prior"]], cellMapping(small[["cluster"]]),
                                          small[["caps"]]), "mrdownscale.cores")
})

test_that("zero totals and zero caps give zero flows", {
  flows <- c("a>b", "a>c", "b>a")
  total <- unitValues(c(0.5, 0, 0.3), "J1", flows)
  cap <- cellValues(1, 5, c("a", "b"))
  cap["c2", , "a"] <- 0
  res <- toolAllocateByCrossEntropy(total, cellValues(1, 5, flows), cellMapping(rep("J1", 5)),
                                    list(list(group = c("a", "a", "b"), cap = cap)), tol = 1e-12)
  expect_equal(as.vector(res[["flowsHr"]][, , "a>c"]), rep(0, 5))
  expect_equal(as.vector(res[["flowsHr"]]["c2", , c("a>b", "a>c")]), c(0, 0))
  expect_equal(as.vector(dimSums(res[["flowsHr"]], dim = 1)), as.vector(total), tolerance = 1e-12)
})

test_that("structural zeros use the fallback", {
  flows <- c("f1", "f2")
  total <- unitValues(c(1, 2, 3, 4), c("J1", "J2"), flows)
  prior <- cellValues(1, 4, flows)
  prior[c("c3", "c4"), , "f2"] <- 0
  cluster <- c("J1", "J1", "J2", "J2")
  fallback <- cellValues(c(1, 1, 1, 1, 1, 1, 1, 3), 4, flows)
  res <- toolAllocateByCrossEntropy(total, prior, cellMapping(cluster), fallbackHr = fallback)
  expect_equal(as.vector(res[["flowsHr"]][c("c3", "c4"), , "f2"]), 4 * c(1, 3) / 4, tolerance = 1e-12)
  expect_equal(as.vector(res[["flowsHr"]][c("c1", "c2"), , "f2"]), c(1.5, 1.5), tolerance = 1e-12)
  expect_error(toolAllocateByCrossEntropy(total, prior, cellMapping(cluster)), "no fallback")
})

test_that("infeasible caps show as residual", {
  flows <- c("a>b", "a>c")
  cluster <- rep(c("J1", "J2"), each = 4)
  total <- unitValues(c(0.5, 0.3, 0.5, 0.3), c("J1", "J2"), flows)
  cap <- cellValues(c(rep(0.2, 4), rep(1, 4)), 8, "a")
  expect_warning(
    res <- toolAllocateByCrossEntropy(total, cellValues(1, 8, flows), cellMapping(cluster),
                                      list(list(group = c("a", "a"), cap = cap)), maxIter = 200),
    "did not converge"
  )
  expect_gte(as.vector(res[["residual"]]["J1", , ]), 1e-8)
  expect_lt(as.vector(res[["residual"]]["J2", , ]), 1e-8)
  expect_equal(as.vector(unitSums(res[["flowsHr"]], cluster)[getCells(total), , ]), as.vector(total),
               tolerance = 1e-12)
})

test_that("soft caps hold where feasible and are exceeded at least cost where not", {
  cluster <- rep("J1", 3)
  prior <- cellValues(1, 3, "f")
  soft <- list(group = "g", cap = cellValues(c(0.2, 0.5, 1), 3, "g"), minFactor = 1e-6)
  total <- unitValues(1.5, "J1", "f")
  res <- toolAllocateByCrossEntropy(total, prior, cellMapping(cluster), list(soft), maxIter = 1e5)
  expect_equal(as.vector(res[["flowsHr"]]), c(0.2, 0.5, 0.8), tolerance = 1e-5)
  # the caps sum to 1.7 < 2: the excess lifts the lowest cells to a common level
  total[, , ] <- 2
  res <- toolAllocateByCrossEntropy(total, prior, cellMapping(cluster), list(soft), maxIter = 1e5)
  expect_equal(as.vector(res[["flowsHr"]]), c(0.5, 0.5, 1), tolerance = 1e-5)
  expect_lt(max(res[["residual"]]), 1e-8)
})

test_that("a cap on a net change equals the closed form", {
  # inflow and outflow of one type in two cells; net change of cell 1 at most -0.5
  cluster <- c("J1", "J1")
  total <- unitValues(1, "J1", c("in", "out"))
  cap <- list(group = c("g", "g"), sign = c(1, -1), cap = cellValues(c(-0.5, Inf), 2, "g"))
  res <- toolAllocateByCrossEntropy(total, cellValues(1, 2, c("in", "out")), cellMapping(cluster), list(cap),
                                    maxIter = 1e4)
  # in = a / (1 + a) and out = 1 / (1 + a) in cell 1, with (a - 1) / (a + 1) = -0.5, so a = 1 / 3
  expect_equal(as.vector(res[["flowsHr"]]["c1", , ]), c(0.25, 0.75), tolerance = 1e-8)
  expect_equal(as.vector(dimSums(res[["flowsHr"]], dim = 1)), c(1, 1), tolerance = 1e-12)
})

test_that("input errors", {
  args <- makeSmall()
  run <- function(a) {
    toolAllocateByCrossEntropy(a[["total"]], a[["prior"]], cellMapping(a[["cluster"]]), a[["caps"]])
  }
  expect_silent(run(args))

  bad <- args
  bad[["total"]] <- as.array(bad[["total"]])[, 1, ]
  expect_error(run(bad), "totalLr must be a magclass object")
  bad <- args
  bad[["total"]]["J1", , "f1"] <- -1
  expect_error(run(bad), "totalLr contains negative values")
  bad <- args
  bad[["prior"]]["c1", , "f1"] <- NA
  expect_error(run(bad), "priorHr contains NA")
  bad <- args
  getNames(bad[["prior"]]) <- c("f1", "x")
  expect_error(run(bad), "flows of totalLr")
  bad <- args
  bad[["cluster"]][1] <- "J9"
  expect_error(run(bad), "not in totalLr")
  bad <- args
  bad[["cluster"]] <- rep("J1", 6)
  expect_error(run(bad), "at least one grid cell")
  bad <- args
  bad[["caps"]][[1]][["group"]] <- "g"
  expect_error(run(bad), "one entry per flow")
  bad <- args
  bad[["caps"]][[1]][["group"]] <- c("g", "h")
  expect_error(run(bad), "needs an item in cap")
  bad <- args
  bad[["caps"]][[1]][["sign"]] <- c(1, 0)
  expect_error(run(bad), "sign needs one entry per flow")
  bad <- args
  bad[["caps"]][[1]][["minFactor"]] <- 1
  expect_error(run(bad), "minFactor must be one value")
})

# helpers ------------------------------------------------------------------------------------

# small matrix with column names
mat <- function(values, nrow, colnames) {
  matrix(values, nrow = nrow, ncol = length(colnames), dimnames = list(paste0("c", seq_len(nrow)), colnames))
}

# one-column cap matrix
capMat <- function(values, group) {
  matrix(values, ncol = 1, dimnames = list(NULL, group))
}

# transitions between three land types in two clusters with source and destination caps
makeTransitions <- function(seed = 2) {
  set.seed(seed)
  types <- c("a", "b", "c")
  flows <- expand.grid(from = types, to = types, stringsAsFactors = FALSE)
  flows <- flows[flows[["from"]] != flows[["to"]], ]
  flowNames <- paste(flows[["from"]], flows[["to"]], sep = ">")
  nCell <- 12
  cluster <- rep(c("J1", "J2"), each = 6)
  state <- matrix(runif(nCell * 3), nCell, 3, dimnames = list(NULL, types))
  total <- matrix(0, 2, length(flowNames), dimnames = list(c("J1", "J2"), flowNames))
  for (j in 1:2) {
    stateSum <- colSums(state[cluster == rownames(total)[j], ])
    total[j, ] <- runif(length(flowNames)) * 0.15 * stateSum[flows[["from"]]]
  }
  prior <- matrix(runif(nCell * length(flowNames)), nCell, dimnames = list(paste0("c", 1:nCell), flowNames))
  dstCap <- matrix(runif(nCell * 3, 0.05, 0.4), nCell, 3, dimnames = list(NULL, types))
  list(total = total, prior = prior, cluster = cluster, flows = flows, types = types,
       srcBlock = list(group = flows[["from"]], cap = state),
       dstBlock = list(group = flows[["to"]], cap = dstCap))
}

# valid small problem for the input checks
makeSmall <- function() {
  flows <- c("f1", "f2")
  list(total = matrix(c(1, 2, 3, 4), 2, dimnames = list(c("J1", "J2"), flows)),
       prior = mat(1, 6, flows), cluster = rep(c("J1", "J2"), each = 3),
       caps = list(list(group = c("g", "g"), cap = capMat(rep(5, 6), "g"))))
}

# tests --------------------------------------------------------------------------------------

test_that("uncapped case equals the closed form", {
  set.seed(1)
  cluster <- c("J1", "J1", "J1", "J2", "J2")
  prior <- mat(runif(15), 5, c("f1", "f2", "f3"))
  total <- matrix(runif(6), 2, dimnames = list(c("J1", "J2"), colnames(prior)))
  res <- toolAllocateByCrossEntropy(total, prior, cellMapping(cluster))
  expected <- prior
  for (j in rownames(total)) {
    cells <- cluster == j
    expected[cells, ] <- sweep(prior[cells, ], 2, colSums(prior[cells, ]), "/") *
      rep(total[j, ], each = sum(cells))
  }
  expect_equal(res[["flowsHr"]], expected, tolerance = 1e-12)
})

test_that("prior scale within a cluster and flow does not matter", {
  tr <- makeTransitions()
  caps <- list(tr[["srcBlock"]])
  res1 <- toolAllocateByCrossEntropy(tr[["total"]], tr[["prior"]], cellMapping(tr[["cluster"]]), caps, tol = 1e-12)
  prior2 <- tr[["prior"]]
  prior2[tr[["cluster"]] == "J1", 2] <- prior2[tr[["cluster"]] == "J1", 2] * 7.3
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
  res <- toolAllocateByCrossEntropy(matrix(total, 1, 1, dimnames = list("A", "f")), mat(prior, nCell, "f"),
                                    cellMapping(rep("A", nCell)), list(list(group = "g", cap = capMat(cap, "g"))),
                                    tol = 1e-12)
  expect_equal(unname(res[["flowsHr"]][, 1]), reference, tolerance = 1e-10)
})

test_that("two binding cap blocks satisfy the optimality conditions", {
  tr <- makeTransitions()
  tol <- 1e-12
  res <- toolAllocateByCrossEntropy(tr[["total"]], tr[["prior"]], cellMapping(tr[["cluster"]]),
                                    list(tr[["srcBlock"]], tr[["dstBlock"]]),
                                    tol = tol)
  x <- res[["flowsHr"]]
  beta1 <- res[["capFactor"]][[1]]
  beta2 <- res[["capFactor"]][[2]]
  expect_lt(max(res[["residual"]]), tol)
  expect_gte(sum(beta1 < 1 - 1e-6), 2)
  expect_gte(sum(beta2 < 1 - 1e-6), 2)

  # (a) cluster totals
  clusterIdx <- match(tr[["cluster"]], rownames(tr[["total"]]))
  expect_equal(unname(rowsum(x, clusterIdx)), unname(tr[["total"]]), tolerance = 1e-12)
  # (b) caps
  srcSum <- vapply(tr[["types"]], function(k) rowSums(x[, tr[["flows"]][["from"]] == k, drop = FALSE]),
                   numeric(nrow(x)))
  dstSum <- vapply(tr[["types"]], function(k) rowSums(x[, tr[["flows"]][["to"]] == k, drop = FALSE]), numeric(nrow(x)))
  expect_lte(max(srcSum - tr[["srcBlock"]][["cap"]]), 1e-9)
  expect_lte(max(dstSum - tr[["dstBlock"]][["cap"]]), 1e-9)
  # (c) beta in [0, 1]
  expect_true(all(beta1 >= 0 & beta1 <= 1 + 1e-12))
  expect_true(all(beta2 >= 0 & beta2 <= 1 + 1e-12))
  # (d) complementary slackness
  expect_lte(max(abs((tr[["srcBlock"]][["cap"]] - srcSum)[beta1 < 1 - 1e-9])), 1e-9)
  expect_lte(max(abs((tr[["dstBlock"]][["cap"]] - dstSum)[beta2 < 1 - 1e-9])), 1e-9)
  # (e) implied alpha is constant within cluster and flow
  implied <- x / (tr[["prior"]] * beta1[, tr[["flows"]][["from"]]] * beta2[, tr[["flows"]][["to"]]])
  spread <- unlist(lapply(rownames(tr[["total"]]), function(j) {
    apply(implied[tr[["cluster"]] == j, , drop = FALSE], 2, function(v) {
      v <- v[is.finite(v)]
      diff(range(v)) / mean(v)
    })
  }))
  expect_lt(max(spread), 1e-8)
})

test_that("result does not depend on the order of blocks and cells", {
  tr <- makeTransitions()
  res1 <- toolAllocateByCrossEntropy(tr[["total"]], tr[["prior"]], cellMapping(tr[["cluster"]]),
                                     list(tr[["srcBlock"]], tr[["dstBlock"]]),
                                     tol = 1e-12)
  set.seed(2)
  perm <- sample(nrow(tr[["prior"]]))
  capsPermuted <- list(list(group = tr[["dstBlock"]][["group"]], cap = tr[["dstBlock"]][["cap"]][perm, ]),
                       list(group = tr[["srcBlock"]][["group"]], cap = tr[["srcBlock"]][["cap"]][perm, ]))
  res2 <- toolAllocateByCrossEntropy(tr[["total"]], tr[["prior"]][perm, ], cellMapping(tr[["cluster"]]), capsPermuted,
                                     tol = 1e-12)
  expect_equal(res2[["flowsHr"]][order(perm), ], res1[["flowsHr"]], tolerance = 1e-10)
})

test_that("zero totals and zero caps give zero flows", {
  flows <- c("a>b", "a>c", "b>a")
  total <- matrix(c(0.5, 0, 0.3), 1, dimnames = list("J1", flows))
  prior <- mat(1, 5, flows)
  cap <- matrix(1, 5, 2, dimnames = list(NULL, c("a", "b")))
  cap[2, "a"] <- 0
  res <- toolAllocateByCrossEntropy(total, prior, cellMapping(rep("J1", 5)),
                                    list(list(group = c("a", "a", "b"), cap = cap)), tol = 1e-12)
  expect_equal(unname(res[["flowsHr"]][, "a>c"]), rep(0, 5))
  expect_equal(unname(res[["flowsHr"]][2, c("a>b", "a>c")]), c(0, 0))
  expect_equal(unname(colSums(res[["flowsHr"]])), unname(total[1, ]), tolerance = 1e-12)
})

test_that("structural zeros use the fallback", {
  flows <- c("f1", "f2")
  total <- matrix(c(1, 2, 3, 4), 2, dimnames = list(c("J1", "J2"), flows))
  prior <- mat(1, 4, flows)
  prior[3:4, "f2"] <- 0
  cluster <- c("J1", "J1", "J2", "J2")
  fallback <- mat(c(1, 1, 1, 1, 1, 1, 1, 3), 4, flows)
  res <- toolAllocateByCrossEntropy(total, prior, cellMapping(cluster), fallbackHr = fallback)
  expect_equal(unname(res[["flowsHr"]][3:4, "f2"]), 4 * c(1, 3) / 4, tolerance = 1e-12)
  expect_equal(unname(res[["flowsHr"]][1:2, "f2"]), c(1.5, 1.5), tolerance = 1e-12)
  expect_error(toolAllocateByCrossEntropy(total, prior, cellMapping(cluster)), "no fallback")
})

test_that("infeasible caps show as residual", {
  flows <- c("a>b", "a>c")
  cluster <- rep(c("J1", "J2"), each = 4)
  total <- matrix(c(0.5, 0.3, 0.5, 0.3), 2, dimnames = list(c("J1", "J2"), flows))
  cap <- capMat(c(rep(0.2, 4), rep(1, 4)), "a")
  expect_warning(
    res <- toolAllocateByCrossEntropy(total, mat(1, 8, flows), cellMapping(cluster),
                                      list(list(group = c("a", "a"), cap = cap)), maxIter = 200),
    "did not converge"
  )
  expect_gte(res[["residual"]][["J1"]], 1e-8)
  expect_lt(res[["residual"]][["J2"]], 1e-8)
  clusterIdx <- match(cluster, rownames(total))
  expect_equal(unname(rowsum(res[["flowsHr"]], clusterIdx)), unname(total), tolerance = 1e-12)
})

test_that("soft caps hold where feasible and are exceeded at least cost where not", {
  cluster <- rep("J1", 3)
  prior <- mat(1, 3, "f")
  soft <- list(group = rep("g", 1), cap = capMat(c(0.2, 0.5, 1), "g"), minFactor = 1e-6)
  total <- matrix(1.5, 1, dimnames = list("J1", "f"))
  res <- toolAllocateByCrossEntropy(total, prior, cellMapping(cluster), list(soft), maxIter = 1e5)
  expect_equal(unname(res[["flowsHr"]][, 1]), c(0.2, 0.5, 0.8), tolerance = 1e-5)
  # the caps sum to 1.7 < 2: the excess lifts the lowest cells to a common level
  total[1, 1] <- 2
  res <- toolAllocateByCrossEntropy(total, prior, cellMapping(cluster), list(soft), maxIter = 1e5)
  expect_equal(unname(res[["flowsHr"]][, 1]), c(0.5, 0.5, 1), tolerance = 1e-5)
  expect_lt(res[["residual"]][["J1"]], 1e-8)
})

test_that("a cap on a net change equals the closed form", {
  # inflow and outflow of one type in two cells; net change of cell 1 at most -0.5
  cluster <- c("J1", "J1")
  total <- matrix(1, 1, 2, dimnames = list("J1", c("in", "out")))
  cap <- list(group = c("g", "g"), sign = c(1, -1), cap = capMat(c(-0.5, Inf), "g"))
  res <- toolAllocateByCrossEntropy(total, mat(1, 2, c("in", "out")), cellMapping(cluster), list(cap), maxIter = 1e4)
  # in = a / (1 + a) and out = 1 / (1 + a) in cell 1, with (a - 1) / (a + 1) = -0.5, so a = 1 / 3
  expect_equal(unname(res[["flowsHr"]][1, ]), c(0.25, 0.75), tolerance = 1e-8)
  expect_equal(unname(colSums(res[["flowsHr"]])), c(1, 1), tolerance = 1e-12)
})

test_that("input errors", {
  args <- makeSmall()
  run <- function(a) toolAllocateByCrossEntropy(a[["total"]], a[["prior"]], cellMapping(a[["cluster"]]), a[["caps"]])
  expect_silent(run(args))

  bad <- args
  bad[["total"]][1, 1] <- -1
  expect_error(run(bad), "totalLr contains negative values")
  bad <- args
  bad[["prior"]][1, 1] <- NA
  expect_error(run(bad), "priorHr contains NA")
  bad <- args
  colnames(bad[["prior"]]) <- c("f1", "x")
  expect_error(run(bad), "column names of totalLr")
  bad <- args
  bad[["cluster"]][1] <- "J9"
  expect_error(run(bad), "not row names of totalLr")
  bad <- args
  bad[["cluster"]] <- rep("J1", 6)
  expect_error(run(bad), "at least one grid cell")
  bad <- args
  bad[["caps"]][[1]][["group"]] <- "g"
  expect_error(run(bad), "one entry per flow")
  bad <- args
  bad[["caps"]][[1]][["group"]] <- c("g", "h")
  expect_error(run(bad), "needs a column in cap")
  bad <- args
  bad[["caps"]][[1]][["sign"]] <- c(1, 0)
  expect_error(run(bad), "sign needs one entry per flow")
  bad <- args
  bad[["caps"]][[1]][["minFactor"]] <- 1
  expect_error(run(bad), "minFactor must be one value")
})

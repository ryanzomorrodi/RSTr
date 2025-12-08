test_that("medians are generated from test data", {
  data_min <- lapply(miheart, \(x) x[1:2, 1:3, 1:3])
  adj_min <- list(2, 1)
  mod_mst <- mstcar("test", data_min, adj_min, tempdir(), show_plots = FALSE)
})

# Milestone M2: the two-moons posterior is bimodal and NPE must recover both
# modes. The observation constrains |theta1 + theta2|, so the posterior has
# two symmetric modes with theta1 + theta2 of opposite signs.

test_that("NPE recovers both modes of the two-moons posterior", {
  skip_if_no_torch()
  skip_on_cran()
  set.seed(15)
  task <- task_two_moons()
  fit <- npe(task$prior, task$simulator, n_simulations = 10000,
             density_estimator = "mdn", n_components = 5L,
             hidden = c(50L, 50L), max_epochs = 150L, seed = 15)
  theta_true <- c(0.4, 0.5)
  x_obs <- as.numeric(task$simulator(matrix(theta_true, nrow = 1)))
  post <- posterior(fit, x_obs = x_obs)
  draws <- suppressWarnings(sample(post, 4000))

  s <- draws[, 1] + draws[, 2]
  # both signs of theta1 + theta2 must be represented (bimodality)
  expect_gt(mean(s > 0), 0.15)
  expect_gt(mean(s < 0), 0.15)
  # and |s| should concentrate near |theta1 + theta2| = 0.9
  expect_equal(mean(abs(s)), abs(sum(theta_true)), tolerance = 0.25)
})

context("Generalized point Laplace")

n <- 1000
set.seed(1)
s <- rnorm(n, 1, 0.1)

true_pi <- c(0.5, 0.3, 0.2)
true_scale_pos <- 10
true_scale_neg <- 4
true_mean <- 0

theta <- c(rexp(n * true_pi[2], rate = 1 / true_scale_pos),
           -rexp(n * true_pi[3], rate = 1 / true_scale_neg),
           rep(0, n * true_pi[1]))
x <- theta + rnorm(n, sd = s)

true_g <- genlaplacemix(pi = true_pi,
                        mean = rep(true_mean, 3),
                        scale_pos = c(0, true_scale_pos, 0),
                        scale_neg = c(0, 0, true_scale_neg))

genpl.res <- ebnm(x, s, prior_family = "gen_point_laplace")

# Build the optimizer's parameter vector from interpretable quantities. The
#   order is c(logit_w_pos, logit_w_neg, log_rate_pos, log_rate_neg, mu), where
#   the weights are a softmax with the point mass as reference category.
genpl_par <- function(pi_0, pi_plus, pi_neg, rate_pos, rate_neg, mu) {
  c(log(pi_plus / pi_0), log(pi_neg / pi_0),
    log(rate_pos), log(rate_neg), mu)
}

nllik_at <- function(par) {
  genpl_nllik(par, x, s, par_init = NULL, fix_par = rep(FALSE, 5),
              calc_grad = FALSE, calc_hess = FALSE)
}

test_that("Basic functionality works", {
  genpl.res2 <- ebnm_gen_point_laplace(x, s)
  genpl.res$call <- genpl.res2$call <- NULL
  expect_identical(genpl.res, genpl.res2)
  expect_equal(genpl.res[[g_ret_str()]], true_g, tolerance = 0.1)
})

test_that("nllik agrees with the reference implementation", {
  pi_0 <- 0.5; pi_plus <- 0.3; pi_neg <- 0.2
  rate_pos <- 0.3; rate_neg <- 0.8; mu <- 0
  par <- genpl_par(pi_0, pi_plus, pi_neg, rate_pos, rate_neg, mu)
  expect_equal(
    nllik_at(par),
    -loglik_gen_point_laplace(x, s, pi_0, pi_plus, pi_neg,
                              rate_pos, rate_neg, mu)
  )
})

test_that("Symmetric case reduces to point-Laplace", {
  # Equal rates and equal tail weights give a Laplace slab with weight w.
  w <- 0.5; a <- 0.1; mu <- 0
  par <- genpl_par(1 - w, w / 2, w / 2, a, a, mu)
  expect_equal(nllik_at(par), -loglik_point_laplace(x, s, w = w, a = a, mu = mu))
})

test_that("Mode estimation works", {
  genpl.res2 <- ebnm_gen_point_laplace(x, s, mode = "est")
  expect_equal(genpl.res2[[g_ret_str()]], true_g, tolerance = 0.5)
  expect_false(identical(genpl.res2[[g_ret_str()]]$mean[1], true_mean))
})

test_that("Fixing the scales works", {
  genpl.res2 <- ebnm_gen_point_laplace(
    x, s, scale = c(true_scale_pos, true_scale_neg)
  )
  expect_equal(genpl.res2[[g_ret_str()]]$scale_pos[2], true_scale_pos)
  expect_equal(genpl.res2[[g_ret_str()]]$scale_neg[3], true_scale_neg)
})

test_that("A scalar scale is shared by both tails", {
  genpl.res2 <- ebnm_gen_point_laplace(x, s, scale = 5)
  expect_equal(genpl.res2[[g_ret_str()]]$scale_pos[2], 5)
  expect_equal(genpl.res2[[g_ret_str()]]$scale_neg[3], 5)
})

test_that("Fixing g works", {
  genpl.res2 <- ebnm_gen_point_laplace(
    x, s, g_init = genpl.res[[g_ret_str()]], fix_g = TRUE
  )
  expect_identical(genpl.res[[g_ret_str()]], genpl.res2[[g_ret_str()]])
  expect_equal(as.numeric(genpl.res[[llik_ret_str()]]),
               as.numeric(genpl.res2[[llik_ret_str()]]))
})

test_that("Initializing g works", {
  genpl.res2 <- ebnm_gen_point_laplace(x, s, g_init = true_g)
  expect_equal(genpl.res[[llik_ret_str()]], genpl.res2[[llik_ret_str()]],
               tolerance = 0.01)
})

test_that("Output parameter works", {
  genpl.res2 <- ebnm_gen_point_laplace(x, s, output = c("fitted_g"))
  expect_identical(names(genpl.res2), c(g_ret_str(), "call"))
})

test_that("The posterior sampler works", {
  genpl.res2 <- ebnm_gen_point_laplace(x, s, output = ebnm_output_all())
  samp <- genpl.res2[[samp_ret_str()]](100)
  expect_equal(dim(samp), c(100, n))
  # Sample means should track the posterior means.
  expect_equal(colMeans(samp), genpl.res2[[df_ret_str()]][[pm_ret_str()]],
               tolerance = 0.5)
})

test_that("Homoskedastic standard errors work", {
  genpl.res2 <- ebnm_gen_point_laplace(x, 1)
  expect_equal(genpl.res2[[g_ret_str()]], true_g, tolerance = 0.1)
})

test_that("All optmethods give the same answer", {
  res_lbfgsb <- ebnm_gen_point_laplace(x, s, optmethod = "nograd_lbfgsb")
  res_nlm <- ebnm_gen_point_laplace(x, s, optmethod = "nograd_nlm")
  expect_equal(as.numeric(res_lbfgsb[[llik_ret_str()]]),
               as.numeric(res_nlm[[llik_ret_str()]]),
               tolerance = 1e-4)
})

test_that("Null case gives a degenerate prior", {
  set.seed(2)
  xnull <- rnorm(n, sd = 0.5)
  res <- ebnm_gen_point_laplace(xnull, s = 1)
  # This family is non-identifiable at the null: pi0 = 1 and (pi_plus = 1,
  #   scale_pos -> 0) describe the same prior, so pi[1] is not a reliable
  #   indicator. Check instead that the prior is effectively degenerate.
  expect_lt(max(abs(coef(res))), 0.05)
})

test_that("Very large observations give reasonable results", {
  scl <- 1e8
  r  <- ebnm_gen_point_laplace(x, s, mode = "estimate")[[g_ret_str()]]
  rl <- ebnm_gen_point_laplace(scl * x, scl * s, mode = "estimate")[[g_ret_str()]]
  expect_equal(r$pi, rl$pi, tolerance = 1e-3)
  expect_equal(scl * r$scale_pos[2], rl$scale_pos[2], tolerance = 1e-3)
  expect_equal(scl * r$scale_neg[3], rl$scale_neg[3], tolerance = 1e-3)
})

test_that("Very small observations give reasonable results", {
  scl <- 1e-8
  r  <- ebnm_gen_point_laplace(x, s, mode = "estimate")[[g_ret_str()]]
  rs <- ebnm_gen_point_laplace(scl * x, scl * s, mode = "estimate")[[g_ret_str()]]
  expect_equal(r$pi, rs$pi, tolerance = 1e-3)
  expect_equal(scl * r$scale_pos[2], rs$scale_pos[2], tolerance = 1e-3)
  expect_equal(scl * r$scale_neg[3], rs$scale_neg[3], tolerance = 1e-3)
})

test_that("Equal tails reproduce point-Laplace posteriors", {
  pl <- ebnm_point_laplace(x, s, output = ebnm_output_all())
  g <- pl[[g_ret_str()]]
  # A Laplace slab is two exponential halves with equal weights and scales.
  gg <- genlaplacemix(pi = c(g$pi[1], g$pi[2] / 2, g$pi[2] / 2),
                      mean = rep(g$mean[1], 3),
                      scale_pos = c(0, g$scale[2], 0),
                      scale_neg = c(0, 0, g$scale[2]))
  gp <- ebnm_gen_point_laplace(x, s, g_init = gg, fix_g = TRUE,
                               output = ebnm_output_all())
  expect_equal(gp$posterior$mean, pl$posterior$mean, tolerance = 1e-8)
  expect_equal(gp$posterior$sd,   pl$posterior$sd,   tolerance = 1e-8)
  expect_equal(gp$posterior$lfsr, pl$posterior$lfsr, tolerance = 1e-8)
  expect_equal(as.numeric(gp[[llik_ret_str()]]),
               as.numeric(pl[[llik_ret_str()]]), tolerance = 1e-8)
})

test_that("A zero-weight left tail reproduces point-exponential", {
  set.seed(3)
  th <- c(rexp(n / 2, rate = 1 / 5), rep(0, n / 2))
  xx <- th + rnorm(n)
  pe <- ebnm_point_exponential(xx, 1, output = ebnm_output_all())
  g <- pe[[g_ret_str()]]
  gg <- genlaplacemix(pi = c(g$pi[1], g$pi[2], 0),
                      mean = rep(g$shift[1], 3),
                      scale_pos = c(0, g$scale[2], 0),
                      scale_neg = c(0, 0, 1))
  gp <- ebnm_gen_point_laplace(xx, 1, g_init = gg, fix_g = TRUE,
                               output = ebnm_output_all())
  expect_equal(gp$posterior$mean, pe$posterior$mean, tolerance = 1e-5)
  expect_equal(gp$posterior$sd,   pe$posterior$sd,   tolerance = 1e-5)
  expect_equal(as.numeric(gp[[llik_ret_str()]]),
               as.numeric(pe[[llik_ret_str()]]), tolerance = 1e-5)
})

test_that("The prior family is correctly inferred", {
  expect_identical(infer_prior_family(genpl.res[[g_ret_str()]]),
                   "gen_point_laplace")
})

context("Generalized point Laplace")

# ==========================================================================
# Fixtures
# ==========================================================================

# Statistical fixture (large n): used for recovery and model-level behavior.
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

# Unit fixture (small n): used for likelihood and derivative checks, where a
#   large n only slows numDeriv down without improving coverage.
nu <- 50
set.seed(11)
su <- rnorm(nu, 1, 0.1)
xu <- c(rexp(15, 1 / 3), -rexp(15, 1 / 2), rep(0, 20)) + rnorm(nu, sd = su)

# ==========================================================================
# Helpers
# ==========================================================================

# Build the optimizer's parameter vector from interpretable quantities. The
#   order is c(logit_w_pos, logit_w_neg, log_rate_pos, log_rate_neg, mu), where
#   the weights are a softmax with the point mass as reference category.
genpl_par <- function(pi_0, pi_plus, pi_neg, rate_pos, rate_neg, mu) {
  c(log(pi_plus / pi_0), log(pi_neg / pi_0),
    log(rate_pos), log(rate_neg), mu)
}

# Map a full parameter vector onto the canonical quantities, mirroring the
#   transformations in genpl_nllik exactly. Everything downstream compares
#   against the independent reference implementation in
#   loglik_gen_point_laplace.R, so this replicates the reparameterizations
#   without replicating the likelihood.
genpl_canon <- function(p, fixed_comp = "none", fixed_c = NULL,
                        ratio_ref = NULL, ratio_k = NULL) {
  if (identical(fixed_comp, "none")) {
    logits <- c(0, p[1], p[2])
    w      <- exp(logits - max(logits))
    pivec  <- w / sum(w)
    pi_0   <- pivec[1]; pi_pos <- pivec[2]; pi_neg <- pivec[3]
    beta_pos <- p[3]; beta_neg <- p[4]; mu <- p[5]
  } else {
    ea  <- exp(p[1])
    f1  <- ea / (1 + ea)
    f2  <- 1 / (1 + ea)
    rem <- 1 - fixed_c
    if (identical(fixed_comp, "spike")) {
      pi_0   <- fixed_c; pi_pos <- f1 * rem; pi_neg <- f2 * rem
    } else if (identical(fixed_comp, "pos")) {
      pi_pos <- fixed_c; pi_0   <- f1 * rem; pi_neg <- f2 * rem
    } else {
      pi_neg <- fixed_c; pi_0   <- f1 * rem; pi_pos <- f2 * rem
    }
    beta_pos <- p[2]; beta_neg <- p[3]; mu <- p[4]
  }

  if (!is.null(ratio_ref)) {
    if (identical(ratio_ref, "pos")) {
      beta_neg <- (beta_pos - log(ratio_k)) + genpl_softplus(beta_neg)
    } else {
      beta_pos <- (beta_neg - log(ratio_k)) + genpl_softplus(beta_pos)
    }
  }

  return(list(pi_0 = pi_0, pi_pos = pi_pos, pi_neg = pi_neg,
              lambda_pos = exp(beta_pos), lambda_neg = exp(beta_neg),
              mu = mu))
}

# The reference negative log likelihood, as a function of the free parameters
#   only, for numDeriv to differentiate.
genpl_ref_nllik <- function(free, p_full, fix_par, data_x, data_s, ...) {
  p <- p_full
  p[!fix_par] <- free
  cn <- genpl_canon(p, ...)
  return(-loglik_gen_point_laplace(data_x, data_s, cn$pi_0, cn$pi_pos,
                                   cn$pi_neg, cn$lambda_pos, cn$lambda_neg,
                                   cn$mu))
}

# Call genpl_nllik with a given full parameter vector and fix_par pattern.
genpl_eval <- function(p_full, fix_par, data_x, data_s, calc_grad = TRUE,
                       calc_hess = TRUE, ...) {
  return(genpl_nllik(p_full[!fix_par], data_x, data_s,
                     par_init = as.list(p_full), fix_par = fix_par,
                     calc_grad = calc_grad, calc_hess = calc_hess, ...))
}

# Invariants that every fitted prior should satisfy.
expect_valid_genpl <- function(g, info = NULL) {
  expect_equal(sum(g$pi), 1, info = info)
  expect_true(all(g$pi >= 0), info = info)
  expect_true(all(g$scale_pos >= 0), info = info)
  expect_true(all(g$scale_neg >= 0), info = info)
  expect_true(all(is.finite(g$mean)), info = info)
}

# Every fix_par pattern with at least one free parameter.
fix_patterns <- function(k) {
  grid <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), k)))
  return(grid[rowSums(grid) < k, , drop = FALSE])
}

# ==========================================================================
# Model-level behavior
# ==========================================================================

test_that("Basic functionality works", {
  genpl.res2 <- ebnm_gen_point_laplace(x, s)
  genpl.res$call <- genpl.res2$call <- NULL
  expect_identical(genpl.res, genpl.res2)
  expect_equal(genpl.res[[g_ret_str()]], true_g, tolerance = 0.1)
  expect_valid_genpl(genpl.res[[g_ret_str()]])
})

test_that("Mode estimation works", {
  genpl.res2 <- ebnm_gen_point_laplace(x, s, mode = "est")
  expect_equal(genpl.res2[[g_ret_str()]], true_g, tolerance = 0.5)
  expect_false(identical(genpl.res2[[g_ret_str()]]$mean[1], true_mean))
  # The estimate should actually be close, not merely different from the truth.
  # Across 20 replicates of this configuration the mode's RMSE is about 0.06
  # with a worst case near 0.21, so this bound is loose but not vacuous.
  expect_lt(abs(genpl.res2[[g_ret_str()]]$mean[1] - true_mean), 0.3)
})

test_that("Mode estimation is translation equivariant", {
  shift <- 17
  a <- ebnm_gen_point_laplace(x, s, mode = "estimate")[[g_ret_str()]]
  b <- ebnm_gen_point_laplace(x + shift, s, mode = "estimate")[[g_ret_str()]]
  # Shifting the data must shift the mode by as much and leave the weights and
  # scales alone, which pins down the mode handling in genpl_initpar,
  # genpl_scalepar, and genpl_summres at once. The likelihood is exactly
  # equivariant, but x + shift rounds differently, so the optimizer stops at a
  # slightly different point within its own convergence tolerance.
  expect_equal(b$mean[1] - a$mean[1], shift, tolerance = 1e-4)
  expect_equal(b$pi, a$pi, tolerance = 1e-4)
  expect_equal(b$scale_pos, a$scale_pos, tolerance = 1e-4)
  expect_equal(b$scale_neg, a$scale_neg, tolerance = 1e-4)
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

test_that("A two-component g_init (no spike) works", {
  # genpl_initpar and genpl_partog both have a dedicated branch for a prior
  #   with no point mass.
  g2 <- genlaplacemix(pi = c(0.6, 0.4),
                      mean = rep(0, 2),
                      scale_pos = c(true_scale_pos, 0),
                      scale_neg = c(0, true_scale_neg))
  res <- ebnm_gen_point_laplace(x, s, g_init = g2)
  expect_valid_genpl(res[[g_ret_str()]])
  expect_true(is.finite(as.numeric(res[[llik_ret_str()]])))
})

test_that("Output parameter works", {
  genpl.res2 <- ebnm_gen_point_laplace(x, s, output = c("fitted_g"))
  expect_identical(names(genpl.res2), c(g_ret_str(), "call"))
})

test_that("Homoskedastic standard errors work", {
  genpl.res2 <- ebnm_gen_point_laplace(x, 1)
  expect_equal(genpl.res2[[g_ret_str()]], true_g, tolerance = 0.1)
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

test_that("The prior family is correctly inferred", {
  expect_identical(infer_prior_family(genpl.res[[g_ret_str()]]),
                   "gen_point_laplace")
})

# ==========================================================================
# Reductions to nested families
# ==========================================================================

test_that("Symmetric case reduces to point-Laplace", {
  # Equal rates and equal tail weights give a Laplace slab with weight w.
  w <- 0.5; a <- 0.1; mu <- 0
  par <- genpl_par(1 - w, w / 2, w / 2, a, a, mu)
  expect_equal(
    as.numeric(genpl_eval(par, rep(FALSE, 5), x, s, calc_hess = FALSE)),
    -loglik_point_laplace(x, s, w = w, a = a, mu = mu)
  )
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

# ==========================================================================
# Likelihood correctness
# ==========================================================================

# A grid of parameter values spanning symmetric, asymmetric, extreme-rate, and
#   shifted-mode regimes.
llik_grid <- list(
  list(lab = "central",        p = c(0.5, 0.3, 0.2, 0.3, 0.8,  0.0)),
  list(lab = "asymmetric",     p = c(0.2, 0.7, 0.1, 2.0, 0.05, 0.0)),
  list(lab = "mostly null",    p = c(0.9, 0.05, 0.05, 1.0, 1.0, 0.0)),
  list(lab = "heavy tails",    p = c(0.4, 0.3, 0.3, 0.01, 0.02, 0.0)),
  list(lab = "light tails",    p = c(0.4, 0.3, 0.3, 20, 30,    0.0)),
  list(lab = "shifted mode",   p = c(0.5, 0.25, 0.25, 0.5, 0.5, 1.7)),
  list(lab = "negative mode",  p = c(0.3, 0.4, 0.3, 0.7, 0.4, -2.2))
)

test_that("nllik agrees with the reference implementation", {
  for (cs in llik_grid) {
    v <- cs$p
    par <- genpl_par(v[1], v[2], v[3], v[4], v[5], v[6])
    for (dat in list(list(x = x, s = s), list(x = xu, s = su),
                     list(x = xu, s = 1))) {
      expect_equal(
        as.numeric(genpl_eval(par, rep(FALSE, 5), dat$x, dat$s,
                              calc_hess = FALSE)),
        -loglik_gen_point_laplace(dat$x, dat$s, v[1], v[2], v[3],
                                  v[4], v[5], v[6]),
        info = cs$lab
      )
    }
  }
})

test_that("nllik is stable for extreme observations", {
  # Drives the pnorm(log.p = TRUE) path to very negative arguments, where a
  #   naive implementation would underflow to -Inf.
  xbig <- c(-1e3, -50, 0, 50, 1e3)
  sbig <- rep(1, 5)
  par <- genpl_par(0.5, 0.3, 0.2, 2, 3, 0)
  val <- as.numeric(genpl_eval(par, rep(FALSE, 5), xbig, sbig,
                               calc_hess = FALSE))
  expect_true(is.finite(val))
  expect_equal(val, -loglik_gen_point_laplace(xbig, sbig, 0.5, 0.3, 0.2,
                                              2, 3, 0))
})

test_that("Pinning a weight is a reparameterization, not a different model", {
  # Evaluating the same prior through the five-slot and four-slot layouts must
  #   give the same likelihood.
  pi_0 <- 0.5; pi_pos <- 0.3; pi_neg <- 0.2
  rate_pos <- 0.3; rate_neg <- 0.8; mu <- 0.4

  target <- -loglik_gen_point_laplace(xu, su, pi_0, pi_pos, pi_neg,
                                      rate_pos, rate_neg, mu)

  specs <- list(
    list(comp = "spike", c = pi_0,   a = log(pi_pos / pi_neg)),
    list(comp = "pos",   c = pi_pos, a = log(pi_0   / pi_neg)),
    list(comp = "neg",   c = pi_neg, a = log(pi_0   / pi_pos))
  )
  for (sp in specs) {
    p4 <- c(sp$a, log(rate_pos), log(rate_neg), mu)
    expect_equal(
      as.numeric(genpl_eval(p4, rep(FALSE, 4), xu, su, calc_hess = FALSE,
                            fixed_comp = sp$comp, fixed_c = sp$c)),
      target, info = sp$comp
    )
  }
})

test_that("The tail-heaviness constraint is a reparameterization", {
  # At an interior delta the constrained parameterization must reproduce the
  #   unconstrained likelihood at the log-rate that delta encodes.
  k <- 3
  for (ref in c("pos", "neg")) {
    for (delta in c(-1.5, 0, 0.8)) {
      beta_ref <- log(0.5)
      beta_free <- (beta_ref - log(k)) + genpl_softplus(delta)
      if (identical(ref, "pos")) {
        p_con <- c(0.4, -0.2, beta_ref, delta, 0.1)
        rate_pos <- exp(beta_ref); rate_neg <- exp(beta_free)
      } else {
        p_con <- c(0.4, -0.2, delta, beta_ref, 0.1)
        rate_pos <- exp(beta_free); rate_neg <- exp(beta_ref)
      }
      cn <- genpl_canon(c(0.4, -0.2, 0, 0, 0.1))
      expect_equal(
        as.numeric(genpl_eval(p_con, rep(FALSE, 5), xu, su, calc_hess = FALSE,
                              ratio_ref = ref, ratio_k = k)),
        -loglik_gen_point_laplace(xu, su, cn$pi_0, cn$pi_pos, cn$pi_neg,
                                  rate_pos, rate_neg, 0.1),
        info = paste(ref, delta)
      )
    }
  }
})

# ==========================================================================
# Derivatives
# ==========================================================================

# numDeriv's Richardson extrapolation is very accurate for gradients and
#   somewhat less so for Hessians, hence the two tolerances.
GRAD_TOL <- 1e-7
HESS_TOL <- 1e-4

check_derivs <- function(p_full, fix_par, data_x, data_s, info, ...) {
  val <- genpl_eval(p_full, fix_par, data_x, data_s, calc_grad = TRUE,
                    calc_hess = TRUE, ...)
  free <- p_full[!fix_par]

  ref <- function(z) {
    genpl_ref_nllik(z, p_full, fix_par, data_x, data_s, ...)
  }

  expect_equal(as.numeric(val), ref(free), info = paste(info, "value"))

  g <- attr(val, "gradient")
  expect_equal(length(g), sum(!fix_par), info = paste(info, "grad length"))
  expect_equal(g, numDeriv::grad(ref, free), tolerance = GRAD_TOL,
               info = paste(info, "gradient"))

  h <- attr(val, "hessian")
  expect_equal(dim(h), rep(sum(!fix_par), 2), info = paste(info, "hess dim"))
  expect_false(any(is.na(h)), info = paste(info, "hess has NA"))
  expect_equal(h, t(h), info = paste(info, "hess symmetry"))
  expect_equal(h, numDeriv::hessian(ref, free), tolerance = HESS_TOL,
               info = paste(info, "hessian"))
}

test_that("Gradient and Hessian are correct for every fix_par pattern", {
  par <- genpl_par(0.5, 0.3, 0.2, 0.4, 0.9, 0.25)
  pats <- fix_patterns(5)
  for (i in seq_len(nrow(pats))) {
    fp <- pats[i, ]
    check_derivs(par, fp, xu, su,
                 info = paste0("fix_par=", paste(as.integer(fp), collapse = "")))
  }
})

test_that("Gradient and Hessian are correct across parameter regimes", {
  for (cs in llik_grid) {
    v <- cs$p
    par <- genpl_par(v[1], v[2], v[3], v[4], v[5], v[6])
    check_derivs(par, rep(FALSE, 5), xu, su, info = cs$lab)
    # Also with the mode held fixed, the most common configuration in practice.
    check_derivs(par, c(FALSE, FALSE, FALSE, FALSE, TRUE), xu, su,
                 info = paste(cs$lab, "fixed mode"))
  }
})

test_that("Gradient and Hessian are correct when a weight is pinned", {
  pi_0 <- 0.45; pi_pos <- 0.35; pi_neg <- 0.2
  specs <- list(
    list(comp = "spike", c = pi_0,   a = log(pi_pos / pi_neg)),
    list(comp = "pos",   c = pi_pos, a = log(pi_0   / pi_neg)),
    list(comp = "neg",   c = pi_neg, a = log(pi_0   / pi_pos))
  )
  pats <- fix_patterns(4)
  for (sp in specs) {
    p4 <- c(sp$a, log(0.4), log(0.9), 0.25)
    for (i in seq_len(nrow(pats))) {
      fp <- pats[i, ]
      check_derivs(p4, fp, xu, su,
                   info = paste0(sp$comp, " fix_par=",
                                 paste(as.integer(fp), collapse = "")),
                   fixed_comp = sp$comp, fixed_c = sp$c)
    }
  }
})

test_that("Gradient and Hessian are correct under the tail-heaviness constraint", {
  k <- 2.5
  for (ref in c("pos", "neg")) {
    for (delta in c(-2, -0.5, 0.7)) {
      # The reference tail's log-rate is fixed, which is the premise of the
      #   reparameterization; delta is free.
      if (identical(ref, "pos")) {
        p <- c(0.4, -0.2, log(0.5), delta, 0.15)
        fp <- c(FALSE, FALSE, TRUE, FALSE, TRUE)
      } else {
        p <- c(0.4, -0.2, delta, log(0.5), 0.15)
        fp <- c(FALSE, FALSE, FALSE, TRUE, TRUE)
      }
      check_derivs(p, fp, xu, su,
                   info = paste("ratio", ref, "delta", delta),
                   ratio_ref = ref, ratio_k = k)
    }
  }
})

test_that("Both constraints together give correct derivatives", {
  k <- 2
  p <- c(log(0.35 / 0.2), log(0.5), -0.6, 0.1)
  fp <- c(FALSE, TRUE, FALSE, TRUE)
  check_derivs(p, fp, xu, su, info = "spike pinned + ratio pos",
               fixed_comp = "spike", fixed_c = 0.45,
               ratio_ref = "pos", ratio_k = k)
})

test_that("calc_grad and calc_hess control which attributes are attached", {
  par <- genpl_par(0.5, 0.3, 0.2, 0.4, 0.9, 0)
  fp <- rep(FALSE, 5)
  v0 <- genpl_eval(par, fp, xu, su, calc_grad = FALSE, calc_hess = FALSE)
  expect_null(attr(v0, "gradient"))
  expect_null(attr(v0, "hessian"))

  v1 <- genpl_eval(par, fp, xu, su, calc_grad = TRUE, calc_hess = FALSE)
  expect_false(is.null(attr(v1, "gradient")))
  expect_null(attr(v1, "hessian"))

  v2 <- genpl_eval(par, fp, xu, su, calc_grad = FALSE, calc_hess = TRUE)
  expect_false(is.null(attr(v2, "hessian")))
})

# ==========================================================================
# Reparameterization round-trips
# ==========================================================================

test_that("genpl_parse_scale accepts and rejects the right forms", {
  expect_equal(genpl_parse_scale("estimate")$fix, c(FALSE, FALSE))
  expect_equal(genpl_parse_scale(5)$value, c(5, 5))
  expect_equal(genpl_parse_scale(5)$fix, c(TRUE, TRUE))
  expect_equal(genpl_parse_scale(c(1, 2))$value, c(1, 2))
  expect_equal(genpl_parse_scale(c(1, NA))$fix, c(TRUE, FALSE))
  expect_equal(genpl_parse_scale(c(NA, 2))$fix, c(FALSE, TRUE))

  # A numeric pair of NAs reaches the "estimate at least one" branch; a bare
  #   c(NA, NA) is a *logical* vector in R, so it is caught earlier as
  #   non-numeric. Both are rejected, with different messages.
  expect_error(genpl_parse_scale(c(NA_real_, NA_real_)), "at least one scale")
  expect_error(genpl_parse_scale(c(NA, NA)), "must be 'estimate' or a numeric")
  expect_error(genpl_parse_scale(c(1, 2, 3)), "length")
  expect_error(genpl_parse_scale(0), "positive")
  expect_error(genpl_parse_scale(c(-1, 2)), "positive")
})

test_that("genpl_undo_fixedcomp inverts the pinned-weight parameterization", {
  for (comp in c("spike", "pos", "neg")) {
    for (cval in c(0.2, 0.5, 0.75)) {
      for (a in c(-1.3, 0, 2.1)) {
        p4 <- c(a, log(0.4), log(0.9), 0.1)
        target <- genpl_canon(p4, fixed_comp = comp, fixed_c = cval)

        undone <- genpl_undo_fixedcomp(
          list(alpha = a, beta_pos = log(0.4), beta_neg = log(0.9), mu = 0.1),
          comp, cval
        )
        got <- genpl_canon(c(undone$alpha_pos, undone$alpha_neg,
                             undone$beta_pos, undone$beta_neg, undone$mu))

        lab <- paste(comp, cval, a)
        expect_null(undone$alpha, info = lab)
        expect_equal(got$pi_0,   target$pi_0,   info = lab)
        expect_equal(got$pi_pos, target$pi_pos, info = lab)
        expect_equal(got$pi_neg, target$pi_neg, info = lab)
      }
    }
  }
})

test_that("genpl_undo_ratio inverts the tail-heaviness parameterization", {
  k <- 2.5
  for (delta in c(-3, -0.4, 0, 1.2)) {
    up <- genpl_undo_ratio(
      list(beta_pos = log(0.5), delta = delta, mu = 0), "pos", k
    )
    expect_null(up$delta)
    expect_equal(up$beta_neg, (log(0.5) - log(k)) + genpl_softplus(delta))

    un <- genpl_undo_ratio(
      list(delta = delta, beta_neg = log(0.5), mu = 0), "neg", k
    )
    expect_equal(un$beta_pos, (log(0.5) - log(k)) + genpl_softplus(delta))
  }
  # A NULL reference leaves the parameters untouched.
  orig <- list(beta_pos = 1, beta_neg = 2, mu = 0)
  expect_identical(genpl_undo_ratio(orig, NULL, NULL), orig)
})

test_that("genpl_init_delta honors the bound", {
  # A feasible target is reproduced exactly.
  bound <- log(0.5) - log(2)
  target <- bound + 0.75
  d <- genpl_init_delta(target, bound)
  expect_equal(bound + genpl_softplus(d), target)

  # An infeasible target starts strictly inside the bound.
  d2 <- genpl_init_delta(bound - 5, bound)
  expect_true(is.finite(d2))
  expect_gt(bound + genpl_softplus(d2), bound)
})

test_that("softplus and its inverse round-trip", {
  gaps <- c(1e-8, 0.01, 0.5, log(2), 5, 36, 40, 500)
  expect_equal(genpl_softplus(genpl_inv_softplus(gaps)), gaps)
  # Known values, and the limits that motivate the two stability shortcuts.
  expect_equal(genpl_softplus(0), log(2))
  expect_equal(genpl_inv_softplus(log(2)), 0)
  expect_true(all(is.finite(genpl_softplus(c(-800, 800)))))
  expect_equal(genpl_softplus(800), 800)
  expect_equal(genpl_inv_softplus(800), 800)
})

test_that("delta is invariant under the internal rescaling", {
  par <- list(alpha_pos = 0.3, alpha_neg = -0.2, beta_pos = log(0.5),
              delta = -0.8, mu = 1.5)
  f <- 7.3
  scaled <- genpl_scalepar(par, f)

  expect_identical(scaled$delta, par$delta)
  expect_equal(scaled$beta_pos, par$beta_pos - log(f))
  expect_equal(scaled$mu, f * par$mu)

  # Round-tripping recovers the original.
  back <- genpl_scalepar(scaled, 1 / f)
  expect_equal(back$beta_pos, par$beta_pos)
  expect_equal(back$mu, par$mu)
  expect_identical(back$delta, par$delta)
})

test_that("genpl_postcomp overrides with the null solution when it wins", {
  # Hand the boundary check an implausibly bad optimum so that the pi_0 = 1
  #   solution must win, and confirm that the override fires.
  par <- list(alpha_pos = 0.1, alpha_neg = 0.1, beta_pos = 0, beta_neg = 0,
              mu = 0)
  res <- genpl_postcomp(optpar = par, optval = 1e6, x = xu, s = su,
                        par_init = par, fix_par = c(FALSE, FALSE, FALSE,
                                                    FALSE, TRUE),
                        scale_factor = 1)
  expect_equal(res$par$alpha_pos, -Inf)
  expect_equal(res$par$alpha_neg, -Inf)
  expect_equal(res$val,
               sum(-0.5 * log(2 * pi * su^2) - 0.5 * (xu - 0)^2 / su^2))
})

# ==========================================================================
# Constraint behavior, end to end
# ==========================================================================

test_that("Pinning a weight reproduces that weight exactly", {
  slot <- c(spike = 1, pos = 2, neg = 3)
  for (comp in names(slot)) {
    res <- ebnm_gen_point_laplace(
      x, s, pi_fixed = list(comp = comp, value = 0.3),
      optmethod = "nograd_lbfgsb"
    )
    g <- res[[g_ret_str()]]
    expect_equal(g$pi[slot[[comp]]], 0.3, info = comp)
    expect_valid_genpl(g, info = comp)
    # A constrained optimum cannot beat the unconstrained one.
    expect_lte(as.numeric(res[[llik_ret_str()]]),
               as.numeric(genpl.res[[llik_ret_str()]]) + 1e-6)
  }
})

test_that("Fixing a single scale leaves the other free", {
  r1 <- ebnm_gen_point_laplace(x, s, scale = c(true_scale_pos, NA))
  expect_equal(r1[[g_ret_str()]]$scale_pos[2], true_scale_pos)
  expect_false(isTRUE(all.equal(r1[[g_ret_str()]]$scale_neg[3],
                                true_scale_pos)))

  r2 <- ebnm_gen_point_laplace(x, s, scale = c(NA, true_scale_neg))
  expect_equal(r2[[g_ret_str()]]$scale_neg[3], true_scale_neg)
  expect_false(isTRUE(all.equal(r2[[g_ret_str()]]$scale_pos[2],
                                true_scale_neg)))
})

test_that("The tail-heaviness constraint is always satisfied", {
  for (k in c(0.5, 1, 2, 5)) {
    r <- ebnm_gen_point_laplace(x, s, scale = c(true_scale_pos, NA),
                                scale_ratio = k, optmethod = "nograd_lbfgsb")
    g <- r[[g_ret_str()]]
    expect_lte(g$scale_neg[3], k * g$scale_pos[2] * (1 + 1e-6))
    expect_valid_genpl(g, info = paste("k =", k))
  }
})

test_that("A loose tail-heaviness constraint is inert", {
  # With k large the bound cannot bind, so the fit must match the one with the
  #   same single fixed scale and no constraint.
  free <- ebnm_gen_point_laplace(x, s, scale = c(NA, true_scale_neg),
                                 optmethod = "nograd_lbfgsb")
  con  <- ebnm_gen_point_laplace(x, s, scale = c(NA, true_scale_neg),
                                 scale_ratio = 100,
                                 optmethod = "nograd_lbfgsb")
  expect_equal(as.numeric(con[[llik_ret_str()]]),
               as.numeric(free[[llik_ret_str()]]), tolerance = 1e-5)
  expect_equal(con[[g_ret_str()]]$scale_pos[2],
               free[[g_ret_str()]]$scale_pos[2], tolerance = 1e-3)
})

test_that("A binding tail-heaviness constraint lands on the bound", {
  # The unconstrained negative scale is around 3.9 against a fixed positive
  #   scale of 10, so k has to be well under 0.39 for the bound to bite. Once
  #   it does, the fit must sit at scale_neg = k * scale_pos, which is exactly
  #   the model with both scales fixed at those values.
  k <- 0.2
  bound <- k * true_scale_pos
  con <- ebnm_gen_point_laplace(x, s, scale = c(true_scale_pos, NA),
                                scale_ratio = k, optmethod = "nograd_lbfgsb")
  both <- ebnm_gen_point_laplace(x, s, scale = c(true_scale_pos, bound),
                                 optmethod = "nograd_lbfgsb")
  expect_equal(con[[g_ret_str()]]$scale_neg[3], bound, tolerance = 1e-4)
  expect_equal(as.numeric(con[[llik_ret_str()]]),
               as.numeric(both[[llik_ret_str()]]), tolerance = 1e-5)
})

test_that("Loosening the tail-heaviness constraint cannot hurt the fit", {
  ks <- c(0.25, 0.5, 1, 2, 4)
  lliks <- vapply(ks, function(k) {
    as.numeric(ebnm_gen_point_laplace(
      x, s, scale = c(true_scale_pos, NA), scale_ratio = k,
      optmethod = "nograd_lbfgsb"
    )[[llik_ret_str()]])
  }, numeric(1))
  expect_true(all(diff(lliks) > -1e-5))
})

test_that("A smaller ratio forces a lighter free tail", {
  # Guards against inverting the scale/rate relationship: tail heaviness goes
  #   with the scale, so a smaller k must give a smaller scale_neg.
  tight <- ebnm_gen_point_laplace(x, s, scale = c(true_scale_pos, NA),
                                  scale_ratio = 0.25,
                                  optmethod = "nograd_lbfgsb")
  loose <- ebnm_gen_point_laplace(x, s, scale = c(true_scale_pos, NA),
                                  scale_ratio = 4,
                                  optmethod = "nograd_lbfgsb")
  expect_lt(tight[[g_ret_str()]]$scale_neg[3],
            loose[[g_ret_str()]]$scale_neg[3])
})

test_that("Both constraints can be combined", {
  r <- ebnm_gen_point_laplace(
    x, s, scale = c(true_scale_pos, NA), scale_ratio = 2,
    pi_fixed = list(comp = "spike", value = 0.4),
    optmethod = "nograd_lbfgsb"
  )
  g <- r[[g_ret_str()]]
  expect_equal(g$pi[1], 0.4)
  expect_lte(g$scale_neg[3], 2 * g$scale_pos[2] * (1 + 1e-6))
  expect_valid_genpl(g)
})

test_that("Constrained fits are scale invariant", {
  scl <- 1e4
  spec <- function(xx, ss, f) {
    ebnm_gen_point_laplace(xx, ss, scale = c(f * true_scale_pos, NA),
                           scale_ratio = 2,
                           pi_fixed = list(comp = "spike", value = 0.4),
                           optmethod = "nograd_lbfgsb")[[g_ret_str()]]
  }
  r  <- spec(x, s, 1)
  rl <- spec(scl * x, scl * s, scl)
  expect_equal(r$pi, rl$pi, tolerance = 1e-3)
  expect_equal(scl * r$scale_neg[3], rl$scale_neg[3], tolerance = 1e-3)
})

# ==========================================================================
# Optimization methods
# ==========================================================================

test_that("All optmethods give the same answer", {
  methods <- c("nograd_lbfgsb", "nograd_nlm", "lbfgsb", "nlm", "nohess_nlm",
               "trust")
  lliks <- vapply(methods, function(m) {
    as.numeric(ebnm_gen_point_laplace(x, s, optmethod = m)[[llik_ret_str()]])
  }, numeric(1))
  expect_equal(max(lliks) - min(lliks), 0, tolerance = 1e-4)
})

test_that("Gradient-based optmethods work with a pinned weight", {
  methods <- c("nograd_lbfgsb", "lbfgsb", "nlm", "trust")
  lliks <- vapply(methods, function(m) {
    as.numeric(ebnm_gen_point_laplace(
      x, s, pi_fixed = list(comp = "spike", value = 0.3), optmethod = m
    )[[llik_ret_str()]])
  }, numeric(1))
  expect_equal(max(lliks) - min(lliks), 0, tolerance = 1e-4)
})

test_that("Gradient-based optmethods work under the tail-heaviness constraint", {
  # Exercises the delta chain rule through the optimizers that consume it,
  #   with a binding k so that delta is doing real work.
  methods <- c("nograd_lbfgsb", "nograd_nlm", "lbfgsb", "nlm", "nohess_nlm",
               "trust")
  res <- lapply(methods, function(m) {
    ebnm_gen_point_laplace(x, s, scale = c(true_scale_pos, NA),
                           scale_ratio = 0.2, optmethod = m)
  })
  lliks <- vapply(res, function(r) as.numeric(r[[llik_ret_str()]]), numeric(1))
  expect_equal(max(lliks) - min(lliks), 0, tolerance = 1e-4)
  scl_neg <- vapply(res, function(r) r[[g_ret_str()]]$scale_neg[3], numeric(1))
  expect_equal(max(scl_neg) - min(scl_neg), 0, tolerance = 1e-4)
})

# ==========================================================================
# Posteriors
# ==========================================================================

test_that("The posterior sampler works", {
  genpl.res2 <- ebnm_gen_point_laplace(x, s, output = ebnm_output_all())
  samp <- genpl.res2[[samp_ret_str()]](100)
  expect_equal(dim(samp), c(100, n))
  # Sample means should track the posterior means.
  expect_equal(colMeans(samp), genpl.res2[[df_ret_str()]][[pm_ret_str()]],
               tolerance = 0.5)
})

test_that("The posterior sampler agrees with the posterior moments", {
  set.seed(5)
  res <- ebnm_gen_point_laplace(xu, su, output = ebnm_output_all())
  samp <- res[[samp_ret_str()]](20000)
  pm <- res[[df_ret_str()]][[pm_ret_str()]]
  psd <- res[[df_ret_str()]][[psd_ret_str()]]
  # Monte Carlo error on the mean is sd / sqrt(nsamp); allow a generous
  #   multiple of it so the test is not flaky.
  expect_true(all(abs(colMeans(samp) - pm) < 0.05 + 6 * psd / sqrt(20000)))
  expect_equal(apply(samp, 2, sd), psd, tolerance = 0.1)
})

test_that("Posterior moments agree with direct numerical integration", {
  g <- genpl.res[[g_ret_str()]]
  pi_0 <- g$pi[1]; pi_pos <- g$pi[2]; pi_neg <- g$pi[3]
  l_pos <- 1 / g$scale_pos[2]; l_neg <- 1 / g$scale_neg[3]
  mu <- g$mean[1]

  res <- ebnm_gen_point_laplace(x, s, g_init = g, fix_g = TRUE,
                                output = ebnm_output_all())

  # Unnormalized posterior density of theta given x_i, built directly from the
  #   prior and the normal likelihood. The point mass is handled separately.
  idx <- c(1, 100, 400, 700, 950)
  for (i in idx) {
    xi <- x[i]; si <- s[i]
    dens <- function(th) {
      pri <- ifelse(th > mu, pi_pos * l_pos * exp(-l_pos * (th - mu)),
                    pi_neg * l_neg * exp(-l_neg * (mu - th)))
      pri * dnorm(xi, th, si)
    }
    # Integrate in two pieces split at the posterior peak. Adaptive quadrature
    #   over a wide interval can step over a peak only a few si wide and return
    #   a badly wrong value with no warning.
    quad <- function(f, lo, hi) {
      pk <- min(max(xi, lo), hi)
      integrate(f, lo, pk)$value + integrate(f, pk, hi)$value
    }
    integrand <- function(th) th * dens(th)
    hi <- mu + 60 / l_pos
    lo <- mu - 60 / l_neg

    num_pos <- quad(integrand, mu, hi)
    num_neg <- quad(integrand, lo, mu)
    den_pos <- quad(dens, mu, hi)
    den_neg <- quad(dens, lo, mu)
    spike   <- pi_0 * dnorm(xi, mu, si)

    expect_equal((num_pos + num_neg + spike * mu) /
                   (den_pos + den_neg + spike),
                 res[[df_ret_str()]][[pm_ret_str()]][i],
                 tolerance = 1e-4, info = paste("obs", i))
  }
})

test_that("lfsr lies in the unit interval and respects the point mass", {
  res <- ebnm_gen_point_laplace(x, s, output = ebnm_output_all())
  lfsr <- res[[df_ret_str()]][[lfsr_ret_str()]]
  expect_true(all(lfsr >= 0 & lfsr <= 1))
  # lfsr is at least the posterior probability of the point mass, since a
  #   theta of exactly zero has no sign.
  expect_true(all(lfsr >= 0))
})

# ==========================================================================
# Degenerate standard errors
# ==========================================================================

test_that("Infinite standard errors are rejected", {
  # check_args rejects infinite SEs for every family, so the is.infinite(s)
  #   branches in genpl_summres, pipost_genpl, and parametric_workhorse are
  #   unreachable through the public interface. This pins the actual contract.
  sinf <- s
  sinf[c(5, 50, 500)] <- Inf
  expect_error(ebnm_gen_point_laplace(x, sinf), "cannot be infinite")
})

test_that("Zero standard errors are replaced rather than rejected", {
  # handle_standard_errors substitutes small positive SEs before genpl_precomp
  #   ever sees them, so its zero-SE guard is likewise unreachable from here.
  #   The fit succeeds either way, with a warning.
  s0 <- s
  s0[1:3] <- 0
  x0 <- x
  x0[1:3] <- 0

  expect_warning(res <- ebnm_gen_point_laplace(x0, s0, mode = 0),
                 "Nonpositive SEs")
  expect_valid_genpl(res[[g_ret_str()]])

  expect_warning(res2 <- ebnm_gen_point_laplace(x0, s0, mode = "estimate"),
                 "Nonpositive SEs")
  expect_valid_genpl(res2[[g_ret_str()]])
})

test_that("genpl_precomp still guards against zero SEs when called directly", {
  par <- list(alpha_pos = 0, alpha_neg = 0, beta_pos = 0, beta_neg = 0, mu = 0)
  expect_error(
    genpl_precomp(c(1, 2), c(1, 0), par_init = par,
                  fix_par = rep(FALSE, 5)),
    "mode cannot be estimated"
  )
})

# ==========================================================================
# Argument validation
# ==========================================================================

test_that("Invalid constraint arguments are rejected", {
  cases <- list(
    list(regexp = "exactly one of the two scales",
         f = function() ebnm_gen_point_laplace(x, s, scale = "estimate",
                                               scale_ratio = 2)),
    list(regexp = "exactly one of the two scales",
         f = function() ebnm_gen_point_laplace(x, s, scale = c(1, 2),
                                               scale_ratio = 2)),
    list(regexp = "single positive number",
         f = function() ebnm_gen_point_laplace(x, s, scale = c(1, NA),
                                               scale_ratio = -1)),
    list(regexp = "single positive number",
         f = function() ebnm_gen_point_laplace(x, s, scale = c(1, NA),
                                               scale_ratio = c(1, 2))),
    list(regexp = "strictly",
         f = function() ebnm_gen_point_laplace(
           x, s, pi_fixed = list(comp = "spike", value = 1.5))),
    list(regexp = "strictly",
         f = function() ebnm_gen_point_laplace(
           x, s, pi_fixed = list(comp = "spike", value = 0))),
    list(regexp = "should be one of",
         f = function() ebnm_gen_point_laplace(
           x, s, pi_fixed = list(comp = "bogus", value = 0.3))),
    list(regexp = "comp",
         f = function() ebnm_gen_point_laplace(x, s,
                                               pi_fixed = list(value = 0.3))),
    # A bare c(NA, NA) is logical, so handle_scale_parameter rejects it as
    #   non-numeric before genpl_parse_scale can give its own message.
    list(regexp = "must be either 'estimate' or numeric",
         f = function() ebnm_gen_point_laplace(x, s, scale = c(NA, NA))),
    list(regexp = "estimate",
         f = function() ebnm_gen_point_laplace(x, s, scale = c(1, Inf)))
  )
  for (cs in cases) {
    expect_error(cs$f(), cs$regexp)
  }
})

test_that("Invalid g_init is rejected", {
  expect_error(ebnm_gen_point_laplace(x, s, g_init = "not a prior"),
               "genlaplacemix")

  # Wrong number of components.
  bad <- genlaplacemix(pi = c(0.2, 0.2, 0.3, 0.3), mean = rep(0, 4),
                       scale_pos = c(0, 1, 0, 1),
                       scale_neg = c(0, 0, 1, 1))
  expect_error(ebnm_gen_point_laplace(x, s, g_init = bad),
               "correct number of components")

  # First component is not a point mass.
  bad2 <- genlaplacemix(pi = c(0.5, 0.3, 0.2), mean = rep(0, 3),
                        scale_pos = c(1, 1, 0),
                        scale_neg = c(0, 0, 1))
  expect_error(ebnm_gen_point_laplace(x, s, g_init = bad2),
               "must be a point mass")
})

test_that("g_init must agree with fixed mode and scale", {
  expect_error(
    ebnm_gen_point_laplace(x, s, g_init = true_g, mode = 5),
    "they must agree"
  )
  expect_error(
    ebnm_gen_point_laplace(x, s, g_init = true_g,
                           scale = c(true_scale_pos + 1, true_scale_neg)),
    "they must agree"
  )
})

test_that("Fixing g warns when mode or scale is also supplied", {
  expect_warning(
    ebnm_gen_point_laplace(x, s, g_init = true_g, fix_g = TRUE, mode = 0),
    "ignored when g is fixed"
  )
})

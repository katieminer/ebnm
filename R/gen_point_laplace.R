#' Constructor for genlaplacemix class
#'
#' Creates a finite mixture of Generalized Laplace distributions.
#'
#' @param pi A vector of mixture proportions.
#'
#' @param mean A vector of means.
#'
#' @param scale_pos A vector of positive scale parameters.
#' 
#' @param scale_neg A vector of negative scale parameters.
#'
#' @return An object of class \code{genlaplacemix} (a list with elements
#'   \code{pi}, \code{mean}, \code{scale_pos}, and \code{scale_neg}, described above).
#'
#' @export
#'
genlaplacemix <- function(pi, mean, scale_pos, scale_neg) {
  structure(data.frame(pi, mean, scale_pos, scale_neg), class="genlaplacemix")
}

#' @importFrom stats pexp
#' @importFrom ashr comp_cdf
#'
#' @method comp_cdf genlaplacemix
#'
#' @export
#'
comp_cdf.genlaplacemix = function (m, y, lower.tail = TRUE) {
  
  # Internal function to evaluate the CDF of the generalized components
  pgenlaplace <- function(q, location = 0, scale_pos = 1, scale_neg = 1, lower.tail = TRUE) {
    q <- q - location
    res <- numeric(length(q))
    
    # Identify the type of each component in the vectorized inputs
    pos_idx <- (scale_pos > 0 & scale_neg == 0)
    neg_idx <- (scale_pos == 0 & scale_neg > 0)
    pt_idx  <- (scale_pos == 0 & scale_neg == 0)
    
    if (lower.tail) {
      # Lower tail CDF: P(X <= q)
      #want to add components from each tail
      if (any(pos_idx)) { #if q < 0, pos slab has no weight. else, we use CDF of exp
        res[pos_idx] <- ifelse(q[pos_idx] < 0, 0, pexp(q[pos_idx], 1 / scale_pos[pos_idx])) #rate
      }
      if (any(neg_idx)) { #if q > 0, we alr have all of negative slab, else, do cdf
        res[neg_idx] <- ifelse(q[neg_idx] > 0, 1, pexp(abs(q[neg_idx]), 1 / scale_neg[neg_idx], lower.tail = FALSE))
      }
      if (any(pt_idx)) { #if q >= 0, we have all of spike, else we have no weight here
        res[pt_idx]  <- ifelse(q[pt_idx] >= 0, 1, 0)
      }
    } else { #Calculate this way to avoid errors
      # Upper tail CDF: P(X > q)
      #want to add components from each tail
      if (any(pos_idx)) { #if q < 0, pos slab is included entirely. else, we use CDF of exp, with lower_tail FALSE since we are interested in Up. Tail
        res[pos_idx] <- ifelse(q[pos_idx] < 0, 1, pexp(q[pos_idx], 1 / scale_pos[pos_idx], lower.tail = FALSE))
      }
      if (any(neg_idx)) { #if q > 0, neg slab has no weight. else, we use CDF of exp, with lower_tail
        res[neg_idx] <- ifelse(q[neg_idx] > 0, 0, pexp(abs(q[neg_idx]), 1 / scale_neg[neg_idx]))
      }
      if (any(pt_idx)) { #if we are looking at on neg slab or spike, spike is included in cdf
        res[pt_idx]  <- ifelse(q[pt_idx] <= 0, 1, 0)
      }
    }
    
    return(res)
  }
  
  return(vapply(y, pgenlaplace, m$mean, m$mean, m$scale_pos, m$scale_neg, lower.tail))
}

# The gen-point-Laplace family uses the above ebnm class genlaplacemix.
# We cannot use check_g_init like all other priors, as check_g_init assumes a 2-component prior with a single scale column
# this family has 3 components (point mass + two tails) and separate scale_pos/scale_neg columns, so we duplicate the checks here.

genpl_checkg <- function(g_init, fix_g, mode, scale, pointmass, call) {
  if (is.null(g_init)) {
    return(invisible(NULL))
  }
  
  if (!inherits(g_init, "genlaplacemix")) {
    stop("g_init must be NULL or an object of class genlaplacemix.")
  }
  
  # Components are (point mass, positive slab, negative slab); 
  # the point mass is omitted when its weight is zero.
  ncomp <- length(g_init$pi)
  if (!(ncomp == 2 || (pointmass && ncomp == 3))) {
    stop("g_init does not have the correct number of components.")
  }
  if (ncomp == 3 && !(g_init$scale_pos[1] == 0 && g_init$scale_neg[1] == 0)) { #gives a non-zero scale to the point mass from either pos/neg 
    stop("The first component of g_init must be a point mass.")
  }
  
  if (fix_g && (!is.null(call$mode) || !is.null(call$scale))) { #ignore instructions about mode or scale because we have them already decided in fixing g
    warning("mode and scale parameters are ignored when g is fixed.")
  }
  
  if (!fix_g) {
    if (!is.null(call$mode)
        && !identical(mode, "estimate")
        && !isTRUE(all.equal(g_init$mean, rep(mode, length(g_init$mean))))) {
      stop("If mode is fixed and g_init is supplied, they must agree.")
    }
    if (!is.null(call$scale) && !identical(scale, "estimate")) {
      parsed  <- genpl_parse_scale(scale)
      g_scale <- c(g_init$scale_pos[ncomp - 1], g_init$scale_neg[ncomp])
      # Only the scales that are actually fixed have to agree; the other is
      # free to be estimated away from its starting value.
      fixed <- parsed$fix
      if (!isTRUE(all.equal(g_scale[fixed], parsed$value[fixed]))) {
        stop("If scale is fixed and g_init is supplied, they must agree.")
      }
    }
  }
  
  return(invisible(NULL))
}

#  Parameters are alpha_pos and alpha_neg (softmax of pi weights), beta_pos, beta_neg, and mu. 
#  Note that lambda is the reciprocal of the scale: beta = log(lambda) = -log(scale).
#
#' @importFrom stats var
#'
genpl_initpar <- function(g_init, mode, scale, pointmass, x, s,
                          fixed_comp = "none", fixed_c = NULL,
                          ratio_ref = NULL, ratio_k = NULL) {
  eps <- 1e-8

  # Recover starting weights, scales, and mode in canonical (untransformed)
  # form, either from g_init or from the data.
  if (!is.null(g_init)) {
    # Components are ordered (point mass, positive slab, negative slab);
    # the point mass is omitted when its weight is zero.
    if (length(g_init$pi) == 3) {
      pi_vec <- g_init$pi
      if (any(pi_vec <= 0)) {
        pi_vec[pi_vec <= 0] <- eps
        pi_vec <- pi_vec / sum(pi_vec)
      }
      scale_pos <- g_init$scale_pos[2]
      scale_neg <- g_init$scale_neg[3]
    }
    else {
      pi_vec    <- c(eps, g_init$pi) / sum(c(eps, g_init$pi)) #add in incredibly small pi_0 weight with no spike
      scale_pos <- g_init$scale_pos[1]
      scale_neg <- g_init$scale_neg[2]
    }
    mu <- g_init$mean[1]
  } else {
    if (pointmass) {
      pi_vec <- c(0.50, 0.25, 0.25)
    } else {
      pi_vec <- c(eps, 0.50, 0.50)
    }

    # Moment match on the variance in excess of the noise, used for whichever
    # scales are not supplied (variance of exp = scale^2).
    parsed     <- genpl_parse_scale(scale)
    excess_var <- max(var(x) - mean(s^2), 1e-4)
    default    <- sqrt(excess_var)
    scale_pos  <- if (is.na(parsed$value[1])) default else parsed$value[1]
    scale_neg  <- if (is.na(parsed$value[2])) default else parsed$value[2]

    mu <- if (identical(mode, "estimate")) mean(x) else mode
  }

  # Weights. Pinning one weight leaves the other two sharing a single free
  # log-ratio; otherwise both logits are free, relative to the spike. The
  # ordering of the free pair matches the f_1 / f_2 split in genpl_nllik.
  par <- list()
  if (identical(fixed_comp, "none")) {
    par$alpha_pos <- log(pi_vec[2] / pi_vec[1])
    par$alpha_neg <- log(pi_vec[3] / pi_vec[1])
  } else {
    free <- switch(fixed_comp,
                   spike = c(2, 3),   # pi_pos over pi_neg
                   pos   = c(1, 3),   # pi_0   over pi_neg
                   neg   = c(1, 2))   # pi_0   over pi_pos
    par$alpha <- log(pi_vec[free[1]] / pi_vec[free[2]])
  }

  # Rates, as beta = log(lambda) = -log(scale). Under the tail-heaviness
  # constraint the free tail's slot instead carries delta, and is named so that
  # genpl_scalepar can tell the two apart. Slot order has to stay
  # (alpha(s), positive, negative, mu) to match fix_par.
  beta_pos <- -log(scale_pos)
  beta_neg <- -log(scale_neg)

  if (is.null(ratio_ref)) {
    par$beta_pos <- beta_pos
    par$beta_neg <- beta_neg
  } else if (identical(ratio_ref, "pos")) {
    par$beta_pos <- beta_pos
    par$delta    <- genpl_init_delta(beta_neg, beta_pos - log(ratio_k))
  } else {
    par$delta    <- genpl_init_delta(beta_pos, beta_neg - log(ratio_k))
    par$beta_neg <- beta_neg
  }

  par$mu <- mu

  return(par)
}

# log(1 + exp(d)), mapping delta onto the gap above the bound. Not exp(d),
#   which would make the rate exp(B + exp(d)) and overflow by d = 6.
#
genpl_softplus <- function(d) {
  # log1p(exp(d)) overflows for large d; this form evaluates exp at -|d|.
  return(pmax(d, 0) + log1p(exp(-abs(d))))
}

# log(exp(gap) - 1), the inverse, for gap > 0.
#
genpl_inv_softplus <- function(gap) {
  # expm1 overflows above ~709, and past 37 the result equals gap to machine
  # precision anyway.
  return(ifelse(gap > 37, gap, log(expm1(gap))))
}

# Starting delta, where beta = bound + softplus(delta): reproduce the target
#   log-rate if it is feasible, else start just inside the bound.
#
genpl_init_delta <- function(beta_target, bound) {
  gap <- beta_target - bound
  if (is.finite(gap) && gap > 0) {
    return(genpl_inv_softplus(gap))
  }
  return(genpl_inv_softplus(0.1))
}

#Scaling function for numerical stability
genpl_scalepar <- function(par, scale_factor) {
  # Adjust the positive log-rate
  if (!is.null(par$beta_pos)) {
    par$beta_pos <- par$beta_pos - log(scale_factor)
  }
  # Adjust the negative log-rate
  if (!is.null(par$beta_neg)) {
    par$beta_neg <- par$beta_neg - log(scale_factor)
  }

  # A delta slot is deliberately left alone. It measures the gap between the
  # constrained log-rate and its bound, and since both shift by
  # -log(scale_factor) the gap is invariant.

  # Adjust the mode
  if (!is.null(par$mu)) {
    par$mu <- scale_factor * par$mu
  }
  return(par)
}
    
# No precomputations are done for generalized point-Laplace, but this is what
#   carries the constraint specification through to genpl_nllik: whatever is
#   returned here is merged into the parameters that mle_parametric passes to
#   the negative log likelihood on every call.
genpl_precomp <- function(x, s, par_init, fix_par, fixed_comp = "none",
                          fixed_c = NULL, ratio_ref = NULL, ratio_k = NULL) {
  # Check if the mode (mu) is fixed and any standard error is 0 as this prevents dl/dmu
  # Parameter vector is ordered
  # c(alpha_pos, alpha_neg, beta_pos, beta_neg, mu), or
  # c(alpha, beta_pos, beta_neg, mu) when a weight is pinned, so mu is last
  # either way.

  fix_mu  <- fix_par[length(fix_par)]
  if (!fix_mu && any(s == 0)) {
    stop("The mode cannot be estimated if any SE is zero (the gradient does ",
         "not exist).")
  }
  return(list(fixed_comp = fixed_comp,
              fixed_c    = fixed_c,
              ratio_ref  = ratio_ref,
              ratio_k    = ratio_k))
}
  
# genpl_nllik is defined in genpl_likelihood.R.
# It evaluates the negative log-likelihood, gradient, and Hessian for the
# generalized point-Laplace prior, and additionally supports a fixed-component
# parameterization (fixed_comp = "spike"/"pos"/"neg") for constrained estimation.
# parametric_workhorse finds it by name from the shared package namespace.

# Undo the pinned-weight reparameterization, recovering the canonical pair of
#   softmax logits from the single free log-ratio and the pinned weight. The
#   weights below are built exactly as in genpl_nllik.
#
genpl_undo_fixedcomp <- function(par, fixed_comp, fixed_c) {
  if (identical(fixed_comp, "none")) {
    return(par)
  }

  exp_a <- exp(par$alpha)
  f_1   <- exp_a / (1 + exp_a)
  f_2   <- 1 / (1 + exp_a)
  rem_c <- 1 - fixed_c

  if (identical(fixed_comp, "spike")) {
    pi_0   <- fixed_c;      pi_pos <- f_1 * rem_c; pi_neg <- f_2 * rem_c
  } else if (identical(fixed_comp, "pos")) {
    pi_pos <- fixed_c;      pi_0   <- f_1 * rem_c; pi_neg <- f_2 * rem_c
  } else {
    pi_neg <- fixed_c;      pi_0   <- f_1 * rem_c; pi_pos <- f_2 * rem_c
  }

  par$alpha     <- NULL
  par$alpha_pos <- log(pi_pos / pi_0)
  par$alpha_neg <- log(pi_neg / pi_0)

  return(par)
}

# Undo the tail-heaviness reparameterization, replacing the delta slot with the
#   log-rate that it encodes.
#
genpl_undo_ratio <- function(par, ratio_ref, ratio_k) {
  if (is.null(ratio_ref)) {
    return(par)
  }

  if (identical(ratio_ref, "pos")) {
    par$beta_neg <- (par$beta_pos - log(ratio_k)) + genpl_softplus(par$delta)
  } else {
    par$beta_pos <- (par$beta_neg - log(ratio_k)) + genpl_softplus(par$delta)
  }
  par$delta <- NULL

  return(par)
}

# Postcomputations: restore the canonical parameterization, then check the
#   boundary solution of solely spike.
genpl_postcomp <- function(optpar, optval, x, s, par_init, fix_par, scale_factor,
                           fixed_comp = "none", fixed_c = NULL,
                           ratio_ref = NULL, ratio_k = NULL) {
  llik <- -optval #move back to positive log likelihood

  # Put the parameters back into the canonical
  # c(alpha_pos, alpha_neg, beta_pos, beta_neg, mu) form, so that everything
  # downstream of here is spared the constrained parameterizations.
  optpar <- genpl_undo_fixedcomp(optpar, fixed_comp, fixed_c)
  optpar <- genpl_undo_ratio(optpar, ratio_ref, ratio_k)

  retlist <- list(par = optpar, val = llik)

  # Check the solution pi_0 = 1. Only reachable when the weights are free: a
  # pinned weight lies strictly inside (0, 1), so the spike cannot take all of
  # the mass.
  weights_estimated <- identical(fixed_comp, "none") &&
    (!fix_par[1] || !fix_par[2])
  fix_mu  <- fix_par[length(fix_par)]
  if (weights_estimated && fix_mu) { #want to check that we can modify the weights, and mu shouldn't change
    pi_0_llik <- sum(-0.5 * log(2 * pi * s^2) - 0.5 * (x - par_init$mu)^2 / s^2)
    pi_0_llik <- pi_0_llik + sum(is.finite(x)) * log(scale_factor) #add back scaling factor
    if (pi_0_llik > llik) { #Check if our manual solution is still better than optimizer, then override
      
      # Squash both slabs to exactly 0 weight
      retlist$par$alpha_pos <- -Inf
      retlist$par$alpha_neg <- -Inf
      
      #Reset rates to 0
      retlist$par$beta_pos <- 0
      retlist$par$beta_neg <- 0
      
      retlist$val <- pi_0_llik #update manual log likelihood to our calculation
    }
  }
  return(retlist)
}

genpl_summres <- function(x, s, optpar, output) {
  
  # Unpack and transform the weights (Softmax)
  logits <- c(0, optpar$alpha_pos, optpar$alpha_neg) 
  max_logit <- max(logits)
  w <- exp(logits - max_logit)
  pi_vec <- w / sum(w)
  
  pi_0    <- pi_vec[1]
  pi_pos <- pi_vec[2]
  pi_neg  <- pi_vec[3]
  
  # Unpack and transform the rates
  lambda_pos <- exp(optpar$beta_pos)
  lambda_neg <- exp(optpar$beta_neg)
  
  # Unpack the mode
  mu <- optpar$mu
  
  # Pass the clean parameters to the engine
  return(genpl_summres_untransformed(x, s, pi_0, pi_pos, pi_neg, 
                                     lambda_pos, lambda_neg, mu, output))
}

#' @importFrom ashr my_etruncnorm my_e2truncnorm
#'
genpl_summres_untransformed <- function(x, s, pi_0, pi_pos, pi_neg, lambda_pos, lambda_neg, mu, output) {
  # Center the data
  x <- x - mu
  
  # Posterior Inclusion Probabilities (PIPs)
  gammas <- pipost_genpl(x, s, pi_0, pi_pos, pi_neg, lambda_pos, lambda_neg)
  gamma_pos <- gammas$pos
  gamma_neg  <- gammas$neg
  
  post <- list()
  
  if (result_in_output(output)) {
    # Conditional Means (E[theta | Y_+]) and (E[theta | Y_-])
    # The negative slab integrates from -Inf to 0
    E_pos <- my_etruncnorm(0, Inf, x - s^2 * lambda_pos, s)
    E_neg  <- my_etruncnorm(-Inf, 0, x + s^2 * lambda_neg, s)
    
    post$mean <- (gamma_pos * E_pos) + (gamma_neg * E_neg)
    
    # Conditional Second Moments (E[theta^2 | Y])
    E2_pos <- my_e2truncnorm(0, Inf, x - s^2 * lambda_pos, s)
    E2_neg  <- my_e2truncnorm(-Inf, 0, x + s^2 * lambda_neg, s)
    
    post$mean2 <- (gamma_pos * E2_pos) + (gamma_neg * E2_neg)
    
    post$mean2 <- pmax(post$mean2, post$mean^2) 
    # We mathematically require E[theta^2] >= (E[theta])^2 because Var(theta) >= 0.
    # Numerical errors can cause the computed second moment to dip below the squared mean. 
    # This step bounds it to ensure variance is never negative.
    
    #Handle Infinite Standard Errors, not needed but kept for consistency
    if (any(is.infinite(s))) {
      post$mean[is.infinite(s)]  <- (pi_pos / lambda_pos) - (pi_neg / lambda_neg)
      post$mean2[is.infinite(s)] <- (2 * pi_pos / lambda_pos^2) + (2 * pi_neg / lambda_neg^2)
    }
    
    # Calculate Standard Deviation and Un-center
    post$sd <- sqrt(pmax(0, post$mean2 - post$mean^2)) #not needed but kept for consistency
    post$mean2 <- post$mean2 + mu^2 + 2 * mu * post$mean #reshifted
    post$mean  <- post$mean + mu 
  }
  
  # Local False Sign Rate (lfsr)
  if ("lfsr" %in% output) {
    # Probability we are wrong is 1 minus the probability of the most likely sign
    post$lfsr <- 1 - pmax(gamma_pos, gamma_neg)
    
    if (any(is.infinite(s))) {
      post$lfsr[is.infinite(s)] <- 1 - pmax(pi_pos, pi_neg)
    }
  }
  
  return(post)
}

# Calculate posterior weights for Generalized Point-Laplace
#
#' @importFrom stats dnorm pnorm
#'
pipost_genpl <- function(x, s, pi_0, pi_pos, pi_neg, lambda_pos, lambda_neg) {
  
  # Boundary Checks
  if (pi_0 == 1) {
    return(list(pos = rep(0, length(x)), neg = rep(0, length(x))))
  }
  if (pi_pos == 1) {
    return(list(pos = rep(1, length(x)), neg = rep(0, length(x))))
  }
  if (pi_neg == 1) {
    return(list(pos = rep(0, length(x)), neg = rep(1, length(x))))
  }
  
  # Marginal Log-Likelihoods
  # Spike
  ll_0 <- dnorm(x, 0, s, log = TRUE)
  
  # Positive Slab (m_+)
  lg_pos <- log(lambda_pos) + (s^2 * lambda_pos^2 / 2) - lambda_pos * x + pnorm(x / s - s * lambda_pos, log.p = TRUE)
  
  # Negative Slab (m_-)
  lg_neg <- log(lambda_neg) + (s^2 * lambda_neg^2 / 2) + lambda_neg * x + pnorm(-x / s - s * lambda_neg, log.p = TRUE)
  
  l_0   <- log(pi_0) + ll_0
  l_pos <- log(pi_pos) + lg_pos
  l_neg <- log(pi_neg) + lg_neg
  
  # Calculate p(x_i), using log-sum-exp
  max_l <- pmax(l_0, l_pos, l_neg)
  l_tot <- max_l + log(exp(l_0 - max_l) + exp(l_pos - max_l) + exp(l_neg - max_l))
  
  # Exponentiate back to pi posterior
  pipost_pos <- exp(l_pos - l_tot)
  pipost_neg <- exp(l_neg - l_tot)
  
  return(list(pos = pipost_pos, neg = pipost_neg))
}

#Convert optimization parameters to a genlaplacemix prior object
genpl_partog <- function(par) {
  
  # Transform weights back to probabilities (Softmax)
  logits <- c(0, par$alpha_pos, par$alpha_neg) 
  max_logit <- max(logits)
  w <- exp(logits - max_logit)
  pi_vec <- w / sum(w)
  
  pi_0    <- pi_vec[1]
  pi_pos <- pi_vec[2]
  pi_neg  <- pi_vec[3]
  
  # Transform log-rates back to scales 
  scale_pos <- exp(-par$beta_pos)
  scale_neg <- exp(-par$beta_neg)
  
  # Extract the mode
  mu <- par$mu
  
  # Construct the genlaplacemix object
  if (pi_0 == 0) {
    # If there is no spike, only build the two exponential slabs
    g <- genlaplacemix(pi = c(pi_pos, pi_neg),
                       mean = rep(mu, 2),
                       scale_pos = c(scale_pos, 0),
                       scale_neg = c(0, scale_neg))
  } else {
    # Otherwise, build the full 3-component mixture
    g <- genlaplacemix(pi = c(pi_0, pi_pos, pi_neg),
                       mean = rep(mu, 3),
                       scale_pos = c(0, scale_pos, 0),
                       scale_neg = c(0, 0, scale_neg))
  }
  return(g)
}

genpl_postsamp <- function(x, s, optpar, nsamp) {
    
  # Unpack and transform the weights (Softmax)
  logits <- c(0, optpar$alpha_pos, optpar$alpha_neg) 
  max_logit <- max(logits)
  w <- exp(logits - max_logit)
  pi_vec <- w / sum(w)
  
  pi_0    <- pi_vec[1]
  pi_pos <- pi_vec[2]
  pi_neg  <- pi_vec[3]
  
  # Unpack and transform the rates
  lambda_pos <- exp(optpar$beta_pos)
  lambda_neg <- exp(optpar$beta_neg)
  
  # Extract mode
  mu <- optpar$mu
  
  return(genpl_postsamp_untransformed(x, s, pi_0, pi_pos, pi_neg, 
                                      lambda_pos, lambda_neg, mu, nsamp))
}

#' @importFrom truncnorm rtruncnorm
#' @importFrom stats runif
#'
genpl_postsamp_untransformed <- function(x, s, pi_0, pi_pos, pi_neg, lambda_pos, lambda_neg, mu, nsamp) {
  # Center the data
  x <- x - mu
  
  # Get the Posterior Inclusion Probabilities (gammas)
  gammas <- pipost_genpl(x, s, pi_0, pi_pos, pi_neg, lambda_pos, lambda_neg)
  gamma_pos <- gammas$pos
  gamma_neg  <- gammas$neg
  
  nobs <- length(x) 
  
  # Draw a uniform random number between 0 and 1 for every single sample
  U <- runif(nsamp * nobs)
  
  # Expand the gammas so they match the number of samples being drawn
  rep_gamma_pos <- rep(gamma_pos, each = nsamp)
  rep_gamma_neg  <- rep(gamma_neg, each = nsamp)
  
  # If U is less than gamma_pos, it belongs to the positive slab
  is_positive <- U < rep_gamma_pos
  
  # If U is between gamma_pos and (gamma_pos + gamma_neg), it belongs to negative slab (so likelihood of neg gamma)
  is_negative <- (U >= rep_gamma_pos) & (U < (rep_gamma_pos + rep_gamma_neg))
  
  # Remaining values are assigned to spike
  
  if (length(s) == 1) {
    s <- rep(s, nobs) #homoskedastic (makes it repeat for all samples)
  }
  # Generate samples for both slabs for theta
  # Uses that p(theta|x,y+) is truncated normal with z_+,- as mean, s as sd
  negative_samp <- mapply(FUN = function(mean, sd) {
    rtruncnorm(nsamp, -Inf, 0, mean, sd)
  }, mean = x + s^2 * lambda_neg, sd = s)
  
  positive_samp <- mapply(FUN = function(mean, sd) {
    rtruncnorm(nsamp, 0, Inf, mean, sd)
  }, mean = x - s^2 * lambda_pos, sd = s)
  
  # Initialize the final matrix with exactly 0 (representing the spike)
  samp <- matrix(0, nrow = nsamp, ncol = nobs)
  
  # Swap in the positive & negative samples where theta is from positive or negative slab respectively.
  samp[is_positive] <- positive_samp[is_positive]
  samp[is_negative] <- negative_samp[is_negative]
  
  # Un-center the entire matrix
  samp <- samp + mu
  
  return(samp)
}

#Below is one final function for Gen. Point Laplace necessary because we have 5 parameters instead of just 2 or 3

# Translate the ebnm interface (pointmass/scale/mode) into the generalized
#   optimization interface. The other parametric families have three
#   parameters, so parametric_workhorse builds fix_par itself; the generalized
#   point-Laplace family has five:
#
#     c(alpha_pos, alpha_neg, beta_pos, beta_neg, mu)
#
#   The weights are parameterized as a softmax with the point mass as the
#   reference category, so slots 1 and 2 are both weight parameters and slots
#   3 and 4 are both scale parameters.
#

# Normalize the scale argument into the two scale values (NA where a scale is
#   to be estimated) and a logical vector marking which are fixed. Shared by
#   genpl_fixpar and genpl_initpar so that the two cannot disagree about what a
#   given scale argument means.
#
genpl_parse_scale <- function(scale) {
  if (identical(scale, "estimate")) {
    return(list(value = c(NA_real_, NA_real_), fix = c(FALSE, FALSE)))
  }
  if (!is.numeric(scale) || !(length(scale) %in% c(1, 2))) {
    stop("Argument 'scale' must be 'estimate' or a numeric vector of length ",
         "one or two (positive scale, negative scale), in which NA marks a ",
         "scale that is to be estimated.")
  }
  value <- rep(as.numeric(scale), length.out = 2)
  fix   <- !is.na(value)
  if (!any(fix)) {
    stop("Argument 'scale' must fix at least one scale; use ",
         "scale = 'estimate' to estimate both.")
  }
  if (any(value[fix] <= 0)) {
    stop("Fixed scales must be positive.")
  }
  return(list(value = value, fix = fix))
}

#   Unlike the other families, the constraints here are rich enough that this
#   returns a list rather than a bare fix_par vector: the weight constraint
#   changes the number of parameters, and the tail-heaviness constraint changes
#   what one of the slots means. parametric_workhorse unpacks it and forwards
#   the extra fields to genpl_initpar and genpl_precomp.
#
genpl_fixpar <- function(pointmass, scale, mode, pi_fixed = NULL,
                         scale_ratio = NULL) {
  # Weights. Pinning one weight drops a parameter: the two that remain remain
  # share a single free log-ratio, so there is one alpha rather than two.
  if (is.null(pi_fixed)) {
    fixed_comp <- "none"
    fixed_c    <- NULL
    fix_alpha  <- c(!pointmass, !pointmass)
  } else {
    if (!all(c("comp", "value") %in% names(pi_fixed))) {
      stop("Argument 'pi_fixed' must be a list with elements 'comp' and ",
           "'value'.")
    }
    fixed_comp <- match.arg(pi_fixed$comp, c("spike", "pos", "neg"))
    fixed_c    <- pi_fixed$value
    if (!is.numeric(fixed_c) || length(fixed_c) != 1
        || fixed_c <= 0 || fixed_c >= 1) {
      stop("Element 'value' of 'pi_fixed' must be a single number strictly ",
           "between zero and one.")
    }
    fix_alpha <- FALSE
  }

  fix_beta <- genpl_parse_scale(scale)$fix

  # Tail-heaviness constraint. The fixed tail is the reference; the estimated
  # one is bounded by scale_ratio times its scale, and so stays free (its slot
  # carries the unconstrained delta of genpl_nllik).
  if (is.null(scale_ratio)) {
    ratio_ref <- NULL
    ratio_k   <- NULL
  } else {
    if (!is.numeric(scale_ratio) || length(scale_ratio) != 1
        || scale_ratio <= 0) {
      stop("Argument 'scale_ratio' must be a single positive number.")
    }
    if (all(fix_beta) || !any(fix_beta)) {
      stop("Argument 'scale_ratio' requires exactly one of the two scales to ",
           "be fixed, which 'scale' specifies by giving one number and one NA.")
    }
    ratio_ref <- if (fix_beta[1]) "pos" else "neg"
    ratio_k   <- scale_ratio
  }

  return(list(fix_par    = c(fix_alpha,                       # alpha(s)
                             fix_beta[1],                     # beta_pos
                             fix_beta[2],                     # beta_neg
                             !identical(mode, "estimate")),   # mu
              fixed_comp = fixed_comp,
              fixed_c    = fixed_c,
              ratio_ref  = ratio_ref,
              ratio_k    = ratio_k))
}

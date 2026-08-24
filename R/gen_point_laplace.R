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
    } else { #Calculate this way to avoid errors?
      # Upper tail CDF: P(X > q)
      #want to add components from each tail
      if (any(pos_idx)) { #if q < 0, pos slab has all weight. else, we use CDF of exp, with lower_tail FALSe since we go other way
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
#
genpl_checkg <- function(g_init, fix_g, mode, scale, pointmass, call) {
  check_g_init(g_init = g_init,
               fix_g = fix_g,
               mode = mode,
               scale = scale,
               pointmass = pointmass,
               call = call,
               class_name = "genlaplacemix",
               scale_name = c("rate_pos", "rate_neg"))
}

#' @param x A vector of observations.
#' @param s A vector of standard errors.
#' @param g_init An optional initial genlaplacemix object.
#' @param fix_g Boolean indicating if the prior is completely fixed.
#' @param pointmass Boolean indicating if a point mass at zero is allowed.
#' @param scale_pos The scale for the positive tail ("estimate" or numeric).
#' @param scale_neg The scale for the negative tail ("estimate" or numeric).
#'
#' @return A numeric vector of unconstrained parameters for the optimizer.

genpl_initpar <- function(x, s, g_init, fix_g, pointmass, scale_pos, scale_neg) {
  par <- list()
  
  # SCENARIO A: The user provided an initial prior (g_init) and it has exactly 3 components
  if (!is.null(g_init) && length(g_init$pi) == 3) {
    # Extract raw weights and regularize
    # Add a tiny epsilon to ensure no weight is exactly 0, preventing log(0)
    eps <- 1e-8
    pi_vec <- c(g_init$pi[1], g_init$pi[2], g_init$pi[3]) + eps
    pi_vec <- pi_vec / sum(pi_vec) # Re-normalize so they sum exactly to 1
    
    pi_0    <- pi_vec[1]
    pi_plus <- pi_vec[2]
    pi_neg  <- pi_vec[3]
    
    # Map weights to be unconstrained (anchored to the spike pi_0)
    par$logit_w_pos <- log(pi_plus / pi_0)
    par$logit_w_neg <- log(pi_neg / pi_0)
    
    # Rates: Convert finite rates to unconstrained log-rates
    # In genlaplacemix, Component 2 is the positive slab, Component 3 is the negative slab
    par$log_rate_pos <- log(g_init$rate_pos[2])
    par$log_rate_neg <- log(g_init$rate_neg[3])
    
    par$mu <- g_init$mean[1]
    
  } else {
    # SCENARIO B: Heuristic Empirical Initialization
    
    # Handle Weights (Pointmass check)
    if (!pointmass) {
      # No point mass means we split the weight 50/50 between the positive and negative slabs. 
      # We use a tiny epsilon for pi_0 to prevent log(0) = Inf
      pi_0    <- 1e-8 
      pi_plus <- 0.50
      pi_neg  <- 0.50
      
    } else {
      # Default: 50% spike, split the rest, same as point laplace
      pi_0    <- 0.50
      pi_plus <- 0.25
      pi_neg  <- 0.25
    } 
    
    par$logit_w_pos <- log(pi_plus / pi_0)
    par$logit_w_neg <- log(pi_neg / pi_0)
    
    # Rates initialization
    # Estimate empirical rate by moment matching
    excess_var <- max(var(x, na.rm = TRUE) - mean(s^2, na.rm = TRUE), 1e-4) 
    empirical_rate <- 1 / sqrt(excess_var)
    
    if (!identical(scale_pos, "estimate")) {
      if (length(scale_pos) != 1) stop("Argument 'scale_pos' must be a scalar.")
      par$log_rate_pos <- -log(scale_pos)
    } else {
      par$log_rate_pos <- log(empirical_rate)
    }
    
    # Check Negative Tail
    if (!identical(scale_neg, "estimate")) {
      if (length(scale_neg) != 1) stop("Argument 'scale_neg' must be a scalar.")
      par$log_rate_neg <- -log(scale_neg)
    } else {
      par$log_rate_neg <- log(empirical_rate)
    }
    # Mode
    if (!identical(mode, "estimate")) {
      par$mu <- mode
    } else {
      par$mu <- mean(x) # Center of the data for a two-sided distribution
    }
  }
  
  # Return as a flat numeric vector for L-BFGS-B
  return(unlist(par))
}

genpl_scalepar <- function(par, scale_factor) {
  # Adjust the positive log-rate
  if (!is.null(par$log_rate_pos)) {
    par$log_rate_pos <- par$log_rate_pos - log(scale_factor)
  }
  # Adjust the negative log-rate
  if (!is.null(par$log_rate_neg)) {
    par$log_rate_neg <- par$log_rate_neg - log(scale_factor)
  }
  
  # Adjust the mode (if you choose to estimate it)
  if (!is.null(par$mu)) {
    par$mu <- scale_factor * par$mu
  }
  return(par)
}
    
# No precomputations are done for generalized point-Laplace. ??
# Q: What does this function do?

genpl_precomp <- function(x, s, par_init, fix_par) {
  # Check if the mode (mu) is fixed. 
  # Assuming your parameter vector is ordered: 
  # c(logit_w_pos, logit_w_neg, log_rate_pos, log_rate_neg, mu)
  # If fix_par is a named vector, fix_par["mu"] is even safer
  
  fix_mu  <- fix_par[5]
  
  if (!fix_mu && any(s == 0)) {
    stop("The mode cannot be estimated if any SE is zero (the gradient does ",
         "not exist).")
  }
  
  return(NULL)
}    
  
# The negative log likelihood.
#
#' @importFrom stats pnorm
#'
genpl_nllik <- function(par, x, s, par_init, fix_par) {
  
  # Unpack and Transform Parameters
  p <- unlist(par_init)
  p[!fix_par] <- par
  
  # Unconstrained parameters
  logit_w_pos  <- p[1]
  logit_w_neg  <- p[2]
  log_rate_pos <- p[3]
  log_rate_neg <- p[4]
  mu           <- p[5]
  
  # Transform weights (Softmax)
  logits <- c(0, logit_w_pos, logit_w_neg) 
  max_logit <- max(logits)
  w <- exp(logits - max_logit)
  pi_vec <- w / sum(w)
  
  pi_0    <- pi_vec[1]
  pi_plus <- pi_vec[2]
  pi_neg  <- pi_vec[3]
  
  # Transform rates
  rate_pos <- exp(log_rate_pos)
  rate_neg <- exp(log_rate_neg)
  
  # Log-Likelihoods
  ll_0 <- -0.5 * log(2 * pi * s^2) - 0.5 * (x - mu)^2 / s^2 #Spike
  
  # The Positive Slab (Right Tail)
  xright <- (x - mu) / s - s * rate_pos
  lpnormright <- pnorm(xright, log.p = TRUE)
  ll_plus <- log(rate_pos) + (s^2 * rate_pos^2 / 2) - rate_pos * (x - mu) + lpnormright
  
  #The Negative Slab
  xleft <- -(x - mu) / s - s * rate_neg
  lpnormleft <- pnorm(xleft, log.p = TRUE)
  ll_neg <- log(rate_neg) + (s^2 * rate_neg^2 / 2) + rate_neg * (x - mu) + lpnormleft
  
  #Use Log Sum Exp to add all 3 arguments
  ll_total_slab <- logscale_add(log(pi_plus) + ll_plus, log(pi_neg) + ll_neg)
  llik <- logscale_add(log(pi_0) + ll_0, ll_total_slab)
  nllik <- -sum(llik)
  
  
  #want to do gradients and hessians for these next?? can i just use optim?
  return(nllik)
}



logscale_add <- function(log.x, log.y) {
  C <- pmax(log.x, log.y)
  return(log(exp(log.x - C) + exp(log.y - C)) + C)
  }
  
# Postcomputations: check boundary solutions.
genpl_postcomp <- function(optpar, optval, x, s, par_init, fix_par, scale_factor) {
  llik <- -optval #move back to positive log likelihood
  retlist <- list(par = optpar, val = llik)
  
  # Check the solution pi0 = 1.
  weights_estimated <- !fix_par[1] || !fix_par[2]
  fix_mu  <- fix_par[5]
  if (!fix_pi0 && fix_mu) {
    pi0_llik <- sum(-0.5 * log(2 * pi * s^2) - 0.5 * (x - par_init$mu)^2 / s^2)
    pi0_llik <- pi0_llik + sum(is.finite(x)) * log(scale_factor) #add back scaling factor
    if (pi0_llik > llik) { #The override (why do we need this?)
      
      # Squash both slabs to exactly 0 weight
      retlist$par$logit_w_pos <- -Inf
      retlist$par$logit_w_neg <- -Inf
      
      #Reset rates to 0
      retlist$par$log_rate_pos <- 0
      retlist$par$log_rate_neg <- 0
      
      retlist$val <- pi0_llik #update manual log likelihood to our calculation
    }
  }
  return(retlist)
}

genpl_summres <- function(x, s, optpar, output) {
  
  # Unpack and transform the weights (Softmax)
  logits <- c(0, optpar$logit_w_pos, optpar$logit_w_neg) 
  max_logit <- max(logits)
  w <- exp(logits - max_logit)
  pi_vec <- w / sum(w)
  
  pi_0    <- pi_vec[1]
  pi_plus <- pi_vec[2]
  pi_neg  <- pi_vec[3]
  
  # Unpack and transform the rates
  rate_pos <- exp(optpar$log_rate_pos)
  rate_neg <- exp(optpar$log_rate_neg)
  
  # Unpack the mode
  mu <- optpar$mu
  
  # Pass the clean parameters to the engine
  return(genpl_summres_untransformed(x, s, pi_0, pi_plus, pi_neg, 
                                     rate_pos, rate_neg, mu, output))
}

#' @importFrom ashr my_etruncnorm my_e2truncnorm
#'
genpl_summres_untransformed <- function(x, s, pi_0, pi_plus, pi_neg, rate_pos, rate_neg, mu, output) {
  # Center the data
  x <- x - mu
  
  # Get the Posterior Inclusion Probabilities (PIPs)
  gammas <- wpost_genpl(x, s, pi_0, pi_plus, pi_neg, rate_pos, rate_neg)
  gamma_plus <- gammas$pos #post inclusion prob for pos slab
  gamma_neg  <- gammas$neg #post inclusion prob for neg slab
  
  post <- list()
  
  if (result_in_output(output)) {
    # Calculate Conditional Means (E[theta | Y_+]) and (E[theta | Y_-])
    # The negative slab integrates from -Inf to 0
    E_plus <- my_etruncnorm(0, Inf, x - s^2 * rate_pos, s)
    E_neg  <- my_etruncnorm(-Inf, 0, x + s^2 * rate_neg, s)
    
    post$mean <- (gamma_plus * E_plus) + (gamma_neg * E_neg)
    
    # Calculate Conditional Second Moments (E[theta^2 | Y])
    E2_plus <- my_e2truncnorm(0, Inf, x - s^2 * rate_pos, s)
    E2_neg  <- my_e2truncnorm(-Inf, 0, x + s^2 * rate_neg, s)
    
    post$mean2 <- (gamma_plus * E2_plus) + (gamma_neg * E2_neg)
    post$mean2 <- pmax(post$mean2, post$mean^2) #why do we do this?
    
    #Handle Infinite Standard Errors
    if (any(is.infinite(s))) {
      post$mean[is.infinite(s)]  <- (pi_plus / rate_pos) - (pi_neg / rate_neg)
      post$mean2[is.infinite(s)] <- (2 * pi_plus / rate_pos^2) + (2 * pi_neg / rate_neg^2)
    }
    
    # Calculate Standard Deviation and Un-center
    post$sd <- sqrt(pmax(0, post$mean2 - post$mean^2)) #Just definition of variance
    post$mean2 <- post$mean2 + mu^2 + 2 * mu * post$mean #reshifted
    post$mean  <- post$mean + mu 
  }
  
  # Local False Sign Rate (LFSR)
  if ("lfsr" %in% output) {
    # Probability we are wrong is 1 minus the probability of the most likely sign
    post$lfsr <- 1 - pmax(gamma_plus, gamma_neg)
    
    if (any(is.infinite(s))) {
      post$lfsr[is.infinite(s)] <- 1 - pmax(pi_plus, pi_neg)
    }
  }
  
  return(post)
}

#' Calculate posterior weights for Generalized Point-Laplace
#'
#' @importFrom stats dnorm pnorm
#'
pipost_genpl <- function(x, s, pi_pos, pi_neg, lambda_pos, lambda_neg) {
  
  pi_0 <- 1 - pi_pos - pi_neg
  
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
  lf <- dnorm(x, 0, s, log = TRUE)
  
  # Positive Slab (m_+)
  lg_pos <- log(lambda_pos) + (s^2 * lambda_pos^2 / 2) - lambda_pos * x + pnorm(x / s - s * lambda_pos, log.p = TRUE)
  
  # Negative Slab (m_-)
  lg_neg <- log(lambda_neg) + (s^2 * lambda_neg^2 / 2) + lambda_neg * x + pnorm(-x / s - s * lambda_neg, log.p = TRUE)
  
  l_0   <- log(pi_0) + lf
  l_pos <- log(pi_pos) + lg_pos
  l_neg <- log(pi_neg) + lg_neg
  
  # Calculate p(x_i), using log-sum-exp
  max_l <- pmax(l_0, l_pos, l_neg)
  l_tot <- max_l + log(exp(l_0 - max_l) + exp(l_pos - max_l) + exp(l_neg - max_l))
  
  # Exponentiate back to pi posterior
  pipost_pos <- exp(l_pos - l_tot)
  pipost_neg <- exp(l_neg - l_tot)
  
  return(list(pos = lambdapost_pos, neg = lambdapost_neg))
}

#Convert optimization parameters to a genlaplacemix prior object
genpl_partog <- function(par) {
  
  # Transform weights back to probabilities (Softmax)
  logits <- c(0, par$logit_w_pos, par$logit_w_neg) 
  max_logit <- max(logits)
  w <- exp(logits - max_logit)
  pi_vec <- w / sum(w)
  
  pi_0    <- pi_vec[1]
  pi_plus <- pi_vec[2]
  pi_neg  <- pi_vec[3]
  
  # Transform log-rates back to scales 
  scale_pos <- exp(-par$log_rate_pos)
  scale_neg <- exp(-par$log_rate_neg)
  
  # Extract the mode
  mu <- par$mu
  
  # Construct the genlaplacemix object
  if (pi_0 == 0) {
    # If there is no spike, only build the two exponential slabs
    g <- genlaplacemix(pi = c(pi_plus, pi_neg),
                       mean = rep(mu, 2),
                       scale_pos = c(scale_pos, 0),
                       scale_neg = c(0, scale_neg))
  } else {
    # Otherwise, build the full 3-component mixture
    g <- genlaplacemix(pi = c(pi_0, pi_plus, pi_neg),
                       mean = rep(mu, 3),
                       scale_pos = c(0, scale_pos, 0),
                       scale_neg = c(0, 0, scale_neg))
  
  return(g)
}
}

genpl_postsamp <- function(x, s, optpar, nsamp) {
    
  # Unpack and transform the weights (Softmax)
  logits <- c(0, optpar$logit_w_pos, optpar$logit_w_neg) 
  max_logit <- max(logits)
  w <- exp(logits - max_logit)
  pi_vec <- w / sum(w)
  
  pi_0    <- pi_vec[1]
  pi_plus <- pi_vec[2]
  pi_neg  <- pi_vec[3]
  
  # Unpack and transform the rates
  lambda_pos <- exp(optpar$log_rate_pos)
  lambda_neg <- exp(optpar$log_rate_neg)
  
  # Extract mode
  mu <- optpar$mu
  
  return(genpl_postsamp_untransformed(x, s, pi_0, pi_plus, pi_neg, 
                                      lambda_pos, lambda_neg, mu, nsamp))
}

#' @importFrom truncnorm rtruncnorm
#' @importFrom stats runif
#'
genpl_postsamp_untransformed <- function(x, s, pi_0, pi_plus, pi_neg, lambda_pos, lambda_neg, mu, nsamp) {
  # Center the data
  x <- x - mu
  
  # Get the Posterior Inclusion Probabilities (gammas)
  gammas <- wpost_genpl(x, s, pi_plus, pi_neg, lambda_pos, lambda_neg)
  gamma_plus <- gammas$pos
  gamma_neg  <- gammas$neg
  
  nobs <- length(gamma_plus)
  
  # ROLL THE 3-SIDED DIE
  # Draw a uniform random number between 0 and 1 for every single sample
  U <- runif(nsamp * nobs)
  
  # Expand the gammas so they match the number of samples being drawn
  rep_gamma_plus <- rep(gamma_plus, each = nsamp)
  rep_gamma_neg  <- rep(gamma_neg, each = nsamp)
  
  # If U is less than gamma_plus, it belongs to the positive slab
  is_positive <- U < rep_gamma_plus
  
  # If U is between gamma_plus and (gamma_plus + gamma_neg), negative slab (so likelihood of neg gamma)
  is_negative <- (U >= rep_gamma_plus) & (U < (rep_gamma_plus + rep_gamma_neg))
  
  # Remaining values are assigned to spike
  
  if (length(s) == 1) {
    s <- rep(s, nobs) #homoskedastic makes it repeat for all samples
  }
  
  # Generate the truncated random samples for both slabs
  # Uses the exact same shifted base means we derived for the conditional expectations, rtruncnorm is not vectorized.
  negative_samp <- mapply(FUN = function(mean, sd) {
    rtruncnorm(nsamp, -Inf, 0, mean, sd)
  }, mean = x + s^2 * lambda_neg, sd = s)
  
  positive_samp <- mapply(FUN = function(mean, sd) {
    rtruncnorm(nsamp, 0, Inf, mean, sd)
  }, mean = x - s^2 * lambda_pos, sd = s)
  
  # Initialize the final matrix with exactly 0 (representing the spike)
  samp <- matrix(0, nrow = nsamp, ncol = nobs)
  
  # Swap in the positive & negative samples where the die rolled positive & negative respectively.
  samp[is_positive] <- positive_samp[is_positive]
  samp[is_negative] <- negative_samp[is_negative]
  
  # Un-center the entire matrix
  samp <- samp + mu
  
  return(samp)
}
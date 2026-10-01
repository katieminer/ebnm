#' Negative Log-Likelihood for Generalized Point-Laplace for all cases
#'
#' @description 
#' Evaluates the negative log-likelihood, gradient, and Hessian for an EBNM 
#' generalized point-Laplace prior. Automatically branches between the 
#' full 5-slot parameterization (all weights free) and the 4-slot 
#' parameterization (one weight fixed).
#' 
#' @inheritParams ebnm
#' @param par_init List of all initial parameter values.
#' @param fix_par Logical vector indicating which parameters are fixed.
#' @param fixed_comp Character: "none" (default), "spike", "pos", or "neg".
#' @param fixed_c Numeric constant in [0, 1] representing the fixed weight. Required if fixed_comp != "none".
#' @param ratio_ref Character or NULL. Which tail has its scale fixed and so
#'   serves as the reference for the tail-heaviness constraint: "pos" or "neg".
#'   NULL (default) means unconstrained.
#' @param ratio_k Positive numeric. Bounds the free tail's scale by ratio_k
#'   times the fixed tail's scale, i.e. the free tail can be at most ratio_k
#'   times as heavy. Required if ratio_ref is non-NULL.
#'
#' @importFrom stats pnorm plogis
#' @export
genpl_nllik <- function(par, x, s, par_init, fix_par, fixed_comp = "none", fixed_c = NULL, ratio_ref = NULL, ratio_k = NULL, calc_grad = TRUE, calc_hess = FALSE) {

  #Unpack & Transform Parameters
  p <- unlist(par_init)
  p[!fix_par] <- par

  if (fixed_comp == "none") {
    # All weights of the prior are free or all fixed, so either estimate both alphas or none
    fix_alpha_pos <- fix_par[1] 
    fix_alpha_neg <- fix_par[2] #always the same as fix_par[1]
    fix_beta_pos  <- fix_par[3]
    fix_beta_neg  <- fix_par[4]
    fix_mu        <- fix_par[5]
    
    alpha_pos <- p[1]
    alpha_neg <- p[2]
    beta_pos  <- p[3]
    beta_neg  <- p[4]
    mu        <- p[5]
    
    # Standard Softmax
    logits <- c(0, alpha_pos, alpha_neg) 
    max_logit <- max(logits)
    w <- exp(logits - max_logit)
    pi_vec <- w / sum(w)
    
    pi_0   <- pi_vec[1]
    pi_pos <- pi_vec[2]
    pi_neg <- pi_vec[3]
    
  } else {
    # When one component of pi is fixed, we only define one unconstrained alpha
    fix_alpha    <- fix_par[1]
    fix_beta_pos <- fix_par[2]
    fix_beta_neg <- fix_par[3]
    fix_mu       <- fix_par[4]
    
    alpha    <- p[1]
    beta_pos <- p[2]
    beta_neg <- p[3]
    mu       <- p[4]
    
    # Convert the unconstrained log-ratio (alpha) into relative proportions (f_1, f_2) that sum to 1
    exp_a <- exp(alpha)
    f_1 <- exp_a / (1 + exp_a) 
    f_2 <- 1 / (1 + exp_a)     
    
    #Scale these proportions by the remaining probability mass (1 - c) to get the final weights
    rem_c <- 1 - fixed_c
    
    if (fixed_comp == "spike") {
      pi_0   <- fixed_c
      pi_pos <- f_1 * rem_c  #f_+ = pi_+/(1-c)  
      pi_neg <- f_2 * rem_c    
    } else if (fixed_comp == "pos") {
      pi_pos <- fixed_c
      pi_0   <- f_1 * rem_c
      pi_neg <- f_2 * rem_c
    } else if (fixed_comp == "neg") {
      pi_neg <- fixed_c
      pi_0   <- f_1 * rem_c
      pi_pos <- f_2 * rem_c
    }
  }
  
  # Tail-heaviness constraint. The free tail's slot carries delta rather than a
  # log-rate: beta_free = B + softplus(delta), where B = beta_ref - log(ratio_k)
  # bounds its scale at ratio_k times the reference scale. softplus > 0 enforces
  # the bound, and B is constant because the reference log-rate is fixed, so no
  # cross-derivative runs through it.
  if (!is.null(ratio_ref)) {
    if (identical(ratio_ref, "pos")) {
      delta    <- beta_neg
      beta_neg <- (beta_pos - log(ratio_k)) + genpl_softplus(delta)
    } else {
      delta    <- beta_pos
      beta_pos <- (beta_neg - log(ratio_k)) + genpl_softplus(delta)
    }
  }

  # Transform rates
  lambda_pos <- exp(beta_pos)
  lambda_neg <- exp(beta_neg)
  
  #  ==========================================
  #  ==========================================
  
  # Log-Likelihood Evaluation
  ll_0 <- -0.5 * log(2 * pi * s^2) - 0.5 * (x - mu)^2 / s^2 
  
  z_pos <- (x - mu) / s - s * lambda_pos
  lpnorm_pos <- pnorm(z_pos, log.p = TRUE)
  ll_pos <- log(lambda_pos) + (s^2 * lambda_pos^2 / 2) - lambda_pos * (x - mu) + lpnorm_pos
  
  z_neg <- -(x - mu) / s - s * lambda_neg
  lpnorm_neg <- pnorm(z_neg, log.p = TRUE)
  ll_neg <- log(lambda_neg) + (s^2 * lambda_neg^2 / 2) + lambda_neg * (x - mu) + lpnorm_neg
  
  ll_total_slab <- logscale_add(log(pi_pos) + ll_pos, log(pi_neg) + ll_neg)
  llik <- logscale_add(log(pi_0) + ll_0, ll_total_slab)
  nllik <- -sum(llik)
  
  # ==========================================
  # Gradients
  # ==========================================
  gamma_0   <- exp(log(pi_0) + ll_0 - llik)
  gamma_pos <- exp(log(pi_pos) + ll_pos - llik)
  gamma_neg <- exp(log(pi_neg) + ll_neg - llik)
  
  # Robust inverse Mills ratio: exp(dnorm(z, log=T) - pnorm(z, log=T)).
  # Computed unconditionally: the Hessian block needs these via B_c and D_c
  # whichever parameters are free, and they are only two vectorized exp calls.
  dlogpnorm.pos <- exp(-log(2 * pi) / 2 - z_pos^2 / 2 - lpnorm_pos)
  dlogpnorm.neg <- exp(-log(2 * pi) / 2 - z_neg^2 / 2 - lpnorm_neg)
  
  if (calc_grad || calc_hess) {
    grad <- numeric(length(par))
    i <- 1
    
    if (fixed_comp == "none") {
      if (!fix_alpha_pos) {
        grad[i] <- sum(-(gamma_pos - pi_pos))
        i <- i + 1
      }
      if (!fix_alpha_neg) {
        grad[i] <- sum(-(gamma_neg - pi_neg))
        i <- i + 1
      }
    } else {
      if (!fix_alpha) { #fix_alpha should always be False
        if (fixed_comp == "spike") {
          dnllik.dalpha <- -(gamma_pos * f_2 - gamma_neg * f_1)
        } else if (fixed_comp == "pos") {
          dnllik.dalpha <- -(gamma_0 * f_2 - gamma_neg * f_1)
        } else if (fixed_comp == "neg") {
          dnllik.dalpha <- -(gamma_0 * f_2 - gamma_pos * f_1)
        }
        grad[i] <- sum(dnllik.dalpha)
        i <- i + 1
      }
    }
    
    if (!fix_beta_pos) {
      dgpos.dlambda <- exp(ll_pos - llik) * (1 / lambda_pos + lambda_pos * s^2 - (x - mu) - s * dlogpnorm.pos)
      grad[i] <- sum(-pi_pos * dgpos.dlambda * lambda_pos)
      i <- i + 1
    }
    
    if (!fix_beta_neg) {
      dgneg.dlambda <- exp(ll_neg - llik) * (1 / lambda_neg + lambda_neg * s^2 + (x - mu) - s * dlogpnorm.neg)
      grad[i] <- sum(-pi_neg * dgneg.dlambda * lambda_neg)
      i <- i + 1
    }
    
    if (!fix_mu) {
      df.dmu    <- exp(ll_0 - llik) * ((x - mu) / s^2)
      dgpos.dmu <- exp(ll_pos - llik) * (lambda_pos - dlogpnorm.pos / s)
      dgneg.dmu <- exp(ll_neg - llik) * (-lambda_neg + dlogpnorm.neg / s)
      grad[i] <- sum(-pi_0 * df.dmu - pi_pos * dgpos.dmu - pi_neg * dgneg.dmu)
    }
    
    attr(nllik, "gradient") <- grad
  }
  
  # ==========================================
  # Hessians
  # ==========================================
  if (calc_hess) {
    hess <- matrix(0, nrow = length(par), ncol = length(par))
    n_obs <- length(x)
    
    d2logpnorm.pos <- -dlogpnorm.pos * (z_pos + dlogpnorm.pos)
    d2logpnorm.neg <- -dlogpnorm.neg * (z_neg + dlogpnorm.neg)
    
    #B_c = dlog(m_c(x))/dbeta_c
    B_pos <- 1 - lambda_pos * (x - mu) + lambda_pos^2 * s^2 - s * lambda_pos * dlogpnorm.pos
    B_neg <- 1 + lambda_neg * (x - mu) + lambda_neg^2 * s^2 - s * lambda_neg * dlogpnorm.neg
    
    #D_c = dlog(m_c(x))/dmu
    D_0   <- (x - mu) / s^2 
    D_pos <- lambda_pos - dlogpnorm.pos / s
    D_neg <- -lambda_neg + dlogpnorm.neg / s
    D_bar <- gamma_0 * D_0 + gamma_pos * D_pos + gamma_neg * D_neg
    
    i <- 1
    
    if (fixed_comp == "none") {
      # Original 2 Alpha Hessian Block
      if (!fix_alpha_pos) {
        hess[i, i] <- - (sum(gamma_pos * (1 - gamma_pos)) - n_obs * pi_pos * (1 - pi_pos))
        j <- i + 1
        if (!fix_alpha_neg) {
          hess[i, j] <- hess[j, i] <- - (sum(-gamma_pos * gamma_neg) + n_obs * pi_pos * pi_neg)
          j <- j + 1
        }
        if (!fix_beta_pos) {
          hess[i, j] <- hess[j, i] <- - sum(gamma_pos * (1 - gamma_pos) * B_pos)
          j <- j + 1
        }
        if (!fix_beta_neg) {
          hess[i, j] <- hess[j, i] <- - sum(-gamma_pos * gamma_neg * B_neg)
          j <- j + 1
        }
        if (!fix_mu) {
          hess[i, j] <- hess[j, i] <- - sum(gamma_pos * (D_pos - D_bar))
        }
        i <- i + 1
      }
      
      if (!fix_alpha_neg) {
        hess[i, i] <- - (sum(gamma_neg * (1 - gamma_neg)) - n_obs * pi_neg * (1 - pi_neg))
        j <- i + 1
        if (!fix_beta_pos) {
          hess[i, j] <- hess[j, i] <- - sum(-gamma_neg * gamma_pos * B_pos)
          j <- j + 1
        }
        if (!fix_beta_neg) {
          hess[i, j] <- hess[j, i] <- - sum(gamma_neg * (1 - gamma_neg) * B_neg)
          j <- j + 1
        }
        if (!fix_mu) {
          hess[i, j] <- hess[j, i] <- - sum(gamma_neg * (D_neg - D_bar))
        }
        i <- i + 1
      }
      
    } else {
      # Hessian Block for Only One Alpha
      if (!fix_alpha) { #fix_alpha should always be False
        j <- i + 1
        if (fixed_comp == "spike") {
          d2nllik.dalpha2 <- -sum( (gamma_pos * f_2 - gamma_neg * f_1) * (f_2 * (1 - gamma_pos) - f_1 * (1 - gamma_neg)) )
          hess[i, i] <- d2nllik.dalpha2
          
          if (!fix_beta_pos) { hess[i, j] <- hess[j, i] <- - sum(gamma_pos * B_pos * ((1 - gamma_pos) * f_2 + gamma_neg * f_1)); j <- j + 1 }
          if (!fix_beta_neg) { hess[i, j] <- hess[j, i] <-   sum(gamma_neg * B_neg * (gamma_pos * f_2 + (1 - gamma_neg) * f_1)); j <- j + 1 }
          if (!fix_mu)       { hess[i, j] <- hess[j, i] <- - sum(gamma_pos * (D_pos - D_bar) * f_2 - gamma_neg * (D_neg - D_bar) * f_1) }
          
        } else if (fixed_comp == "pos") {
          
          d2nllik.dalpha2 <- -sum( (gamma_0 * f_2 - gamma_neg * f_1) * (f_2 * (1 - gamma_0) - f_1 * (1 - gamma_neg)) )
          hess[i, i] <- d2nllik.dalpha2
          
          if (!fix_beta_pos) { hess[i, j] <- hess[j, i] <- - sum(gamma_pos * B_pos * (gamma_neg * f_1 - gamma_0 * f_2)); j <- j + 1 }
          if (!fix_beta_neg) { hess[i, j] <- hess[j, i] <-   sum(gamma_neg * B_neg * (gamma_0 * f_2 + (1 - gamma_neg) * f_1)); j <- j + 1 }
          if (!fix_mu)       { hess[i, j] <- hess[j, i] <- - sum(gamma_0 * (D_0 - D_bar) * f_2 - gamma_neg * (D_neg - D_bar) * f_1) }
          
        } else if (fixed_comp == "neg") {
          
          d2nllik.dalpha2 <- -sum( (gamma_0 * f_2 - gamma_pos * f_1) * (f_2 * (1 - gamma_0) - f_1 * (1 - gamma_pos)) )
          hess[i, i] <- d2nllik.dalpha2
          
          if (!fix_beta_pos) { hess[i, j] <- hess[j, i] <-   sum(gamma_pos * B_pos * (gamma_0 * f_2 + (1 - gamma_pos) * f_1)); j <- j + 1 }
          if (!fix_beta_neg) { hess[i, j] <- hess[j, i] <- - sum(gamma_neg * B_neg * (gamma_pos * f_1 - gamma_0 * f_2)); j <- j + 1 }
          if (!fix_mu)       { hess[i, j] <- hess[j, i] <- - sum(gamma_0 * (D_0 - D_bar) * f_2 - gamma_pos * (D_pos - D_bar) * f_1) }
        }
        i <- i + 1
      }
    }
    
    #  Shared Beta & Mu Hessian Block
    if (!fix_beta_pos) {
      dB_pos_dbeta <- 2 * lambda_pos^2 * s^2 - lambda_pos * (x - mu) - s * lambda_pos * dlogpnorm.pos + s^2 * lambda_pos^2 * d2logpnorm.pos
      hess[i, i] <- - sum(gamma_pos * (1 - gamma_pos) * B_pos^2 + gamma_pos * dB_pos_dbeta)
      
      j <- i + 1
      if (!fix_beta_neg) {
        hess[i, j] <- hess[j, i] <- - sum(-gamma_pos * gamma_neg * B_pos * B_neg)
        j <- j + 1
      }
      if (!fix_mu) {
        dB_pos_dmu <- lambda_pos + lambda_pos * d2logpnorm.pos
        hess[i, j] <- hess[j, i] <- - sum(gamma_pos * (D_pos - D_bar) * B_pos + gamma_pos * dB_pos_dmu)
      }
      i <- i + 1
    }
    
    if (!fix_beta_neg) {
      dB_neg_dbeta <- 2 * lambda_neg^2 * s^2 + lambda_neg * (x - mu) - s * lambda_neg * dlogpnorm.neg + s^2 * lambda_neg^2 * d2logpnorm.neg
      hess[i, i] <- - sum(gamma_neg * (1 - gamma_neg) * B_neg^2 + gamma_neg * dB_neg_dbeta)
      
      j <- i + 1
      if (!fix_mu) {
        dB_neg_dmu <- -lambda_neg - lambda_neg * d2logpnorm.neg
        hess[i, j] <- hess[j, i] <- - sum(gamma_neg * (D_neg - D_bar) * B_neg + gamma_neg * dB_neg_dmu)
      }
      i <- i + 1
    }
    
    if (!fix_mu) {
      dD_0_dmu   <- -1 / s^2
      dD_pos_dmu <- d2logpnorm.pos / s^2
      dD_neg_dmu <- d2logpnorm.neg / s^2
      
      hess[i, i] <- - sum(gamma_0 * dD_0_dmu + gamma_pos * dD_pos_dmu + gamma_neg * dD_neg_dmu + 
                            gamma_0 * D_0^2 + gamma_pos * D_pos^2 + gamma_neg * D_neg^2 - D_bar^2)
    }
    
    attr(nllik, "hessian") <- hess
  }

  # The blocks above differentiate with respect to the constrained tail's
  # log-rate, but the free parameter is delta. With beta = B + softplus(delta),
  # B constant, and softplus' = plogis:
  #
  #   dnllik/ddelta   = (dnllik/dbeta) * plogis(delta)
  #   d2nllik/ddelta2 = (d2nllik/dbeta2) * plogis(delta)^2
  #                       + (dnllik/dbeta) * plogis(delta) * (1 - plogis(delta))
  #
  # Off-diagonals in delta's row and column take a single factor and no
  # curvature term, since B involves no other free parameter. Scaling the row
  # and then the column gives the diagonal its square automatically.
  if (!is.null(ratio_ref) && (calc_grad || calc_hess)) {
    # delta's slot in the parameter vector, then its position among the free
    # parameters only.
    if (fixed_comp == "none") {
      slot <- if (identical(ratio_ref, "pos")) 4 else 3
    } else {
      slot <- if (identical(ratio_ref, "pos")) 3 else 2
    }
    idx <- sum(!fix_par[seq_len(slot)])

    dbeta        <- plogis(delta)
    dnllik.dbeta <- grad[idx]

    grad[idx] <- dnllik.dbeta * dbeta
    attr(nllik, "gradient") <- grad

    if (calc_hess) {
      hess[idx, ] <- hess[idx, ] * dbeta
      hess[, idx] <- hess[, idx] * dbeta
      hess[idx, idx] <- hess[idx, idx] + dnllik.dbeta * dbeta * (1 - dbeta)
      attr(nllik, "hessian") <- hess
    }
  }

  return(nllik)
}
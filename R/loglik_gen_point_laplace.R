# Functions to compute the log likelihood under the generalized-point-Laplace
# prior. Currently used only for testing/reference.

loglik_gen_point_laplace = function(x, s, pi_0, pi_plus, pi_neg, l_plus, l_neg, mu) {
  return(sum(vloglik_gen_point_laplace(x, s, pi_0, pi_plus, pi_neg, l_plus, l_neg, mu)))
}

# Return log((1 - w)f + wg) as a vector (deal with cases of all slab and all spike separately for stability).
#
#' @importFrom stats dnorm
#'
vloglik_gen_point_laplace = function(x, s, pi_0, pi_plus, pi_neg, l_plus, l_neg, mu) {
  lf <- dnorm(x - mu, sd = s, log = TRUE)
  if (pi_plus <= 0 && pi_neg <= 0)  {
    return(lf) #pi_0 = 1 so don't need to multiply by it
  }
  
  lg_positive <- logg_gen_laplace(x - mu, s, l_plus, minus = FALSE)
  lg_negative <- logg_gen_laplace(x - mu, s, l_neg, minus = TRUE)
  # No point mass: the prior is the two slabs alone.
  if (pi_0 <= 0) {
    term_p <- log(pi_plus) + lg_positive
    term_n <- log(pi_neg) + lg_negative #log sum exp trick
    
    c_max <- pmax(term_p, term_n)
    return(c_max + log(exp(term_p-c_max) + exp(term_n - c_max)))
  }
  #else we have all 3
  lf <- dnorm(x - mu, sd = s, log = TRUE)
  term_0 <- log(pi_0) + lf
  term_p <- log(pi_plus) + lg_positive
  term_n <- log(pi_neg) + lg_negative #log sum exp trick 
  
  c_max <- pmax(term_0, pmax(term_p, term_n))
  return(c_max + log(exp(term_0-c_max) + exp(term_p-c_max) +exp(term_n-c_max)))
}

# This is the log of g: a generalized (asymmetric) Laplace convolved with a
#   normal. Argument x is assumed to be centered (that is, x - mu). The
#   nonpositive component has rate lambda_pos and weight pi_plus; the nonnegative  component has rate lambda_neg and weight pi_neg.
#  
#
#' @importFrom stats pnorm
#'
logg_gen_laplace = function(x, s, lambda, minus = FALSE) {
  #x here is already x - mu 
  # The parameter 'a' represents either lambda_plus or lambda_minus.
  if (!minus) {
    # POSITIVE SLAB
    # Mathematically should be 
    #log(lambda_+) - lambda_+(x_i - mu) + (lambda_+^2 * s_i^2)/2 +     
    #log(Phi((x_i - mu - lambda_+ * s_i^2)/s_i))
    
    lpnorm_pos <- pnorm((x - s^2 * lambda) / s, log.p = TRUE) #log.p stably computes 
    lg <- log(lambda) - (lambda * x) + (0.5 * s^2 * lambda^2) + lpnorm_pos
    
  } else {
    # NEGATIVE SLAB
    # Math: log(lambda_-) + lambda_-(x_i - mu) + (lambda_-^2 * s_i^2)/2 + log(Phi((mu - x_i - lambda_- * s_i^2)/s_i))
    # Again, (mu - x_i) is just -x here (relevant for Phi term)
    # We use lower.tail = FALSE on the positive equivalent for numerical  stability, since Phi(-z) = 1 - Phi(z)
    lpnorm_neg <- pnorm((x + s^2 * lambda) / s, log.p = TRUE, lower.tail = FALSE)
    lg <- log(lambda) + (lambda * x) + (0.5 * s^2 * lambda^2) + lpnorm_neg
  }
  return(lg)
}

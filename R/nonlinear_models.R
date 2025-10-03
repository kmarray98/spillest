#' Rescale nonlinear estimators under assumption 5
#' 
#' This function computes bias-corrected estimates of spillover effects from non-linear models with
#' sampled network data under the assumption that the distribution of sampled
#'and unsampled links are independent. As instruments, use neighbours of up to degree 2.
#' @param y Dependent variable in the regression.
#' @param x Vector of shocks.
#' @param z Matrix of covariates. If you want to include an intercept, you need to add column of ones here.
#' @param H Adjacency matrix of sampled network.
#' @param db Mean missing degree for nodes with at least one incorrectly sampled link. 
#' @return Vector of coefficient estimates, and vector of standard errors (computed using
#' normal lm procedure) assuming no uncertainty in estimate of db. To get coefficients
#' assuming some uncertainty in db, you need to use the bootstrap function.
#' @export
rescaled_est_5 <- function(y,x,z, H, db, nb){
  z1 <- x
  z2 <- t(H) %*% H %*% x
  z3 <- db * H %*% x
  insts <- matrix(c(z1, z2, z3), ncol = 3)
  Hy <- H %*% y
  spill_hat <- insts %*% solve(t(insts) %*% insts) %*% t(insts) %*% Hy
  out_sum <- summary(lm(y~0+spill_hat + x)) # use lm object here to get the full output
  div <- 1 + solve(t(spill_hat) %*% spill_hat) %*% (t(spill_hat) * db) %*% y
  b <- out_sum$coefficients[1]#solve(t(spill) %*% spill) %*% (t(spill) %*% y) # first lets build the standard ols estimates
  # now lets fill back in the coefficients, resulting standard error, and test statistics
  out_sum$coefficients[1] <- b/div
  num_covariates <- 1
  out_sum$coefficients[num_covariates+2] <- (1/(div)) * out_sum$coefficients[num_covariates+2]
  # computing se assuming no uncertainty in eta
  return(out_sum$coefficients[1:(num_covariates+1)], out_sum$coefficients[(num_covariates+2):(2*num_covariates + 2)])
}
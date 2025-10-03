# functions to rescale linear estimates

# example to test this
#rm(list=ls())
#x <- rnorm(100, mean=1, sd=1)
#H <- matrix(rbinom(100*100, 1, 0.05), ncol=100)
#diag(H) <- 0

#B <- matrix(rbinom(100*100, 1, 0.01), ncol=100)
#diag(B) <- 0
#B_sums <- rowSums(B)
#db <- mean(B_sums[B_sums > 0])
#nb <- length(B_sums[B_sums > 0])
#G <- H + B
#G[G > 1] <- 1

#z <- matrix(rnorm(200, 3,2.5), ncol=2)# some arbitrary covariates
#y <-  G %*% x * 0.5 + z %*% c(-0.3, 2.0) + rnorm(100, mean=0, sd=1)


#' Rescale linear estimates when sampling scheme satisfies assumption 4a
#'
#' This function computes bias-corrected OLS estimates of spillover effects from
#' sampled network data under the assumption that the distribution of sampled
#' degrees and unsampled degrees are independent for node with at least one 
#' incorrectly sampled link. In cases where they are dependent, see rescaled_est_4b.
#' @param y Dependent variable in the regression.
#' @param x Vector of shocks.
#' @param z Matrix of covariates. If you want to include an intercept, you need to add column of ones here.
#' @param H Adjacency matrix of sampled network.
#' @param db Mean missing degree for nodes with at least one incorrectly sampled link. 
#' @param nb Number of nodes with at least one incorrectly sampled link.
#' @return Vector of coefficient estimates, and vector of standard errors (computed using
#' normal lm procedure) assuming no uncertainty in estimate of db. To get coefficients
#' assuming some uncertainty in db, you need to use the bootstrap function.
#' @export
rescaled_est_4a <- function(y,x,z, H, db, nb){
  spill <- H %*% x
  out_sum <- summary(lm(y~0+spill + z)) # use lm object here to get the full output
  b <- out_sum$coefficients[1]#solve(t(spill) %*% spill) %*% (t(spill) %*% y) # first lets build the standard ols estimates
  dh <- mean(rowSums(H)) 
  n <- length(y)
  x_bar <- mean(x)
  eta <- ((nb/n) * dh* db*  x_bar^2) / mean(spill^2)
  # now lets fill back in the coefficients, resulting standard error, and test statistics
  out_sum$coefficients[1] <- b/(1+eta)
  num_covariates <- dim(z)[2]
  out_sum$coefficients[num_covariates+2] <- (1/(1+eta)) * out_sum$coefficients[num_covariates+2]
  # computing se assuming no uncertainty in eta
  return(out_sum$coefficients[1:(num_covariates+1)], out_sum$coefficients[(num_covariates+2):(2*num_covariates + 2)])
}





  


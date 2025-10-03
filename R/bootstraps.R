
# need to use -1 for undersampled
one_bootstrap_est_un <- function(out_val, y,x,z,H, Mask_mat, num_missing){
  n = length(y)
  is_masked <- rowSums(Mask_mat)
  one_positions <- which(Mask_mat == 1, arr.ind = TRUE)
  sampled_positions <- one_positions[sample(nrow(one_positions), num_missing), , drop=FALSE]
  
  new_mat <- matrix(0, nrow = nrow(H), ncol = ncol(H))
  new_mat[sampled_positions] <- 1
  spill = H %*% x
  out_sum <- summary(lm(y~0+spill + z))
  new_spill <- new_mat %*% x
  eta_hat <- ((nb/n) * mean(new_spill[is_masked > 0])* dh*x_bar) / mean(spill^2)
  out_val = out_sum$coefficients[1]/(1  + eta_hat)            
}

one_bootstrap_est_ov <- function(out_val, y,x,z,H, Mask_mat, num_missing){
  n = length(y)
  is_masked <- rowSums(Mask_mat)
  one_positions <- which(Mask_mat == 1, arr.ind = TRUE)
  sampled_positions <- one_positions[sample(nrow(one_positions), num_missing), , drop=FALSE]
  
  new_mat <- matrix(0, nrow = nrow(H), ncol = ncol(H))
  new_mat[sampled_positions] <- -1
  spill = H %*% x
  out_sum <- summary(lm(y~0+spill + z))
  new_spill <- new_mat %*% x
  eta_hat <- ((nb/n) * mean(new_spill[is_masked > 0])* dh*x_bar) / mean(spill^2)
  out_val = out_sum$coefficients[1]/(1  + eta_hat)            
}

#' Bootstrap estimators for uncertainty of bias-corrected estimators for spillover
#' effects.
#'
#' This function computes bootstrap estimates of uncertainty for bias-corrected
#' estimators of spillover effects for linear regression estimators when network 
#' is undersampled.
#' @param y Dependent variable in the regression.
#' @param x Vector of shocks.
#' @param z Matrix of covariates. If you want to include an intercept, you need to add column of ones here.
#' @param H Adjacency matrix of sampled network.
#' @param mask_mat Binary matrix of positions of possibly incorrectly sampled links. 
#' @param nos Number of bootstrap estimates.
#' @return Standard error of estimates.
#' @export
uniform_boostrap_un <- function(y,x,z, H, Mask_mat, num_missing, nos){
  sm_lst <- rep(c(0),times=nos)
  out_lst <- unlist(lapply(sm_lst, FUN=one_bootstrap_est_un, y=y, x=x, z=z, H=H, Mask_mat=Mask_mat, num_missing=num_missing))
  return(sd(out_lst))
}


#' Bootstrap estimators for uncertainty of bias-corrected estimators for spillover
#' effects.
#'
#' This function computes bootstrap estimates of uncertainty for bias-corrected
#' estimators of spillover effects for linear regression estimators when network
#' is oversampled.
#' @param y Dependent variable in the regression.
#' @param x Vector of shocks.
#' @param z Matrix of covariates. If you want to include an intercept, you need to add column of ones here.
#' @param H Adjacency matrix of sampled network.
#' @param Mask_mat Binary matrix of positions of possibly incorrectly sampled links. 
#' @param num_missing Number of missing links
#' @param nos Number of bootstrap estimates.
#' @return Standard error of estimates.
#' @export
uniform_boostrap_ov <- function(y,x,z, H, Mask_mat, num_missing, nos){
  sm_lst <- rep(c(0),times=nos)
  out_lst <- unlist(lapply(sm_lst, FUN=one_bootstrap_est_ov, y=y, x=x, z=z, H=H, Mask_mat=Mask_mat, num_missing=num_missing))
  return(sd(out_lst))
}
  
  
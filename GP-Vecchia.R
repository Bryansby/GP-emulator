##################################### GP Vecchia ######################################
GPvecchia <- function(x, Y, kernel_function = 'sexp', scale, ls, nugget, x_star, n){
  #' Function using Vecchia approximation in the construction of GP emulator
  #' This Vecchia approximation do condition only on the observation points
  #' @param x is the input (size of m * d)
  #' @param Y is the corresponding output (size of 1 * d)
  #' @param kernel_function specifies the kernel type to be used
  #' @param scale is the scale parameter
  #' @param ls is the length-scale parameter
  #' @param x_star is the points to be predicted or tested(size of p * d)
  #' @param n is the number of neighbors to return
  #' @export result the list contains mean and variance at the x_star
  
  sexp <- function(x, nugget, ls) {
    #' The function to get the correlation matrix when using the squared-exponential kernel
    #' @param x is the input(size = M * D)
    #' @param nugget is the nugget
    #' @param ls is the length-scale parameter, which have size of 1*d
    #' @export R The correlation matrix R of size M * M
    R <- exp(-as.matrix(dist(x))^2 / ls^2)
    diag(R) <- diag(R) + nugget
    return(R)
  }
  # get nearest neighbor array
  find_nn <- function(obs_locs, pred_locs, n){
    #' Function find the five nearest observations locations to the prediction locations
    #' @param obs_locs matrix of observation locations (size: m * d)
    #' @param pred_locs the prediction locations (size: p * d)
    #' @param n is the number of locations nearest the prediction location
    #' @export NNarray each row contains the nearest n locations index in the obs_locs matrix corresponding to the row index in pred_loc matrix
    NNarray <- matrix(0, nrow = nrow(pred_locs), ncol = n)
    for (i in 1:nrow(pred_locs)){
      prediction_location <- pred_locs[i,]
      # Calculate the Euclidean distances
      distances <- apply(obs_locs, 1, function(row) sqrt(sum((row - prediction_location)^2)))
      # Get the indices of the 5 smallest distances
      nearest_indices <- c(order(distances)[1:n])
      NNarray[i, ] <- nearest_indices
    }
    return(NNarray)
  }
  
  NNarray <- find_nn(as.matrix(x), as.matrix(test_x), n) # Get the nearest test locations index
  m <- nrow(x) # Get the number of designs
  d <- ncol(x) # Get the dimension of the inputs
  t <- nrow(as.matrix(x_star)) # Get the number of test/prediction inputs/locations
  
  mean <- c(1:t)
  variance <- c(1:t)
  
  for (j in 1:t) {
    xj <- as.matrix(x[sort(NNarray[j,]), ])
    Yj <- as.matrix(Y[sort(NNarray[j,]), ])
    # The r(x_star) matrix, which have size M * 1
    r <- matrix(0, nrow = m, ncol = 1)
    for (i in 1: n){
      r[i] <- exp(-sum((xj[i,] - x_star[j,])^2 / (ls^2)) )
    }
    
    # Get the correlation matrix R
    R <- as.matrix(sexp(xj, nugget, ls))
    
    L <- t(chol(R))
    
    mean[j] <- t(forwardsolve(L, r)) %*% forwardsolve(L, Yj)
    variance[j] <- scale * (1 + nugget - crossprod(forwardsolve(L, r)))
  }
  result <- list(
    kernel = "sexp",
    mean = mean,
    sd = sqrt(variance)
  )
  
  return(result)
  
}

####################### GP Vecchia (condition only on the observations) ################################
GPvecchia_obs <- function(x, Y, kernel_function = 'sexp', scale, ls, nugget, x_star, n){
  #' @description Function using Vecchia approximation in the construction of GP emulator
  #' @description This Vecchia approximation do condition only on the observation points
  #' 
  #' @param x is the input (size of m * d)
  #' @param Y is the corresponding output (size of 1 * d)
  #' @param kernel_function specifies the kernel type to be used
  #' @param scale is the scale parameter
  #' @param ls is the length-scale parameter
  #' @param x_star is the points to be predicted or tested(size of p * d)
  #' @param n is the number of neighbors to return
  #' 
  #' @export result the list contains mean and variance at the x_star
  
  kse <- function(distance, ls) {
    #' @description The function to do squared-exponential kernel function
    #' @param distance is the distance
    #' @param ls is the length-scale parameter
    #' @return the value of squared exponential kernel(not include exponential)
    return(- distance ^ 2 / (ls ^ 2))
  }
  
  sexp <- function(x, nugget, ls) {
    #' @description The function to get the correlation matrix when using the squared-exponential kernel
    #' @param x is the input(size = M * D)
    #' @param ls is the length-scale parameter, which have size of 1*d
    #' @param m is the number of design points(M)
    #' @return The correlation matrix R of size M * M
    
    D <- ncol(x) # dimension
    M <- nrow(x) # number of design
    
    R <- matrix(NA, M, M)
    
    for (j in 1:M) {
      for (i in 1:M) {
        # Consider (j, i) entry of the covariance matrix
        first <- NULL
        second <- NULL
        for (d in 1:D) {
          first[d] <- kse(x[i, d] - x[j, d], ls = ls[d])
          second[d] <- ifelse(x[i, d] == x[j, d], 1, 0)
        }
        R[j, i] <- exp(sum(first)) + nugget * prod(second)
      }
    }
    
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


############################# GP Vecchia (Observation & Prediction) ##############################
GPvecchia <- function(x, Y, kernel_function = 'sexp', scale, ls, nugget, x_star, n){
  #' @description Function using Vecchia approximation in the construction of GP emulator
  #' @description This Vecchia approximation do condition only on the observation points
  #' @param x is the input (size of m * d)
  #' @param Y is the corresponding output (size of 1 * d)
  #' @param kernel_function specifies the kernel type to be used
  #' @param scale is the scale parameter
  #' @param ls is the length-scale parameter
  #' @param x_star is the points to be predicted or tested(size of p * d)
  #' @param n is the number of neighbors to return
  #' @export result the list contains mean and variance at the x_star
  
  kse <- function(distance, ls) {
    #' @description The function to do squared-exponential kernel function
    #' @param distance is the distance
    #' @param ls is the length-scale parameter
    #' @return the value of squared exponential kernel(not include exponential)
    return(- distance ^ 2 / (ls ^ 2))
  }
  
  sexp <- function(x, nugget, ls) {
    #' @description The function to get the correlation matrix when using the squared-exponential kernel
    #' @param x is the input(size = M * D)
    #' @param ls is the length-scale parameter, which have size of 1*d
    #' @param m is the number of design points(M)
    #' @return The correlation matrix R of size M * M
    
    D <- ncol(x) # dimension
    M <- nrow(x) # number of design
    
    R <- matrix(NA, M, M)
    
    for (j in 1:M) {
      for (i in 1:M) {
        # Consider (j, i) entry of the covariance matrix
        first <- NULL
        second <- NULL
        for (d in 1:D) {
          first[d] <- kse(x[i, d] - x[j, d], ls = ls[d])
          second[d] <- ifelse(x[i, d] == x[j, d], 1, 0)
        }
        R[j, i] <- exp(sum(first)) + nugget * prod(second)
      }
    }
    
    return(R)
  }
  
  find_nn <- function(obs_locs, pred_locs, n) {
    #' @description Function to find the nearest locations(both the observation locations and prediction locations) to the prediction locations
    #' @description And here we consider the observation locations first, and prediction locations after
    #' @param obs_locs matrix of observation locations (size: m * d)
    #' @param pred_locs the prediction locations (size: p * d)
    #' @param n the number of nearest locations to find
    #' @return NNarray a matrix where each row contains the indices of the nearest n locations in obs_locs for each row in pred_locs
    
    all_locs <- as.matrix(rbind(obs_locs, pred_locs)) # make the observation locations ordered before the prediction locations
    NNarray <- matrix(NA, nrow = nrow(pred_locs), ncol = n) # to store the index of the nearest n locations
    n_obs <- nrow(obs_locs)
    
    for (i in 1:nrow(pred_locs)) {
      prediction_location <- all_locs[n_obs + i, ]
      # Calculate the Euclidean distances
      distances <- sqrt(  rowSums( as.matrix( (all_locs[1:(n_obs + i - 1), ] - prediction_location) ^ 2, ncol = ncol(obs_locs))   ) )
      # Get the indices of the n smallest distances
      nearest_indices <- order(distances)[1:n]
      NNarray[i, ] <- nearest_indices
    }
    
    return(NNarray)
  }
  
  # Function to check if a row in matrix2 exists in matrix1
  row_in_matrix1 <- function(row, matrix1) {
    apply(matrix1, 1, function(x) all(x == row))
  }
  # Logical vector indicating rows to keep
  keep_rows <- !apply(as.matrix(x_star), 1, function(row) any(row_in_matrix1(row, as.matrix(x))))
  # Subset matrix2 to keep only unique rows
  x_star <- as.matrix(x_star)[keep_rows, , drop = FALSE]
  
  x_star <- as.matrix(x_star[sample(nrow(x_star)), ]) # Random ordering the prediction locations
  
  NNarray <- find_nn(as.matrix(x), as.matrix(x_star), n) # Get the nearest test locations index in all locations, which is obs first, pred after
  m <- nrow(x) # Get the number of designs
  d <- ncol(x) # Get the dimension of the inputs
  t <- nrow(as.matrix(x_star)) # Get the number of test/prediction inputs/locations
  all_locs <- rbind(as.matrix(x, nrow = m), as.matrix(x_star, nrow = t)) # combine the observed location and prediction location together
  all_value <- rbind(as.matrix(Y, nrow = m), matrix(NA, nrow = t, ncol = 1))
  variance <- c(1:t)
  
  for (j in 1:t) {
    xj <- as.matrix(all_locs[sort(NNarray[j,]), ])
    Yj <- as.matrix(all_value[sort(NNarray[j,]), ])
    # The r(x_star) matrix, which have size M * 1
    r <- matrix(0, nrow = n, ncol = 1)
    for (i in 1: n){
      r[i] <- exp(-sum((xj[i,] - x_star[j,])^2 / (ls^2)) )
    }
    
    # Get the correlation matrix R
    R <- as.matrix(sexp(xj, nugget, ls))
    
    L <- t(chol(R))
    
    all_value[m + j,] <- t(forwardsolve(L, r)) %*% forwardsolve(L, Yj)
    variance[j] <- scale * (1 + nugget - crossprod(forwardsolve(L, r)))
  }
  result <- list(
    x = all_locs[(m+1):nrow(all_locs), ],
    kernel = "sexp",
    mean = all_value[(m+1):nrow(all_value),],
    sd = sqrt(variance)
  )
  
  return(result)
  
}


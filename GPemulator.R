#################################### Training of GP ##########################################
GPtrainingmle <- function(x, Y, kernel_function = 'sexp', optim_method = "SANN"){
  #' Function to train the GP using MLE
  #'
  #' @param x is the M * D input
  #' @param kernel_function's default is 'sexp'
  #' @param optim_method's default is SANN
  #' @return return three parameters estimation: scale, length-scale, and nugget
  
  m <- nrow(x)
  
  sexp <- function(x, nugget, ls, m) {
    ### The function to get the correlation matrix when using the squared-exponential kernel
    ##Input:
    # x is the input(size = M * D)
    # ls is the length-scale
    # m is the number of design points(M)
    ## Output:
    # The correlation matrix R of size M * M
    kse <- function(d, ls){
      ## The function to do squared-exponential kernel function
      # d is the distance
      # ls is the length-scale parameter
      return(sum(-d^2 / (ls^2)))
    }
    R <- matrix(0, nrow = m, ncol = m)
    for (j in 1:m){
      for (i in 1:m){
        R[i, j] <- exp(kse((x[i,] - x[j,]), ls))
      }
    }
    diag(R) <- diag(R) + nugget
    return(R)
  }
  
  nugget <- 1e-8 # nugget
  log_likelihood <- function(params, x, Y){
    l <- exp(params[1]) # length-scale parameter
    nugget <- 1e-8 # nugget
    
    R <- sexp(x, nugget, l, m)
    
    # Compute the covariance matrix
    m <- nrow(x) # Get the number of designs
    
    # Compute the log-likelihood
    sc <- (1/m) * t(Y) %*% solve(R, Y)
    #log_likelihood_value <- - 0.5 * (log(det(R)) + m * log(2 * pi * sc) + m)
                            
    log_likelihood_value <- 0.5 * (log(det(R) + m * log(sc)))
    
    return(log_likelihood_value)
  }
  
  # Initial guess for parameters
  initial_params <- c(log(1))
  
  # Perform optimization to maximize the log-likelihood
  optimized_params <- optim(par = initial_params, fn = log_likelihood, x = x, Y = Y, method = "SANN")
  
  # Extract the optimized parameters
  l_opt <- exp(optimized_params$par[1])
  R1 <- sexp(x, nugget, l_opt, m)
  sc <- (1/m) * t(Y) %*% solve(R1, Y)
  
  result <- list(
    kernel = "sexp",
    lengthscales = l_opt,
    scale = sc,
    nugget = nugget
  )
  
  return(result)

################################## GP training using reference prior ######################################
GPtrainingrobust <- function(x, Y, zero_mean = 'Yes', kernel = 'pow_exp', alpha = 2){
  #' PACKAGE rgasp IS REQUIRED
  #' Function to train the GP emulator using rgasp
  #'
  #' @param x is the M*D input
  #' @param Y is the corresponding output which have size 1*D
  #' @param zero_mean determines whether the trend zero or not
  #' @param kernel sepecifies the kernel function
  #' @param alpha determines the roughness parameters when using pow_exp kernel
  #' @return the result list contains kernel type, lengthscales(1*D), scale, nugget
  
  m <- rgasp(x, Y, zero.mean = zero_mean, 
             kernel_type = kernel, alpha = alpha)
  sc <- m@sigma2_hat
  l_opt <- 1/m@beta_hat
  
  result <- list(
    kernel = "sexp",
    lengthscales = l_opt,
    scale = sc,
    nugget = 1e-8
  )
  
  return(result)
}

###################################### GP emulator ###################################################
GPemulator <- function(x, Y, kernel_function = 'sexp', scale, ls, nugget, x_star){
  #' @description Function to construct the Gaussian Process emulator
  #' 
  #' @param x: M-point design of D-dimensional inputs to a computer model, which size is M * D
  #' @param Y: the corresponding scalar-valued out with size of M * 1
  #' @param kernel_function: can be selected between Matérn2.5(mat2.5) and squared exponential(sexp)
  #' @param scale: the scale parameter
  #' @param nugget: the nugget added in the covariance matrix diagonal
  #' @param ls: the length-scales parameters in the kernel function
  #' @param x_star: new input positions/locations which have size of T * D
  #'
  #' @return: 
  #' The mean at new input positions x_star
  #' The std at new input positions x_star
  
  m <- nrow(x) # Get the number of designs
  d <- ncol(x) # Get the dimension of the inputs
  n_pred <- nrow(x_star) # Get the number of test/prediction inputs/locations
  
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
  
  mean <- c(1:n_pred)
  variance <- c(1:n_pred)
  
  for (j in 1:n_pred) {
    # The r(x_star) matrix, which have size M * 1
    r <- matrix(0, nrow = m, ncol = 1)
    for (i in 1: m){
      r[i] <- exp(-sum((x[i,] - x_star[j,])^2 / (ls^2)) )
    }
    # Get the correlation matrix R
    R <- as.matrix(sexp(x, nugget, ls))
    L <- t(chol(R))
    mean[j] <- t(forwardsolve(L, r)) %*% forwardsolve(L, Y)
    variance[j] <- scale * (1 + nugget - crossprod(forwardsolve(L, r)))
  }
  
  result <- list(
    kernel = "sexp",
    mean = mean,
    sd = sqrt(variance)
  )
  
  return(result)
  
}

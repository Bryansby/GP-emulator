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
  #' Function to construct the Gaussian Process emulator
  #'
  #' @param x: M-point design of D-dimensional inputs to a computer model, which size is M * D
  #' @param Y: the corresponding scalar-valued out with size of M * 1
  #' @param kernel_function: can be selected between Matérn2.5(mat2.5) and squared exponential(sexp)
  #' @param scale: the scale parameter
  #' @param nugget: the nugget added in the covariance matrix diagonal
  #' @param ls: the length-scales parameter in the kernel function
  #' @param x_star: a new input position which have size of 1 * D
  #'
  #' @return: 
  #' The mean at the new input position x_star
  #' The variance at the new input position x_star
  
  m <- nrow(x) # Get the number of designs
  
  sexp <- function(x, nugget, ls) {
    ### The function to get the correlation matrix when using the squared-exponential kernel
    ##Input:
    # x is the input(size = M * D)
    # ls is the length-scale parameter, which have size of 1*d
    # m is the number of design points(M)
    ## Output:
    # The correlation matrix R of size M * M
    R <- exp(-as.matrix(dist(x))^2 / ls^2)
    diag(R) <- diag(R) + nugget
    return(R)
  }
  
  # The r(x_star) matrix, which have size M * 1
  r <- matrix(0, nrow = m, ncol = 1)
  for (i in 1: m){
    r[i] <- exp(-sum((x[i,] - x_star)^2 / (ls^2)) )
  }
  
  # Get the correlation matrix R
  R <- as.matrix(sexp(x, nugget, ls))
  
  L <- t(chol(R))
  mean <- t(forwardsolve(L, r)) %*% forwardsolve(L, Y)
  #mean <- t(r) %*% solve(R, Y)

  variance <- scale * (1 + nugget - crossprod(forwardsolve(L, r)))
  #variance <- scale * (1 + nugget - t(r) %*% solve(R, r))
  
  result <- list(
    kernel = "sexp",
    mean = mean,
    sd = sqrt(variance)
  )
  
  return(result)
}

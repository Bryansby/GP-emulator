GPtraining <- function(x, Y, kernel_function = 'sexp', optim_method = "L-BFGS-B"){
  # Function to train the GP
  # x is the M * D input
  # kernel_function's default is 'sexp'
  # optim_method's default is L-BFGS-B
  # return three parameters estimation: scale, length-scale, and nugget
  log_likelihood <- function(params, x, Y){
    scale <- exp(params[1]) # scale parameter
    ls <- exp(params[2]) # length-scale parameter
    nugget <- exp(params[3]) # nugget
    
    # Compute the covariance matrix
    m <- nrow(x) # Get the number of designs
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
        return(exp(-d^2 / (2 * ls^2)))
      }
      R <- matrix(0, nrow = m, ncol = m)
      for (j in 1:m){
        for (i in 1:m){
          R[i, j] <- prod(kse(abs(x[i,] - x[j,]), ls))
        }
      }
      diag(R) <- diag(R) + nugget
      return(R)
    }
    
    K <- scale * sexp(x, nugget, ls, m) # Get the covariance matrix
    
    # Compute the log-likelihood
    log_likelihood_value <- 0.5 * (t(Y) %*% solve(K) %*% Y + log(det(K)) + m * log(2 * pi))
    
    return(log_likelihood_value)
    
  }
  
  # Initial guess for parameters
  initial_params <- c(log(1), log(1), log(.1))
  
  # Perform optimization to maximize the log-likelihood
  optimized_params <- optim(par = initial_params, fn = log_likelihood, x = x, Y = Y, method = optim_method)
  
  # Extract the optimized parameters
  sigma_sq_opt <- exp(optimized_params$par[1])
  l_opt <- exp(optimized_params$par[2])
  nugget_opt <- exp(optimized_params$par[3])
  
  return(c(sigma_sq_opt, l_opt, nugget_opt))
}


GPtraining(x, Y)

############################################################################################

GPemulator <- function(x, Y, kernel_function = 'sexp', scale, ls, nugget, x_star){
  ### Function to construct the Gaussian Process emulator
  ## Inputs:
  # x: M-point design of D-dimensional inputs to a computer model, which size is M * D
  # Y: the corresponding scalar-valued out with size of M * 1
  # kernel_function: can be selected between Matérn2.5(mat2.5) and squared exponential(sexp)
  # scale: the scale parameter
  # nugget: the nugget added in the covariance matrix diagonal
  # ls: the length-scales parameter in the kernel function
  # x_star: a new input position which have size of 1 * D
  ## Outputs:
  # The mean at the new input position x_star
  # The variance at the new input position x_star
  
  m <- nrow(x) # Get the number of designs
  
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
      return(exp(-d^2 / (2 * ls^2)))
    }
    R <- matrix(0, nrow = m, ncol = m)
    for (j in 1:m){
      for (i in 1:m){
        R[i, j] <- prod(kse(abs(x[i,] - x[j,]), ls))
      }
    }
    diag(R) <- diag(R) + nugget
    return(R)
  }
  
  # The r(x_star) matrix, which have size M * 1
  r <- matrix(0, nrow = m, ncol = 1)
  for (i in 1: m){
    r[i] <- prod(exp(-abs(x[i,] - x_star)^2 / (2 * ls^2))) + as.numeric(x[i,] == x_star) * nugget
  }
  
  # Get the correlation matrix R
  R <- as.matrix(sexp(x, nugget, ls, m))
  
  mean <- t(r) %*% solve(R) %*% Y
  variance <- sqrt(scale) * (1 + nugget - t(r) %*% solve(R) %*% r)
  
  return(c(mean, sqrt(variance)))
}

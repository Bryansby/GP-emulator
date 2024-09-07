library(plgp)
library(MASS)
library(gridExtra)
library(ggplot2)
library(deepgp)
library(mvtnorm)
################## Training #########################
fit_two_layer3.0_matern <- function(x, Y, u = 2, l = 1, ls_y = 1, ls_w = 1, node, W, n_iteration = 10000, burn_in = 7000, nugget = 1e-6, v){
  #' @description To do the training for simple two layer Deep Gaussian Processes model
  #' @param x is the input value, which have the size of M * D
  #' @param Y is the output value, which have the size of M * 1
  #' @param u is the parameter for the proposals distribution, default to 2 (according to Sauer)
  #' @param l is the parameter for the proposals distribution, default to 1 (according to Sauer)
  #' @param ls_y is the initial value for the length-scale parameter in the first layer, default to 1
  #' @param ls_w is the initial value for the length-scale parameter in the second layer, default to 1
  #' @param node is the number of latent nodes in the model
  #' @param W is the initial output value come from the second layer, is the latent variable, default to 1
  #' @param n_iteration is the number of iteration for the MCMC, default is 10000
  #' @param burn_in the iteration number for MCMC to warm up
  #' @param nugget is the nugget, set to 1e-4
  #' @param v is the smooth parameter for Matern Kernel
  #' @returns result list contain:
  #'                              1. result summary
  #'                              2. samples of theta_y
  #'                              3. samples of theta_w
  #'                              4. samples of latent variable w
  #'                              5. trace plot of theta_y
  
  library(plgp)
  library(MASS)
  library(gridExtra)
  library(ggplot2)
  library(deepgp)
  library(mvtnorm)
  
  scale_sampling <- function(w, ls_y, Y, g = nugget, v){
    n <- nrow(Y)
    R <- deepgp:::Matern(distance(w), 1, ls_y, g, v)
    quadterm <- t(Y) %*% (deepgp:::invdet(R))$Mi %*% (Y)
    scale <- c(quadterm) / n
    return(scale)
  }
  
  loglik_yw <- function(Y, w, ls_y, g = 1e-6, v){
    #' @description To compute the log-likelihood of the first layer
    n <- nrow(Y)
    R <- deepgp:::Matern(distance(w), 1, ls_y, g, v)
    quadterm <- t(Y) %*% (deepgp:::invdet(R))$Mi %*% (Y)
    return(- 0.5 * n * log(2*pi*quadterm/n) - 0.5 * (log(det(R))))
  }
  
  loglik_wx <- function(w, x, ls_w, g = 1e-6, v){
    #' @description To compute the log-likelihood of the second layer
    w <- matrix(w, ncol = 1)
    R <- deepgp:::Matern(distance(x), 1, ls_w, g, v)
    quadterm <- t(w) %*% (deepgp:::invdet(R))$Mi %*% (w)
    log_l <- (- 0.5) * log(det(R)) - 0.5 * (quadterm)
    return(log_l)
  }
  
  MH_1 <- function(ls, Y, w, u = 2, l = 1, alpha = 1.5, beta = 0.65, v){
    #' @description Function to sample the length-scale parameter using Metropolis-Hasting in layer 1, ls_y
    w <- as.matrix(w, ncol = node)
    ls_star <- runif(1, min = l*ls / u, max = u*ls / l) # new value
    log_alpha <- loglik_yw(Y, w, ls_star, 1e-6, v) + 
      dgamma(ls_star - 1.490116e-08, alpha, beta, log = TRUE) +
      log(ls) -
      loglik_yw(Y, w, ls, 1e-6, v) - 
      dgamma(ls - 1.490116e-08, alpha, beta, log = TRUE) - 
      log(ls_star)
    U <- runif(1, 0, 1)
    if(log_alpha > log(U)){ # Accepted
      return(ls_star)
    }
    else{ # Rejected
      return(ls)
    }
  }
  
  MH_2 <- function(ls, w, x, u = 2, l = 1, alpha = 1.5, beta = 0.975, v){
    #' @description Function to sample the length-scale parameter using Metropolis-Hasting in layer 2, ls_w
    w <- as.matrix(w, ncol = 1)
    ls_star <- runif(1, min = l*ls / u, max = u*ls / l) # new value
    log_alpha <- loglik_wx(w, x, ls_star, 1e-6, v) + 
      dgamma(ls_star - 1.490116e-08, alpha, beta, log = TRUE) + 
      log(ls) - 
      loglik_wx(w, x, ls, 1e-6, v) - 
      dgamma(ls - 1.490116e-08, alpha, beta, log = TRUE) -
      log(ls_star)
    U <- runif(1, 0, 1)
    if(log_alpha > log(U)){ # Accepted
      return(ls_star)
    }
    else{ # Rejected
      return(ls)
    }
  }
  
  ESS_w <- function(x, Y, w, ls_y, ls_w, p = node, g = nugget, v){
    #' @description To do the Elliptical Slice Sampling update for latent variable w
    #' @param x is the input of the model
    #' @param Y is the Output value of the model
    #' @param w is the initial value for the latent variable
    #' @param ls_y is the length-scale parameter sampled from MH_1
    #' @param ls_w is the length-scale parameter sampled from MH_2
    #' @param p is the nodes in the latent layer
    #' @param g is the nugget
    #' @return w is the latent variable matrix(size m*p) after the ESS
    m <- nrow(w)
    if (p != ncol(w)){
      stop("Error: The initial value for the latent value doesn't match the nodes number!")
    }
    for (i in 1:p) {
      theta <- runif(1, 0, 2*pi) # angle
      theta_min <- theta - 2*pi  # lower basket
      theta_max <- theta         # upper basket
      nu <- mvtnorm::rmvnorm(1, mean = matrix(0, nrow = nrow(w)), sigma = deepgp:::Matern(distance(x), 1, ls_w[,i], 0, v)) # random draw from the prior, size = m*1
      w_prev <- w[, i]
      #ll_prev <- loglik_yw(Y, w[,i], ls_y, g, v)
      ll_prev <- loglik_yw(Y, w, ls_y, g, v)
      accept <- FALSE
      count <- 0
      while (accept == FALSE){
        count <- count + 1
        w[,i] <- w_prev * cos(theta) + nu * sin(theta) # Proposal
        dw <- deepgp:::sq_dist(w)
        #log_alpha <- loglik_yw(Y, w[,i], ls_y, g, v) - ll_prev # log-alpha
        log_alpha <- loglik_yw(Y, w, ls_y, g, v) - ll_prev # log-alpha
        U <- runif(1, 0, 1)
        # Check if the proposed sample is on the slice
        if (log_alpha > log(U)) # Accepted
        { 
          accept <- TRUE
        } 
        else { # Rejected
          # Shrink the bracket
          if (theta < 0) {
            theta_min <- theta
          } else {
            theta_max <- theta
          }
          # Draw a new angle from the updated bracket
          theta <- runif(1, theta_min, theta_max)
        }
      }
    }
    return(w)
  }
  
  cat("Training ... \n")
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_iteration, style = 3)
  
  theta_y_samples <- c(ls_y, rep(NA, n_iteration - 1))                              # To store theta_y samples (layer 1)
  theta_w_samples <- rbind(ls_w, matrix(NA, ncol = node, nrow = n_iteration - 1))   # To store theta_w samples (layer 2)
  w_samples <- vector("list", n_iteration)                                          # To store w(latent variables) samples
  scale_sample <- rep(NA, n_iteration)   # To store scale
  outer_logl <- rep(NA, n_iteration)     # To store logl
  w_samples[[1]] <- W
  scale_sample[1] <- scale_sampling(W, ls_y, Y, g = nugget, v = v)
  outer_logl[1] <- loglik_yw(Y, W, ls_y, g = 1e-6, v = v)
  
  for (i in 2:n_iteration) {
    w_samples[[i]] <- matrix(NA, nrow = nrow(x), ncol = node)
  }
  
  for (i in 2:n_iteration) {
    theta_y_samples[i] <- MH_1(ls = theta_y_samples[i-1],
                               Y, 
                               w = as.matrix(w_samples[[i-1]], ncol = node),
                               u = 2, l = 1, v = v)
    
    for(j in 1:node){
      theta_w_samples[i, j] <- MH_2(ls = theta_w_samples[i - 1, j],
                                    w = (matrix(w_samples[[i-1]], ncol = node))[,j], x,
                                    u = 2, l = 1, v = v)
    }
    
    w_samples[[i]] <- ESS_w(x, Y, w = as.matrix(w_samples[[i-1]], ncol = node),
                            ls_y = theta_y_samples[i],
                            ls_w = matrix(theta_w_samples[i,], ncol = node), v = v)
    
    scale_sample[i] <- scale_sampling(matrix(w_samples[[i]], ncol = node),
                                      theta_y_samples[i], 
                                      Y, 
                                      g = nugget, v = v)
    
    outer_logl[i] <- loglik_yw(Y, matrix(w_samples[[i]], ncol = node), theta_y_samples[i], g = 1e-6, v = v)
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cat("done\n")
  
  theta_y <- mean(theta_y_samples[burn_in:n_iteration])
  theta_w <- colMeans(as.matrix(theta_w_samples[burn_in:n_iteration, ]))
  
  ########## graph ############### 
  # Initialize the data frame
  df_param <- data.frame(
    iterations = c(burn_in:n_iteration),
    Theta_y = theta_y_samples[burn_in:n_iteration],
    outer_logl = outer_logl[burn_in:n_iteration]
  )
  
  n_dims <- ncol(theta_w_samples)
  # Add columns for each dimension of Theta_w
  for (i in 1:n_dims) {
    col_name <- paste0("Theta_w_", i)
    df_param[[col_name]] <- theta_w_samples[burn_in:n_iteration, i]
  }
  
  # Initialize a list to hold all the plots
  plot_list <- list()
  
  # Add the Theta_y trace plot
  plot_list[[1]] <- ggplot(df_param, aes(x = iterations, y = Theta_y)) +
    geom_line() +
    labs(title = "Trace Plot of theta_y",
         x = "Iteration",
         y = "theta_y") +
    theme_minimal()
  
  # Loop through each dimension of theta_w_samples and create a trace plot
  for (i in 1:n_dims) {
    plot_name <- paste0("Theta_w_", i)
    # Add the plot to the list
    plot_list[[i + 1]] <- ggplot(df_param, aes(x = iterations, y = !!sym(plot_name))) +
      geom_line() +
      labs(title = paste0("Trace Plot of theta_w[", i, "]"),
           x = "Iteration",
           y = paste0("theta_w[", i, "]")) +
      theme_minimal()
  }
  
  # Add the outer_logl plot
  plot_list[[length(plot_list) + 1]] <- ggplot(df_param, aes(x = iterations, y = outer_logl)) +
    geom_line() +
    labs(title = "Trace Plot of outer logl",
         x = "Iteration",
         y = "outerlogl") +
    theme_minimal()
  
  # Arrange all the plots in a grid
  p <- do.call(grid.arrange, c(plot_list, ncol = 2))
  
  ########## result #######
  result <- list(
    theta_y_samples = theta_y_samples[burn_in:n_iteration],
    theta_w_samples = theta_w_samples[burn_in:n_iteration, ],
    w_samples = w_samples[burn_in:n_iteration],
    plot_theta_y = p,
    scale = scale_sample[burn_in:n_iteration]
  )
  return(result)
}
################## Predicting #######################
Two_layer_prediction_matern <- function(list, x, Y, x_star, nugget = 1e-6, g = deepgp:::eps, v = 2.5){
  node <- ncol(x_star)
  theta_y_samples <- c(list[[1]])         # size: iteration * 1
  theta_w_samples <- as.matrix(list[[2]]) # size: iteration * p
  w_sample <- list[[3]]                   # size: a list contain number of iteration matrices, each have size of m*p
  scale_sample <- list[[5]]
  iterations <- length(theta_y_samples)
  
  cat("Predicting ... \n")
  pb <- txtProgressBar(min = 0, max = iterations, style = 3)
  
  W <- matrix(NA, nrow = nrow(x_star), ncol = node) # To store the latent variable come out from the first layer, size: m'*p
  mu <- matrix(NA, nrow = iterations, ncol = nrow(x_star))  # To store the mean, size: T * m'
  Sigma <- vector("list", iterations)  # To store the correlation matrix, is a list, contain T matrices, each have size: 
  
  for(i in 1:iterations){
    ###### layer 1 #######
    for (j in 1:node) {
      dx <- distance(x[,j])
      d_new <- distance(x_star[,j])
      d_cross <- distance(x_star[,j], x[,j])
      theta <- theta_w_samples[i,j]
      C <- deepgp:::Matern(dx, 1, theta, g, v)
      C_cross <- deepgp:::Matern(d_cross, 1, theta, nugget, v)
      C_new <- deepgp:::Matern(d_new, 1, theta, g, v)
      C_inv <- deepgp:::invdet(C)$Mi
      L <- chol(C)
      Z <- forwardsolve(t(L), t(C_cross))
      quadterm <- t(Z) %*% Z
      mean <- C_cross %*% C_inv %*% (as.matrix(w_sample[[i]]))[,j]
      sigma_w <- (C_new - quadterm)
      
      W[,j] <- matrix(mvtnorm:::rmvnorm(1, mean, sigma_w), ncol = 1) # New w, (size: m'*1). don't do this
    }
    ####### Layer 2 #######
    theta <- theta_y_samples[i]
    dw <- distance(as.matrix(w_sample[[i]]))
    dw_new <- distance(W)
    dw_cross <- distance(W, as.matrix(w_sample[[i]]))
    R <- deepgp:::Matern(dw, 1, theta, nugget, v)
    R_cross <- deepgp:::Matern(dw_cross, 1, theta, 0, v)
    R_new <- deepgp:::Matern(dw_new, 1, theta, nugget, v)
    R_inv <- deepgp:::invdet(R)$Mi
    quadterm <- R_cross %*% R_inv %*% t(R_cross)
    mu[i,] <- R_cross %*% R_inv %*% Y
    Sigma[[i]] <- scale_sample[i] * (R_new - quadterm) 
    
    setTxtProgressBar(pb, i)
  }
  
  result <- list(
    mean = colMeans(mu),
    sigma2 = (1/iterations) * Reduce(`+`, Sigma) + cov(mu)
  )
  
  return(result)
}
################## Plot the result ##################
plot_result_1 <- function(list, x_test, y_test, x_train, y_train){
  predf_1 <- data.frame(test_x = x_test, test_y = y_test,
                        pre = list$mean, sigma = sqrt((diag(list$sigma2))))
  
  my_plot <- ggplot(predf_1, aes(x = test_x)) +
    geom_point(aes(y = test_y), color = 'cyan3', size = 2) +
    geom_line(aes(y = pre), color = "red", linewidth = 1) + # Plot predicted values as a line
    geom_ribbon(aes(ymin = pre + qnorm(0.05, 0, sigma), ymax = pre + qnorm(0.95, 0, sigma)), fill = "grey", alpha = 0.5) +
    geom_line(aes(y = pre + qnorm(0.05, 0, sigma)), color = 'black', linewidth = 0.3, linetype = "dashed") +
    geom_line(aes(y = pre + qnorm(0.95, 0, sigma)), color = 'black', linewidth = 0.3, linetype = "dashed") +
    labs(title = "Actual vs Predicted Values",
         x = "X",
         y = "Y / Predicted Y") +
    theme_minimal()
  
  data_train <- dplyr::tibble(x_train = x_train, y_train = y_train)
  myplot + geom_point(data = data_train, mapping = aes(x = x_train, y = y_train), col = 'red', size = 3)
}

plot_ESS_samples1 <- function(matrix_list, x){
  matrix_list <- matrix_list
  # Create a data frame to store all the data for plotting
  data_list <- list()
  
  for (i in 1:length(matrix_list)) {
    # Extract the y-values (from each matrix in the list)
    y <- matrix_list[[i]]
    
    # Combine x and y into a data frame, along with an identifier for each matrix
    data_list[[i]] <- data.frame(x = x[,1], y = y[,1], matrix_id = i)
  }
  
  # Combine all data frames into one large data frame
  plot_data <- do.call(rbind, data_list)
  
  # Plot using ggplot2 and overlay all lines on the same plot
  ggplot(plot_data, aes(x = x, y = y, group = matrix_id)) +
    geom_line(alpha = 0.1, color = "red") +  # Set transparency with alpha
    labs(title = "ESS samples",
         x = "x",
         y = "w")
}

#################### Vecchia #######################
NNRO <- function(input_loc, k) {
  #' @description Function to find the nearest locations for the input locations after random ordering
  #' @description And here we consider the observation locations first, and prediction locations after
  #' @param input_loc matrix of observation locations (size: m * d) (After ordering)
  #' @param k the number of nearest locations to find
  #' @return NNarray a matrix where each row contains the indices of the nearest k locations in the training data
  n_obs <- nrow(input_loc)
  node <-  ncol(input_loc)
  if(k > n_obs){
    stop("Error: The number of nearest number can't be larger than the number of rows in the input matrix.")
  }
  
  NNarray <- matrix(NA, nrow = n_obs, ncol = k) # to store the row index of the nearest k locations
  NN_array <- matrix(NA, nrow = n_obs, ncol = k)
  ro_indices <- sample(n_obs)  # Random ordered row indices
  ro_input <- cbind(as.matrix(input_loc[ro_indices, ]), ro_indices, 1:n_obs)
  
  for (i in 1:nrow(ro_input)) {
    # Calculate the Euclidean distances
    distances <- sqrt(  rowSums( as.matrix( (ro_input[1:(i-1), 1:node] - ro_input[i, 1:node]) ^ 2, ncol = node) ))
    # Get the indices of the n smallest distances
    nearest_indices <- order(distances)[1:k]
    NN_array[i,] <- sort(c((ro_input[1:(i-1), node + 2])[nearest_indices]), na.last = TRUE)
  }
  return(list(k = k,
              NN_array = NN_array,
              ro_indices = ro_indices,
              random_ordered_input = ro_input))
}

MM_NN <- function(input_loc, k){
  n_obs <- nrow(input_loc)
  node <-  ncol(input_loc)
  if(k > n_obs){
    stop("Error: The number of nearest number can't be larger than the number of rows in the input matrix.")
  }
  
  max_min_ordering <- function(data) {
    n <- nrow(data)
    ordering <- integer(n)
    ordering[1] <- sample(1:n, 1)
    dist_to_selected <- rep(Inf, n)
    euclidean_dist <- function(x, y) {
      sqrt(sum((x - y)^2))
    }
    
    for (i in 2:n) {
      for (j in 1:n) {
        dist_to_selected[j] <- pmin(dist_to_selected[j], euclidean_dist(data[ordering[i-1], ], data[j, ]))
      }
      ordering[i] <- which.max(dist_to_selected)
    }
    
    return(ordering)
  }
  
  NNarray <- matrix(NA, nrow = n_obs, ncol = k) # to store the row index of the nearest k locations
  NN_array <- matrix(NA, nrow = n_obs, ncol = k)
  ro_indices <- max_min_ordering(input_loc)  # Maximum minimum ordered row indices
  ro_input <- cbind(as.matrix(input_loc[ro_indices, ]), ro_indices, 1:n_obs)
  
  for (i in 1:nrow(ro_input)) {
    # Calculate the Euclidean distances
    distances <- sqrt(  rowSums( as.matrix( (ro_input[1:(i-1), 1:node] - ro_input[i, 1:node]) ^ 2, ncol = node) ))
    # Get the indices of the n smallest distances
    nearest_indices <- order(distances)[1:k]
    NN_array[i,] <- sort(c((ro_input[1:(i-1), node + 2])[nearest_indices]), na.last = TRUE)
  }
  return(list(k = k,
              NN_array = NN_array,
              ro_indices = ro_indices,
              random_ordered_input = ro_input))
}

create_U <- function(w, NN, nugget, ls_y, g, v){
  k <- NN$k
  order <- NN$ro_indices
  w_order <- matrix(w[order, ], ncol = ncol(w))
  U <- matrix(0, ncol = nrow(w), nrow = nrow(w))  # Size of m*m
  sigma <- rep(NA, nrow(w)) # Size of m*1
  for(i in 1:nrow(w_order)){
    index <- NN$NN_array[i,]     # Get the nearest neighbor index
    index <- index[!is.na(index)]  # Get rid of the NA values
    W <- as.matrix(w_order[index, ], ncol = node) # Get the W_(c(i)) which is the conditional set
    sigma_w <- deepgp:::Matern(distance(w_order[i,]), 1, ls_y, g, v)
    sigma_W <- deepgp:::Matern(distance(W), 1, ls_y, g, v)
    sigma_wW <- deepgp:::Matern(distance(w_order[i,], W), 1, ls_y, g, v)
    B <- sigma_wW %*% deepgp:::invdet(sigma_W)$Mi
    sigma[i] <- sigma_w - B %*% t(sigma_wW)
    U[index, i] <- - (1/sqrt(sigma[i])) %*% (B)
  }
  diag(U) <- diag(U) + (1/sqrt(sigma))
  return(U)
}

fit_two_layer3.0_matern_Vecchia_RO <- function(x, Y, u = 2, l = 1, ls_y = 1, ls_w = 1, node, W, n_iteration = 7000,
                                               burn_in = 5000, nugget = 1e-6, v = 2.5, k, k_2){
  #' @description To do the training for simple two layer Deep Gaussian Processes model
  #' @param x is the input value, which have the size of M * D
  #' @param Y is the output value, which have the size of M * 1
  #' @param u is the parameter for the proposals distribution, default to 2 (according to Sauer)
  #' @param l is the parameter for the proposals distribution, default to 1 (according to Sauer)
  #' @param ls_y is the initial value for the length-scale parameter in the first layer, default to 1
  #' @param ls_w is the initial value for the length-scale parameter in the second layer, default to 1
  #' @param node is the number of latent nodes in the model
  #' @param W is the initial output value come from the second layer, is the latent variable, default to 1
  #' @param n_iteration is the number of iteration for the MCMC, default is 10000
  #' @param burn_in the iteration number for MCMC to warm up
  #' @param nugget is the nugget, set to 1e-4
  #' @param v is the smooth parameter for Matern Kernel
  #' @param k is the number of nearest neighbor for Vecchia conditioning
  #' @returns result list contain:
  #'                              1. result summary
  #'                              2. samples of theta_y
  #'                              3. samples of theta_w
  #'                              4. samples of latent variable w
  #'                              5. trace plot of theta_y
  
  library(plgp)
  library(MASS)
  library(gridExtra)
  library(ggplot2)
  library(deepgp)
  library(mvtnorm)
  
  loglik_yw <- function(Y, w, ls_y, g = 1e-6, v, k, Vecchia = FALSE){
    #' @description To compute the log-likelihood of the first layer
    n <- nrow(Y)
    if(Vecchia == FALSE){
      R <- deepgp:::Matern(distance(w), 1, ls_y, g, v)
      quadterm <- t(Y) %*% (deepgp:::invdet(R))$Mi %*% (Y)
      logl <- - 0.5 * n * log(2*pi*quadterm/n) - 0.5 * (log(det(R)))
      tau2 <- c(quadterm) / n
    } 
    else{
      NN <- NNRO(w, k)
      U <- create_U(w, NN, g, ls_y, g, v) # The upper triangular matrix of the cholesky decomposition from the precision matrix
      Y_order <- as.matrix(Y[NN$ro_indices,], nrow = n) # Let output Y follow the random order
      logdet <- sum(log(diag(U)))
      Uty <- crossprod(U, Y_order)
      ytUUty <- sum(Uty^2)
      logl <- logdet - (n * 0.5) * log(ytUUty)
      tau2 <- c(ytUUty)/n
    }
    
    return(list(logl, tau2))
  }
  
  loglik_wx <- function(w, x, ls_w, g = 1e-6, v, NN){
    #' @description To compute the log-likelihood of the second layer
    log_l <- rep(NA, node)
    n <- nrow(w)
    ro_indices <- NN$ro_indices
    w <- as.matrix(w[ro_indices,], ncol = ncol(w)) 
    U <- create_U(x, NN, g, ls_w, g, v)
    logdet <- sum(log(diag(U)))
    for (i in 1:node) {
      Utw <- crossprod(U, w[,i])
      wtUUtw <- sum(Utw^2)
      log_l[i] <-logdet - (n * 0.5) * log(wtUUtw)
    }
    return(sum(log_l))
  }
  
  
  MH_1 <- function(ls, Y, w, u = 2, l = 1, alpha = 1.5, beta = 0.65, v, k, Vecchia = FALSE){
    #' @description Function to sample the length-scale parameter using Metropolis-Hasting in layer 1, ls_y
    #' @param ls is the initial value of the length-scale parameter
    w <- as.matrix(w, ncol = node)
    ls_star <- runif(1, min = l*ls / u, max = u*ls / l) # new value
    log_alpha <- loglik_yw(Y, w, ls_star, 1e-6, v, k, Vecchia = Vecchia)[[1]] + 
      dgamma(ls_star - 1.490116e-08, alpha, beta, log = TRUE) +
      log(ls) -
      loglik_yw(Y, w, ls, 1e-6, v, k, Vecchia = Vecchia)[[1]] - 
      dgamma(ls - 1.490116e-08, alpha, beta, log = TRUE) - 
      log(ls_star)
    U <- runif(1, 0, 1)
    if(log_alpha > log(U)){ # Accepted
      return(ls_star)
    }
    else{ # Rejected
      return(ls)
    }
  }
  
  MH_2 <- function(ls, w, x, u = 2, l = 1, alpha = 1.5, beta = 0.975, v, NN){
    #' @description Function to sample the length-scale parameter using Metropolis-Hasting in layer 2, ls_w
    #' @param ls is the last value of the length-scale parameter
    w <- as.matrix(w, ncol = node)
    ls_star <- runif(1, min = l*ls / u, max = u*ls / l) # new value
    log_alpha <- loglik_wx(w, x, ls_star, 1e-6, v, NN) + 
      dgamma(ls_star - 1.490116e-08, alpha, beta, log = TRUE) + 
      log(ls) - 
      loglik_wx(w, x, ls, 1e-6, v, NN) - 
      dgamma(ls - 1.490116e-08, alpha, beta, log = TRUE) -
      log(ls_star)
    U <- runif(1, 0, 1)
    if(log_alpha > log(U)){ # Accepted
      return(ls_star)
    }
    else{ # Rejected
      return(ls)
    }
  }
  
  ESS_w <- function(x, Y, w, ls_y, ls_w, p = node, g = nugget, v, k, Vecchia = FALSE){
    #' @description To do the Elliptical Slice Sampling update for latent variable w
    #' @param x is the input of the model
    #' @param Y is the Output value of the model
    #' @param w is the initial value for the latent variable
    #' @param ls_y is the length-scale parameter sampled from MH_1
    #' @param ls_w is the length-scale parameter sampled from MH_2
    #' @param p is the nodes in the latent layer
    #' @param g is the nugget
    #' @return w is the latent variable matrix(size m*p) after the ESS
    
    m <- nrow(w)
    if (p != ncol(w)){
      stop("Error! The initial value for the latent value doesn't match the nodes number!")
    }
    for (i in 1:p) {
      theta <- runif(1, 0, 2*pi) # angle
      theta_min <- theta - 2*pi  # lower basket
      theta_max <- theta         # upper basket
      nu <- mvtnorm::rmvnorm(1, mean = matrix(0, nrow = nrow(w)), sigma = deepgp:::Matern(distance(x), 1, ls_w[i], 0, v)) # random draw from the prior, size = m*1
      w_prev <- w[, i]
      ll_prev <- loglik_yw(Y, w, ls_y, g, v, k, Vecchia = Vecchia)[[1]]
      accept <- FALSE
      count <- 0
      while (accept == FALSE){
        count <- count + 1
        w[,i] <- w_prev * cos(theta) + nu * sin(theta) # Proposal
        dw <- deepgp:::sq_dist(w)
        log_alpha <- loglik_yw(Y, w, ls_y, g, v, k = k, Vecchia = Vecchia)[[1]] - ll_prev # log-alpha
        U <- runif(1, 0, 1)
        # Check if the proposed sample is on the slice
        if (log_alpha > log(U)) # Accepted
        { 
          accept <- TRUE
        } 
        else { # Rejected
          # Shrink the bracket
          if (theta < 0) {
            theta_min <- theta
          } else {
            theta_max <- theta
          }
          # Draw a new angle from the updated bracket
          theta <- runif(1, theta_min, theta_max)
        }
      }
    }
    return(w)
  }
  
  cat("Training ... \n")
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_iteration, style = 3)
  
  theta_y_samples <- c(ls_y, rep(NA, n_iteration - 1))                              # To store theta_y samples (layer 1)
  theta_w_samples <- rbind(ls_w, matrix(NA, ncol = node, nrow = n_iteration - 1))   # To store theta_w samples (layer 2)
  w_samples <- vector("list", n_iteration)                                          # To store w(latent variables) samples
  scale_sample <- matrix(NA, nrow = n_iteration, ncol = ncol(x))                    # To store scale
  outer_logl <- rep(NA, n_iteration)                                                # To store outer layer logl
  w_samples[[1]] <- W
  scale_sample[1,] <- loglik_yw(Y, W, ls_y, g = 1e-6, v = v, k = k)[[2]]
  outer_logl[1] <- loglik_yw(Y, W, ls_y, g = 1e-6, v = v, k = k)[[1]]
  
  
  for (i in 2:n_iteration) {
    w_samples[[i]] <- matrix(NA, nrow = nrow(x), ncol = node)
  }
  NN <- NNRO(x, k_2)
  for (i in 2:n_iteration) {
    if(i %% 2 == 1){
      vec <- FALSE
    } else{
      vec <- TRUE
    }
    theta_y_samples[i] <- MH_1(ls = theta_y_samples[i-1],
                               Y, 
                               w = as.matrix(w_samples[[i-1]], ncol = node),
                               u = 2, l = 1, v = v, k = k, Vecchia = vec)
    for(j in 1:node){
      theta_w_samples[i, j] <- MH_2(ls = theta_w_samples[i - 1, j],
                                    w = as.matrix(w_samples[[i-1]], ncol = node), x,
                                    u = 2, l = 1, v = v, NN = NN)
    }
    
    w_samples[[i]] <- ESS_w(x, Y, w = as.matrix(w_samples[[i-1]], ncol = node),
                            ls_y = theta_y_samples[i],
                            ls_w = theta_w_samples[i,], v = v, k = k, Vecchia = vec)
    
    scale_sample[i,] <- loglik_yw(Y, matrix(w_samples[[i]], ncol = node), theta_y_samples[i], g = 1e-6, v = v, k, Vecchia = vec)[[2]]
    
    outer_logl[i] <- loglik_yw(Y, matrix(w_samples[[i]], ncol = node), theta_y_samples[i], g = 1e-6, v = v, k, Vecchia = vec)[[1]]
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("done\n")
  
  df_param <- data.frame(
    Theta_y = theta_y_samples[burn_in:n_iteration],
    Theta_w_1 = theta_w_samples[burn_in:n_iteration,1],
    outer_logl = outer_logl[burn_in:n_iteration]
  )
  
  theta_y <- mean(theta_y_samples[burn_in:n_iteration])
  theta_w <- colMeans(as.matrix(theta_w_samples[burn_in:n_iteration, ]))
  
  result_summary <- data.frame(
    "Layer No." = c("Layer 1", "Layer 2"),
    "Kernel" = c("Matérn", "Matérn"),
    "Length-scale" = c(theta_w, theta_y), 
    "Variance" = c("1 (fixed)", mean(scale_sample)),
    "Nugget" = c("1e-6 (fixed)", "1e-6 (fixed)"),
    "Input Dims" = c(ncol(W), ncol(x))
  )
  
  ########## graph ############### 
  df_param <- data.frame(
    iterations = c(burn_in:n_iteration),
    Theta_y = theta_y_samples[burn_in:n_iteration],
    outer_logl = outer_logl[burn_in:n_iteration]
  )
  
  n_dims <- ncol(theta_w_samples)
  # Add columns for each dimension of Theta_w
  for (i in 1:n_dims) {
    col_name <- paste0("Theta_w_", i)
    df_param[[col_name]] <- theta_w_samples[burn_in:n_iteration, i]
  }
  
  # Initialize a list to hold all the plots
  plot_list <- list()
  
  # Add the Theta_y trace plot
  plot_list[[1]] <- ggplot(df_param, aes(x = iterations, y = Theta_y)) +
    geom_line() +
    labs(title = "Trace Plot of theta_y",
         x = "Iteration",
         y = "theta_y") +
    theme_minimal()
  
  # Loop through each dimension of theta_w_samples and create a trace plot
  for (i in 1:n_dims) {
    plot_name <- paste0("Theta_w_", i)
    # Add the plot to the list
    plot_list[[i + 1]] <- ggplot(df_param, aes(x = iterations, y = !!sym(plot_name))) +
      geom_line() +
      labs(title = paste0("Trace Plot of theta_w[", i, "]"),
           x = "Iteration",
           y = paste0("theta_w[", i, "]")) +
      theme_minimal()
  }
  
  # Add the outer_logl plot
  plot_list[[length(plot_list) + 1]] <- ggplot(df_param, aes(x = iterations, y = outer_logl)) +
    geom_line() +
    labs(title = "Trace Plot of outer logl",
         x = "Iteration",
         y = "outerlogl") +
    theme_minimal()
  
  # Arrange all the plots in a grid
  p <- do.call(grid.arrange, c(plot_list, ncol = 2))
  ########## result #######
  result <- list(
    theta_y_samples = theta_y_samples[burn_in:n_iteration],
    theta_w_samples = theta_w_samples[burn_in:n_iteration, ],
    w_samples = w_samples[burn_in:n_iteration],
    plot_theta_y = p,
    scale = scale_sample[burn_in:n_iteration]
  )
  
  return(result)
}

fit_two_layer3.0_matern_Vecchia_MM <- function(x, Y, u = 2, l = 1, ls_y = 1, ls_w = 1, node, W, n_iteration = 7000,
                                               burn_in = 5000, nugget = 1e-6, v = 2.5, k, k_2){
  #' @description To do the training for simple two layer Deep Gaussian Processes model
  #' @param x is the input value, which have the size of M * D
  #' @param Y is the output value, which have the size of M * 1
  #' @param u is the parameter for the proposals distribution, default to 2 (according to Sauer)
  #' @param l is the parameter for the proposals distribution, default to 1 (according to Sauer)
  #' @param ls_y is the initial value for the length-scale parameter in the first layer, default to 1
  #' @param ls_w is the initial value for the length-scale parameter in the second layer, default to 1
  #' @param node is the number of latent nodes in the model
  #' @param W is the initial output value come from the second layer, is the latent variable, default to 1
  #' @param n_iteration is the number of iteration for the MCMC, default is 10000
  #' @param burn_in the iteration number for MCMC to warm up
  #' @param nugget is the nugget, set to 1e-4
  #' @param v is the smooth parameter for Matern Kernel
  #' @param k is the number of nearest neighbor for Vecchia conditioning
  #' @returns result list contain:
  #'                              1. result summary
  #'                              2. samples of theta_y
  #'                              3. samples of theta_w
  #'                              4. samples of latent variable w
  #'                              5. trace plot of theta_y
  
  library(plgp)
  library(MASS)
  library(gridExtra)
  library(ggplot2)
  library(deepgp)
  library(mvtnorm)
  
  loglik_yw <- function(Y, w, ls_y, g = 1e-6, v, k, Vecchia = FALSE){
    #' @description To compute the log-likelihood of the first layer
    n <- nrow(Y)
    if(Vecchia == FALSE){
      R <- deepgp:::Matern(distance(w), 1, ls_y, g, v)
      quadterm <- t(Y) %*% (deepgp:::invdet(R))$Mi %*% (Y)
      logl <- - 0.5 * n * log(2*pi*quadterm/n) - 0.5 * (log(det(R)))
      tau2 <- c(quadterm) / n
    } 
    else{
      NN <- MM_NN(w, k)
      U <- create_U(w, NN, g, ls_y, g, v) # The upper triangular matrix of the cholesky decomposition from the precision matrix
      Y_order <- as.matrix(Y[NN$ro_indices,], nrow = n) # Let output Y follow the random order
      logdet <- sum(log(diag(U)))
      Uty <- crossprod(U, Y_order)
      ytUUty <- sum(Uty^2)
      logl <- logdet - (n * 0.5) * log(ytUUty)
      tau2 <- c(ytUUty)/n
    }
    
    return(list(logl, tau2))
  }
  
  loglik_wx <- function(w, x, ls_w, g = 1e-6, v, NN){
    #' @description To compute the log-likelihood of the second layer
    log_l <- rep(NA, node)
    n <- nrow(w)
    ro_indices <- NN$ro_indices
    w <- as.matrix(w[ro_indices,], ncol = ncol(w)) 
    U <- create_U(x, NN, g, ls_w, g, v)
    logdet <- sum(log(diag(U)))
    for (i in 1:node) {
      Utw <- crossprod(U, w[,i])
      wtUUtw <- sum(Utw^2)
      log_l[i] <-logdet - (n * 0.5) * log(wtUUtw)
    }
    return(sum(log_l))
  }
  
  
  MH_1 <- function(ls, Y, w, u = 2, l = 1, alpha = 1.5, beta = 0.65, v, k, Vecchia = FALSE){
    #' @description Function to sample the length-scale parameter using Metropolis-Hasting in layer 1, ls_y
    #' @param ls is the initial value of the length-scale parameter
    w <- as.matrix(w, ncol = node)
    ls_star <- runif(1, min = l*ls / u, max = u*ls / l) # new value
    log_alpha <- loglik_yw(Y, w, ls_star, 1e-6, v, k, Vecchia = Vecchia)[[1]] + 
      dgamma(ls_star - 1.490116e-08, alpha, beta, log = TRUE) +
      log(ls) -
      loglik_yw(Y, w, ls, 1e-6, v, k, Vecchia = Vecchia)[[1]] - 
      dgamma(ls - 1.490116e-08, alpha, beta, log = TRUE) - 
      log(ls_star)
    U <- runif(1, 0, 1)
    if(log_alpha > log(U)){ # Accepted
      return(ls_star)
    }
    else{ # Rejected
      return(ls)
    }
  }
  
  MH_2 <- function(ls, w, x, u = 2, l = 1, alpha = 1.5, beta = 0.975, v, NN){
    #' @description Function to sample the length-scale parameter using Metropolis-Hasting in layer 2, ls_w
    #' @param ls is the last value of the length-scale parameter
    w <- as.matrix(w, ncol = node)
    ls_star <- runif(1, min = l*ls / u, max = u*ls / l) # new value
    log_alpha <- loglik_wx(w, x, ls_star, 1e-6, v, NN) + 
      dgamma(ls_star - 1.490116e-08, alpha, beta, log = TRUE) + 
      log(ls) - 
      loglik_wx(w, x, ls, 1e-6, v, NN) - 
      dgamma(ls - 1.490116e-08, alpha, beta, log = TRUE) -
      log(ls_star)
    U <- runif(1, 0, 1)
    if(log_alpha > log(U)){ # Accepted
      return(ls_star)
    }
    else{ # Rejected
      return(ls)
    }
  }
  
  ESS_w <- function(x, Y, w, ls_y, ls_w, p = node, g = nugget, v, k, Vecchia = FALSE){
    #' @description To do the Elliptical Slice Sampling update for latent variable w
    #' @param x is the input of the model
    #' @param Y is the Output value of the model
    #' @param w is the initial value for the latent variable
    #' @param ls_y is the length-scale parameter sampled from MH_1
    #' @param ls_w is the length-scale parameter sampled from MH_2
    #' @param p is the nodes in the latent layer
    #' @param g is the nugget
    #' @return w is the latent variable matrix(size m*p) after the ESS
    
    m <- nrow(w)
    if (p != ncol(w)){
      stop("Error! The initial value for the latent value doesn't match the nodes number!")
    }
    for (i in 1:p) {
      theta <- runif(1, 0, 2*pi) # angle
      theta_min <- theta - 2*pi  # lower basket
      theta_max <- theta         # upper basket
      nu <- mvtnorm::rmvnorm(1, mean = matrix(0, nrow = nrow(w)), sigma = deepgp:::Matern(distance(x), 1, ls_w[i], 0, v)) # random draw from the prior, size = m*1
      w_prev <- w[, i]
      ll_prev <- loglik_yw(Y, w, ls_y, g, v, k, Vecchia = Vecchia)[[1]]
      accept <- FALSE
      count <- 0
      while (accept == FALSE){
        count <- count + 1
        w[,i] <- w_prev * cos(theta) + nu * sin(theta) # Proposal
        dw <- deepgp:::sq_dist(w)
        log_alpha <- loglik_yw(Y, w, ls_y, g, v, k = k, Vecchia = Vecchia)[[1]] - ll_prev # log-alpha
        U <- runif(1, 0, 1)
        # Check if the proposed sample is on the slice
        if (log_alpha > log(U)) # Accepted
        { 
          accept <- TRUE
        } 
        else { # Rejected
          # Shrink the bracket
          if (theta < 0) {
            theta_min <- theta
          } else {
            theta_max <- theta
          }
          # Draw a new angle from the updated bracket
          theta <- runif(1, theta_min, theta_max)
        }
      }
    }
    return(w)
  }
  
  cat("Training ... \n")
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_iteration, style = 3)
  
  theta_y_samples <- c(ls_y, rep(NA, n_iteration - 1))                              # To store theta_y samples (layer 1)
  theta_w_samples <- rbind(ls_w, matrix(NA, ncol = node, nrow = n_iteration - 1))   # To store theta_w samples (layer 2)
  w_samples <- vector("list", n_iteration)                                          # To store w(latent variables) samples
  scale_sample <- matrix(NA, nrow = n_iteration, ncol = ncol(x))                    # To store scale
  outer_logl <- rep(NA, n_iteration)                                                # To store outer layer logl
  w_samples[[1]] <- W
  scale_sample[1,] <- loglik_yw(Y, W, ls_y, g = 1e-6, v = v, k = k)[[2]]
  outer_logl[1] <- loglik_yw(Y, W, ls_y, g = 1e-6, v = v, k = k)[[1]]
  
  
  for (i in 2:n_iteration) {
    w_samples[[i]] <- matrix(NA, nrow = nrow(x), ncol = node)
  }
  NN <- MM_NN(x, k_2)
  for (i in 2:n_iteration) {
    if(i %% 2 == 1){
      vec <- FALSE
    } else{
      vec <- TRUE
    }
    theta_y_samples[i] <- MH_1(ls = theta_y_samples[i-1],
                               Y, 
                               w = as.matrix(w_samples[[i-1]], ncol = node),
                               u = 2, l = 1, v = v, k = k, Vecchia = vec)
    for(j in 1:node){
      theta_w_samples[i, j] <- MH_2(ls = theta_w_samples[i - 1, j],
                                    w = as.matrix(w_samples[[i-1]], ncol = node), x,
                                    u = 2, l = 1, v = v, NN = NN)
    }
    
    w_samples[[i]] <- ESS_w(x, Y, w = as.matrix(w_samples[[i-1]], ncol = node),
                            ls_y = theta_y_samples[i],
                            ls_w = theta_w_samples[i,], v = v, k = k, Vecchia = vec)
    
    scale_sample[i,] <- loglik_yw(Y, matrix(w_samples[[i]], ncol = node), theta_y_samples[i], g = 1e-6, v = v, k, Vecchia = vec)[[2]]
    
    outer_logl[i] <- loglik_yw(Y, matrix(w_samples[[i]], ncol = node), theta_y_samples[i], g = 1e-6, v = v, k, Vecchia = vec)[[1]]
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("done\n")
  
  df_param <- data.frame(
    Theta_y = theta_y_samples[burn_in:n_iteration],
    Theta_w_1 = theta_w_samples[burn_in:n_iteration,1],
    outer_logl = outer_logl[burn_in:n_iteration]
  )
  
  theta_y <- mean(theta_y_samples[burn_in:n_iteration])
  theta_w <- colMeans(as.matrix(theta_w_samples[burn_in:n_iteration, ]))
  
  result_summary <- data.frame(
    "Layer No." = c("Layer 1", "Layer 2"),
    "Kernel" = c("Matérn", "Matérn"),
    "Length-scale" = c(theta_w, theta_y), 
    "Variance" = c("1 (fixed)", mean(scale_sample)),
    "Nugget" = c("1e-6 (fixed)", "1e-6 (fixed)"),
    "Input Dims" = c(ncol(W), ncol(x))
  )
  
  ########## graph ############### 
  df_param <- data.frame(
    iterations = c(burn_in:n_iteration),
    Theta_y = theta_y_samples[burn_in:n_iteration],
    outer_logl = outer_logl[burn_in:n_iteration]
  )
  
  n_dims <- ncol(theta_w_samples)
  # Add columns for each dimension of Theta_w
  for (i in 1:n_dims) {
    col_name <- paste0("Theta_w_", i)
    df_param[[col_name]] <- theta_w_samples[burn_in:n_iteration, i]
  }
  
  # Initialize a list to hold all the plots
  plot_list <- list()
  
  # Add the Theta_y trace plot
  plot_list[[1]] <- ggplot(df_param, aes(x = iterations, y = Theta_y)) +
    geom_line() +
    labs(title = "Trace Plot of theta_y",
         x = "Iteration",
         y = "theta_y") +
    theme_minimal()
  
  # Loop through each dimension of theta_w_samples and create a trace plot
  for (i in 1:n_dims) {
    plot_name <- paste0("Theta_w_", i)
    # Add the plot to the list
    plot_list[[i + 1]] <- ggplot(df_param, aes(x = iterations, y = !!sym(plot_name))) +
      geom_line() +
      labs(title = paste0("Trace Plot of theta_w[", i, "]"),
           x = "Iteration",
           y = paste0("theta_w[", i, "]")) +
      theme_minimal()
  }
  
  # Add the outer_logl plot
  plot_list[[length(plot_list) + 1]] <- ggplot(df_param, aes(x = iterations, y = outer_logl)) +
    geom_line() +
    labs(title = "Trace Plot of outer logl",
         x = "Iteration",
         y = "outerlogl") +
    theme_minimal()
  
  # Arrange all the plots in a grid
  p <- do.call(grid.arrange, c(plot_list, ncol = 2))
  ########## result #######
  result <- list(
    theta_y_samples = theta_y_samples[burn_in:n_iteration],
    theta_w_samples = theta_w_samples[burn_in:n_iteration, ],
    w_samples = w_samples[burn_in:n_iteration],
    plot_theta_y = p,
    scale = scale_sample[burn_in:n_iteration]
  )
  
  return(result)
}

Vecchia_plot_result_1 <- function(list, x_test, y_test, x_train, y_train){
  predf_1 <- data.frame(test_x = x_test, test_y = y_test,
                        pre = list$mean, sigma = sqrt((diag(list$sigma2))))
  
  myplot <- ggplot(predf_1, aes(x = test_x)) +
    geom_point(aes(y = test_y), color = 'cyan3', size = 2) +
    geom_line(aes(y = pre), color = "red", linewidth = 1) + # Plot predicted values as a line
    geom_ribbon(aes(ymin = pre + qnorm(0.05, 0, sigma), ymax = pre + qnorm(0.95, 0, sigma)), fill = "grey", alpha = 0.5) +
    geom_line(aes(y = pre + qnorm(0.05, 0, sigma)), color = 'black', linewidth = 0.3, linetype = "dashed") +
    geom_line(aes(y = pre + qnorm(0.95, 0, sigma)), color = 'black', linewidth = 0.3, linetype = "dashed") +
    labs(title = "Actual vs Predicted Values(Vecchia)",
         x = "X",
         y = "Y / Predicted Y") +
    theme_minimal()
  
  myplot <- Vecchia_plot_result_1(pre_Vec, x_test, y_test)
  data_train <- dplyr::tibble(x_train = x_train, y_train = y_train)
  myplot + geom_point(data = data_train, mapping = aes(x = x_train, y = y_train), col = 'red', size = 3)
}


####################################################
####################################################
################# One-D example ####################
higdon <- function(x) {
  i <- which(x <= 0.6)
  x[i] <- 2 * sin(pi * 0.8 * x[i] * 4) + 0.4 * cos(pi * 0.8 * x[i] * 16)
  x[-i] <- 2 * x[-i] - 1
  return(x)
}

# Training data
n <- 24
x_train <- matrix(seq(0, 1, length = n), ncol = 1)
y_train <- matrix(higdon(x_train), ncol = 1)

# Testing data
np <- 100
x_test <- matrix(seq(0, 1, length = np), ncol = 1)
y_test <- matrix(higdon(x_test), ncol = 1)

plot(x_test, y_test, type = "l", col = 4, xlab = "X", ylab = "Y", main = "Higdon function")
points(x_train, y_train)

dgp_1 <- fit_two_layer3.0_matern(x = x_train, Y = y_train, u = 2, l = 1, ls_y = 0.1, ls_w = 0.1, node = 1,
                                 W = x_train, n_iteration = 10000, burn_in = 7000, nugget = 1e-6, v = 2.5)
pre_1 <- Two_layer_prediction_matern(dgp_1, x_train, y_train, x_test)

myplot <- plot_result_1(pre_1, x_test, y_test, x_train, y_train)

plot_ESS_samples1(dgp_1$w_samples, x_train)


####################################################
####################################################
############### results comparison #################
higdon <- function(x) {
  i <- which(x <= 0.6)
  x[i] <- 2 * sin(pi * 0.8 * x[i] * 4) + 0.4 * cos(pi * 0.8 * x[i] * 16)
  x[-i] <- 2 * x[-i] - 1
  return(x)
}

# Training data
n <- 30
x_train <- matrix(seq(0, 1, length = n), ncol = 1)
y_train <- matrix(higdon(x_train), ncol = 1)

# Testing data
np <- 300
x_test <- matrix(seq(0, 1, length = np), ncol = 1)
y_test <- matrix(higdon(x_test), ncol = 1)

plot(x_test, y_test, type = "l", col = 4, xlab = "X", ylab = "Y", main = "Higdon function")
points(x_train, y_train)

####### GP emulator ########
gp1 <- GPemulator_MH_training(x_train, y_train, n_iteration = 10000, burn_in = 5000)
pre1 <- GPprediction(gp1, x_train, y_train, x_test, 1e-8)

rmse1 <- rmse(c(y_test), c(pre1$mean))

####### DGP emulator #######
dgp1 <- fit_two_layer3.0_matern(x = x_train, Y = y_train, u = 2, l = 1, ls_y = 0.1, ls_w = 0.1, node = 1,
                                W = x_train, n_iteration = 10000, burn_in = 7000, nugget = 1e-6, v = 2.5)
pre2 <- Two_layer_prediction_matern(dgp1, x_train, y_train, x_test)

rmse2 <- rmse(c(y_test), c(pre2$mean))

###### DGP using Vecchia ######
dgp_Vecchia_RO <- fit_two_layer3.0_matern_Vecchia_RO(x = x_train, Y = y_train, ls_y = 0.1, ls_w = 0.1, node = 1,
                                                     W = x_train, v = 2.5, k = 10, k_2 = 10)

pre_Vec_RO <- Two_layer_prediction_matern(dgp_Vecchia_RO, x_train, y_train, x_test)

rmse3 <- rmse(c(y_test), c(pre_Vec_RO$mean))

dgp_Vecchia_MM <- fit_two_layer3.0_matern_Vecchia_MM(x = x_train, Y = y_train, ls_y = 0.1, ls_w = 0.1, node = 1,
                                                     W = x_train, v = 2.5, k = 10, k_2 = 10)

pre_Vec_MM <- Two_layer_prediction_matern(dgp_Vecchia_MM, x_train, y_train, x_test)

rmse4 <- rmse(c(y_test), c(pre_Vec_MM$mean))

rmse_df <- data.frame(
  Model = c("GP", "Tow-layer DGP", 
            "Tow_layer DGP using Vecchia(Random ordering)", 
            "Tow_layer DGP using Vecchia(Maxmin ordering)"),  # Labels for each model (optional)
  RMSE = c(rmse1, rmse2, rmse3, rmse4)  # RMSE values
)
rmse_df


#################################################################
####### dgpsi version compare of using Vec and not use it #######
library(tidyr)
library(dplyr)
library(ggplot2)
n <- 1000 # number of experiment
rmse_dgpsi <- rep(NA, n)
rmse_dgpsi_vec <- rep(NA, n)
for (i in 1:n) {
  dgpsi2_vec <- dgpsi:::dgp(x_train, y_train, name = 'matern2.5', vecchia = TRUE, M = 10)
  p_vec <- predict(dgpsi2_vec, x = x_test)
  mu_vec <- p_vec$results$mean
  rmse_dgpsi_vec[i] <- deepgp::rmse(y_test, mu_vec)
  dgpsi2 <- dgpsi::dgp(x_train, y_train, name = 'matern2.5')
  p <- predict(dgpsi2, x = x_test)
  mu <- p$results$mean
  rmse_dgpsi[i] <- deepgp::rmse(y_test, mu)
}
dgpsi_rmse_df <- data.frame(
  Model = c("DGP", 
            "DGP using Vecchia"),  # Labels for each model (optional)
  RMSE = c(mean(rmse_dgpsi), mean(rmse_dgpsi_vec))  # RMSE values
)
dgpsi_rmse_df # data frame of DGP using and not using Vecchia

df <- data.frame(
  number_of_experiment = c(1:n),
  rmse_dgpsi = rmse_dgpsi,
  rmse_dgpsi_vec = rmse_dgpsi_vec
)

# Reshape data from wide to long format
df_long <- df %>%
  pivot_longer(cols = c(rmse_dgpsi, rmse_dgpsi_vec), 
               names_to = "Method", 
               values_to = "RMSE")

# Calculate the mean RMSE for each method
mean_rmse <- df_long %>%
  group_by(Method) %>%
  summarize(mean_RMSE = mean(RMSE, na.rm = TRUE))

# Plot using ggplot and add mean lines
p <- ggplot(df_long, aes(x = number_of_experiment, y = RMSE, color = Method)) +
  geom_line() +
  labs(title = "Comparison of RMSE for Two Methods",
       x = "Number of Experiment", 
       y = "RMSE") +
  theme_minimal() +
  # Add mean lines for each method
  geom_hline(data = mean_rmse, aes(yintercept = mean_RMSE, color = Method), 
             linetype = "dashed", linewidth = 1)

# Display the plot
p


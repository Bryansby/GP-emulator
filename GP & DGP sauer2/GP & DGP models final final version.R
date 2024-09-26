######### Build the synthetics functions #########
# Model 1
f1 <- function(x) {
  (sin(7.5*x)+1)/2
}
# Model 2
f2 <- function(x) {
  2/3*sin(2*(2*x - 1))+4/3*exp(-30*(2*(2*x-1))^2)-1/3  
}
# Model 3
f3 <- function(x) {
  x[1]*x[2]^2
}
# Linked Model 
f123 <- function(x) {
  f3(c(f1(x),f2(f1(x))))
}

set.seed(123)

############# Generate Training points using Latin Hyper cube Sampling ##############
n_train <- 30  # Number of training points

# Generate LHS samples in 1D (you can extend this to higher dimensions if needed)
X_gp <- maximinLHS(n_train, 1)  # Latin Hypercube Sampling in 1D, scaled between [0, 1]

# Scale LHS points to [0, 1] range
X_gp <- as.matrix(X_gp)  # Ensure it is in matrix form for further calculations

# Compute the corresponding Y values for the training set using f123
Y_gp <- matrix(sapply(X_gp, f123), ncol = 1)

############ Prediction points ############
test_x <- as.matrix(seq(0, 1, length = 300), ncol = 1)
test_y <- as.matrix(sapply(test_x, f123), ncol = 1)

# Create a data frame for ggplot
df <- data.frame(x = test_x, y = test_y)

# Plot using ggplot2
library(ggplot2)
ggplot(df, aes(x = x, y = y)) +
  geom_line(color = 'red', linewidth = 1) +
  ggtitle('Plot of f123(x)') +
  xlab('x') +
  ylab('f123(x)') +
  theme_minimal()

################## Plot the result ##################
plot_result_1 <- function(list, x_train, y_train, x_test, y_test){
  predf_1 <- data.frame(test_x = x_test, test_y = y_test,
                        pre_mean = list$mean, 
                        s2 = sqrt((list$s2)))
  
  myplot <- ggplot(predf_1, aes(x = test_x)) +
    geom_ribbon(aes(ymin = pre_mean + qnorm(0.05, 0, s2), ymax = pre_mean + qnorm(0.95, 0, s2)), fill = "grey", alpha = 0.5) +
    #geom_line(aes(y = pre_mean), color = "red", linewidth = 1) + # predict mean
    geom_line(aes(y = pre_mean), color = "brown2", linewidth = 0.8) + # predict mean
    geom_line(aes(y = test_y), color = "orange", linewidth = 0.5) +
    geom_line(aes(y = pre_mean + qnorm(0.05, 0, s2)), color = 'black', linewidth = 0.3, linetype = "dashed") +
    geom_line(aes(y = pre_mean + qnorm(0.95, 0, s2)), color = 'black', linewidth = 0.3, linetype = "dashed") +
    labs(title = "Actual vs Predicted Values",
         x = "X",
         y = "Y / Predicted Y") +
    theme_minimal()
  
  data_train <- dplyr::tibble(x_train = x_train, y_train = y_train)
  myplot + geom_point(data = data_train, mapping = aes(x = x_train, y = y_train), col = 'deepskyblue2', size = 2)
  
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
  
  # Plot using ggplot2 and map the line color to matrix_id with a color gradient
  ggplot(plot_data, aes(x = x, y = y, group = matrix_id, color = matrix_id)) +
    geom_line(alpha = 0.6) +  # Set transparency with alpha
    scale_color_gradient(low = "red", high = "yellow") +  # Gradient from red to yellow
    labs(title = "ESS samples",
         x = "x",
         y = "w") +
    theme_minimal()  # Optional: improve plot aesthetics
}

########## GP ###############
###### Training ############
GPemulator_MH_training <- function(x, Y, n_iteration = 5000, burn_in = 3000, ls_y = matrix(rep(0.1, ncol(x)), nrow = 1), nugget = 1e-8, g = deepgp:::eps, v = 2.5){
  dim <- ncol(x)
  n <- nrow(Y)
  ll_store <- rep(NA, length = n_iteration)
  
  loglik_xy <- function(x, y, ls_y, g = nugget){
    #' @description To compute the log-likelihood of the GP emulator
    n <- nrow(y)
    R <- deepgp:::MaternSep(x, x, 1, ls_y, g, v)
    id <- deepgp:::invdet(R)
    quadterm <- t(y) %*% id$Mi %*% (y)
    log_l <- (- n * 0.5) * log(quadterm) - 0.5 * (id$ldet)
    return(log_l)
  }
  
  MH <- function(Ls, x, y, index, u = 2, l = 1, alpha = 1.5, beta = 2.6, g = nugget){
    #' @description Function to sample the length-scale parameter using Metropolis-Hasting, only for single dimension
    eps <- deepgp:::eps
    ls_star <- runif(1, min = l*(Ls[index]) / u, max = u*(Ls[index]) / l) # new value
    ls_updated <- Ls
    ls_updated[index] <- ls_star
    log_alpha <- loglik_xy(x, y, ls_updated, g) + 
      dgamma(ls_star - eps, alpha, beta, log = TRUE) +
      log(Ls[index]) -
      loglik_xy(x, y, Ls, g) - 
      dgamma(Ls[index] - eps, alpha, beta, log = TRUE) - 
      log(ls_star)
    U <- runif(1, 0, 1)
    if(log_alpha > log(U)){ # Accepted
      return(ls_star)
    }
    else{ # Rejected
      return(Ls[index])
    }
  }
  
  cat("Training ... \n")
  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = n_iteration, style = 3)
  
  theta_y_samples <- matrix(NA, nrow = n_iteration, ncol = dim)
  theta_y_samples[1,] <- ls_y
  sc <- rep(NA, n_iteration)
  R <- deepgp:::Exp2Sep(x, x, tau2 = 1, theta = theta_y_samples[1,], g = nugget)
  id <- deepgp:::invdet(R)
  quadterm <- t(Y) %*% id$Mi %*% Y
  sc[1] <- c(quadterm) / n
  
  for (i in 2:n_iteration) {
    for(j in 1:dim){
      theta_y_samples[i, j] <- MH(Ls = theta_y_samples[i-1, ],
                                  x = x, y = Y, index = j,
                                  u = 2, l = 1, g = nugget)
    }
    R <- deepgp:::MaternSep(x, x, tau2 = 1, theta = theta_y_samples[i,], g = nugget, v)
    id <- deepgp:::invdet(R)
    quadterm <- t(Y) %*% id$Mi %*% Y
    sc[i] <- c(quadterm) / n
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("done\n")
  
  if(dim == 1){
    ls <- mean(theta_y_samples[burn_in:n_iteration, ])
  } else{
    ls <- colMeans(theta_y_samples[burn_in:n_iteration, ])
  }
  
  #### graph #####
  # Initialize the data frame
  df_param <- data.frame(
    iterations = c(burn_in:n_iteration)
  )
  
  # Add columns for each dimension of Theta_w
  for (i in 1:dim) {
    col_name <- paste0("Theta_y_", i)
    df_param[[col_name]] <- theta_y_samples[burn_in:n_iteration, i]
  }
  
  # Initialize a list to hold all the plots
  plot_list <- list()
  
  # Loop through each dimension of theta_w_samples and create a trace plot
  for (i in 1:dim) {
    plot_name <- paste0("Theta_y_", i)
    # Add the plot to the list
    plot_list[[i]] <- ggplot(df_param, aes(x = iterations, y = !!sym(plot_name))) +
      geom_line() +
      labs(title = paste0("Trace Plot of theta[", i, "]"),
           x = "Iteration",
           y = paste0("theta_y")) +
      theme_minimal()
  }
  
  # Arrange all the plots in a grid
  p <- do.call(grid.arrange, c(plot_list, ncol = 2))
  
  ###### result #####
  result_summary <- data.frame(
    length_scale = ls,
    scale = mean(sc[burn_in:n_iteration]),
    Nugget = 1e-8
  )
  result <- list(
    result_summary = result_summary,
    length_scale = ls,
    plot = p,
    scale = mean(sc[burn_in:n_iteration]),
    theta_y_samples = matrix(theta_y_samples[burn_in:n_iteration,], ncol = dim),
    scale_samples = sc[burn_in:n_iteration],
    dim = dim
  )
  return(result)
}
###### Prediction ##########
GPprediction <- function(list, x, Y, x_star, nugget, g, v = 2.5){
  theta_y_samples <- list[[5]]  # size : iteration * dim
  scale_samples <- list[[6]]    # size : iteration * 1
  dim <- list[[7]]
  iterations <- length(theta_y_samples)
  
  cat("Predicting ... \n")
  pb <- txtProgressBar(min = 0, max = iterations, style = 3)
  
  mu <- matrix(NA, nrow = iterations, ncol = nrow(x_star))  # To store the mean, size: T * m'
  #Sigma <- vector("list", iterations)  # To store the correlation matrix, is a list, contain T matrices, each have size: 
  #sigma_sum <- matrix(0, nrow = nrow(x_star), ncol = nrow(x_star))
  s2_sum <- rep(0, times = nrow(x_star))
  
  for(i in 1:iterations){
    theta <- theta_y_samples[i,]
    tau2 <- scale_samples[i]
    C <- deepgp:::MaternSep(x, x, 1, theta = theta, nugget, v)
    C_cross <- deepgp:::MaternSep(x_star, x, 1, theta, 0, v)
    C_new <- rep(1 + nugget, times = nrow(x_star))
    C_inv <- deepgp:::invdet(C)$Mi
    quadterm <- C_cross %*% C_inv %*% t(C_cross)
    mu[i,] <- C_cross %*% C_inv %*% Y
    s2 <- tau2 * (C_new - diag(quadterm))
    s2_sum <- s2_sum + s2
    #sigma_sum <- sigma_sum + Sigma
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("done")
  
  result <- list(
    mean = colMeans(mu),
    s2 = s2_sum / iterations + diag(cov(mu))
  )
  return(result)
}

########## DGP #############
###### Training ###########
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
    return( (- 0.5 * n) * log(quadterm) - 0.5 * log(det(R)) )
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
      U <- runif(1, 0, 1)
      while (accept == FALSE){
        count <- count + 1
        w[,i] <- w_prev * cos(theta) + nu * sin(theta) # Proposal
        dw <- deepgp:::sq_dist(w)
        #log_alpha <- loglik_yw(Y, w[,i], ls_y, g, v) - ll_prev # log-alpha
        log_alpha <- loglik_yw(Y, w, ls_y, g, v) - ll_prev # log-alpha
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
###### Prediction ########
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
  #Sigma <- vector("list", iterations)  # To store the correlation matrix, is a list, contain T matrices, each have size: 
  s2_sum <- rep(0, times = nrow(x_star))
  
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
    #R_new <- deepgp:::Matern(dw_new, 1, theta, nugget, v)
    R_new <- rep(1 + nugget, times = nrow(x_star))
    R_inv <- deepgp:::invdet(R)$Mi
    quadterm <- R_cross %*% R_inv %*% t(R_cross)
    mu[i,] <- R_cross %*% R_inv %*% Y
    
    #Sigma[[i]] <- scale_sample[i] * (R_new - quadterm)
    s2 <- scale_sample[i] * (R_new - diag(quadterm))
    s2_sum <- s2_sum + s2
    
    
    setTxtProgressBar(pb, i)
  }
  
  result <- list(
    mean = colMeans(mu),
    #sigma2 = (1/iterations) * Reduce(`+`, Sigma) + cov(mu)
    s2 = s2_sum / iterations + diag(cov(mu))
  )
  
  return(result)
}



gp1 <- GPemulator_MH_training(X_gp, Y_gp, n_iteration = 10000, burn_in = 5000)
pre1 <- GPprediction(gp1, X_gp, Y_gp, test_x, 1e-6)
plot_result_1(pre1, x_train = X_gp, y_train = Y_gp, x_test = test_x, y_test = test_y)

dgp1 <- fit_two_layer3.0_matern(x = X_gp, Y = Y_gp, u = 2, l = 1, ls_y = 0.1, ls_w = 0.1, node = 1,
                                W = X_gp, n_iteration = 10000, burn_in = 7000, nugget = 1e-6, v = 2.5)
pre2 <- Two_layer_prediction_matern(dgp1, X_gp, Y_gp, test_x)
plot_result_1(pre2, X_gp, Y_gp, test_x, test_y)
plot_ESS_samples1(dgp1$w_samples, X_gp)


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
      #log_l[i] <-logdet - (n * 0.5) * log(wtUUtw)
      log_l[i] <- logdet - 0.5 * wtUUtw
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
      U <- runif(1, 0, 1)
      while (accept == FALSE){
        count <- count + 1
        w[,i] <- w_prev * cos(theta) + nu * sin(theta) # Proposal
        dw <- deepgp:::sq_dist(w)
        log_alpha <- loglik_yw(Y, w, ls_y, g, v, k = k, Vecchia = Vecchia)[[1]] - ll_prev # log-alpha
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













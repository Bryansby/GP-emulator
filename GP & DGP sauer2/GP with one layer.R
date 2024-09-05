###### Training ############
GPemulator_MH_training <- function(x, Y, n_iteration = 5000, burn_in = 3000, ls_y = matrix(rep(0.1, ncol(x)), nrow = 1), nugget = 1e-8, g = deepgp:::eps, v = 2.5){
  dim <- ncol(x)
  n <- nrow(Y)
  ll_store <- rep(NA, length = n_iteration)
  
  loglik_xy <- function(x, y, ls_y, g = nugget){
    #' @description To compute the log-likelihood of the GP emulator
    n <- nrow(y)
    R <- deepgp:::Exp2Sep(x, x, 1, ls_y, g)
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
    R <- deepgp:::Exp2Sep(x, x, tau2 = 1, theta = theta_y_samples[i,], g = nugget)
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

gp <- GPemulator_MH_training(x_tray, y_tray, n_iteration = 10000, burn_in = 1)
fit <- deepgp:::fit_one_layer(x_tray, c(y_tray), sep = TRUE, true_g = 1e-8)
deepgp:::plot.gp(fit)

######## Prediction ###########
GPprediction <- function(list, x, Y, x_star, nugget){
  theta_y_samples <- list[[5]]  # size : iteration * dim
  scale_samples <- list[[6]]    # size : iteration * 1
  dim <- list[[7]]
  iterations <- length(theta_y_samples)
  
  cat("Predicting ... \n")
  pb <- txtProgressBar(min = 0, max = iterations, style = 3)

  mu <- matrix(NA, nrow = iterations, ncol = nrow(x_star))  # To store the mean, size: T * m'
  Sigma <- vector("list", iterations)  # To store the correlation matrix, is a list, contain T matrices, each have size: 
  
  for(i in 1:iterations){
    dx <- distance(x)
    d_new <- distance(x_star)
    d_cross <- distance(x_star, x)
    theta <- theta_y_samples[i,]
    C <- deepgp:::Matern(dx, 1, theta, g, v)
    C_cross <- deepgp:::Matern(d_cross, 1, theta, nugget, v)
    C_new <- deepgp:::Matern(d_new, 1, theta, g, v)
    C_inv <- deepgp:::invdet(C)$Mi
    L <- chol(C)
    Z <- forwardsolve(t(L), t(C_cross))
    quadterm <- t(Z) %*% Z
    mu[i,] <- C_cross %*% C_inv %*% Y
    sigma_w <- (C_new - quadterm)
    
    # Update progress bar
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  cat("done")
  
  return(colMeans(mu))
}

pre <- GPprediction(gp, X_gp, Y_gp, test_x, 1e-8)

plot(x_test, pre)

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

####################
# Function for Maximum-Minimum Ordering
max_min_ordering <- function(data) {
  n <- nrow(data)
  # Initialize the ordered indices vector
  ordering <- integer(n)
  
  # Start by selecting a random point
  ordering[1] <- sample(1:n, 1)
  
  # Track the distances from selected points
  dist_to_selected <- rep(Inf, n)
  
  # Function to compute Euclidean distance between two points
  euclidean_dist <- function(x, y) {
    sqrt(sum((x - y)^2))
  }
  
  for (i in 2:n) {
    # Update the minimum distances to the already selected points
    for (j in 1:n) {
      dist_to_selected[j] <- pmin(dist_to_selected[j], euclidean_dist(data[ordering[i-1], ], data[j, ]))
    }
    
    # Select the point with the maximum minimum distance
    ordering[i] <- which.max(dist_to_selected)
  }
  
  return(ordering)
}

# Example usage
set.seed(123)

# Create some sample data (e.g., 2D points)
data <- matrix(runif(20), ncol = 1)

# Get the ordering
ordering <- max_min_ordering(data)

# Print the ordering of points
print(ordering)

# Optionally, plot the points and their order
plot(data, pch = 19, col = "blue", main = "Maximum-Minimum Ordering")
text(data, labels = ordering, pos = 3, col = "red")

###############

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







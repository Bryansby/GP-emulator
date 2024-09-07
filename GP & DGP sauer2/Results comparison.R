################ One-D example #################
library(plgp)
library(MASS)
library(gridExtra)
library(ggplot2)
library(deepgp)
library(mvtnorm)
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



# Initialize an empty data frame to store k, k_2, and the corresponding RMSE
results <- data.frame(k = integer(), k_2 = integer(), RMSE = numeric())

# Loop through values of k and k_2 from 5 to 18
for (k_value in 5:18) {
  # Run the model for each value of k and k_2
  dgp_Vecchia_RO <- fit_two_layer3.0_matern_Vecchia_RO(
    x = x_train, Y = y_train, ls_y = 0.1, ls_w = 0.1, node = 1,
    W = x_train, v = 2.5, k = k_value, k_2 = k_value
  )
  
  # Make predictions using the model
  pre_Vec_RO <- Two_layer_prediction_matern(dgp_Vecchia_RO, x_train, y_train, x_test)
  
  # Calculate RMSE for the current values of k and k_2
  current_rmse <- rmse(c(y_test), c(pre_Vec_RO$mean))
  
  # Append the current values of k, k_2, and RMSE to the results data frame
  results <- rbind(results, data.frame(k = k_value, k_2 = k_value, RMSE = current_rmse))
}

# Print or view the results
print(results)


# Initialize an empty data frame to store k, k_2, and the corresponding RMSE
results2 <- data.frame(k = integer(), k_2 = integer(), RMSE = numeric())

# Loop through values of k and k_2 from 5 to 18
for (k_value in 5:18) {
  # Run the model for each value of k and k_2
  dgp_Vecchia_MM <- fit_two_layer3.0_matern_Vecchia_MM(
    x = x_train, Y = y_train, ls_y = 0.1, ls_w = 0.1, node = 1,
    W = x_train, v = 2.5, k = k_value, k_2 = k_value
  )
  
  # Make predictions using the model
  pre_Vec_MM <- Two_layer_prediction_matern(dgp_Vecchia_MM, x_train, y_train, x_test)
  
  # Calculate RMSE for the current values of k and k_2
  current_rmse <- rmse(c(y_test), c(pre_Vec_MM$mean))
  
  # Append the current values of k, k_2, and RMSE to the results data frame
  results2 <- rbind(results2, data.frame(k = k_value, k_2 = k_value, RMSE = current_rmse))
}

# Print or view the results
print(results2)

p1 <- ggplot(data = results, aes(x = k, y = RMSE)) + geom_line() + labs(title = "RMSE by NN", x = "Number of neaerst neighbors", y = "RMSE")
p2 <- ggplot(data = results2, aes(x = k, y = RMSE)) + geom_line() + labs(title = "RMSE by MM", x = "Number of neaerst neighbors", y = "RMSE")
grid.arrange(p1, p2, ncol = 2)


###############################
rmse_dgpsi <- rep(NA, 1000)
rmse_dgpsi_vec <- rep(NA, 1000)
for (i in 1:1000) {
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
dgpsi_rmse_df

df <- data.frame(
  number_of_experiment = c(1:1000),
  rmse_dgpsi = rmse_dgpsi,
  rmse_dgpsi_vec = rmse_dgpsi_vec
)


library(tidyr)
library(dplyr)
library(ggplot2)
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
#############################






dgpsi_rmse_df <- data.frame(
  Model = c("DGP", 
            "DGP using Vecchia"),  # Labels for each model (optional)
  RMSE = c(rmse_dgpsi, rmse_dgpsi_vec)  # RMSE values
)
dgpsi_rmse_df

library(ggplot2)
library(patchwork)
library(dgpsi)
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
set.seed(999)

############# Training points ##############
X_gp <- seq(0, 1, length = 15)
Y_gp <- sapply(X_gp, f123)
############ Prediction points ############
test_x <- as.matrix(seq(0, 1, length = 300), ncol = 1)
test_y <- as.matrix(sapply(test_x, f123), ncol = 1)

##### Train using dgpsi ######
m_gp <- gp(X_gp, Y_gp, name = 'sexp') 
nugget <- m_gp$specs$nugget
lengthscales <- m_gp$specs$lengthscales
scale <- m_gp$specs$scale

###### Prediction using dgpsi #######
p <- predict(m_gp, x = test_x)
mu <- p$results$mean
sd <- sqrt(p$results$var)

###### Prediction using gp emulator #######
gp <- GPemulator(x = as.matrix(X_gp), Y = as.matrix(Y_gp), 
                 ls = lengthscales, nugget = nugget, kernel_function = "sexp", 
                 scale = scale, x_star = as.matrix(test_x))
mu1 <- gp$mean
testr_std <- gp$sd

############ Prediction using Vecchia ##########
gpv <- GPvecchia(as.matrix(X_gp), as.matrix(Y_gp), scale = scale, ls = lengthscales, nugget = nugget, x_star = as.matrix(test_x), n = 5)

############### Draw graph to do the comparison ###########
plot(test_x, test_y, col = 'red', pch = 19)
points(test_x, mu, col = 'blue', pch = 19) # dgpsi
points(test_x, mu1, col = 'darkgreen', pch = 19) # Mine emulator
points(gpv$x, gpv$mean, col = 'orange', pch = 19) # using Vecchia approximation 
points(X_gp, Y_gp, col = 'green', pch = 20)

lines(test_x, mu + sd, col = 'blue', lty = 'dashed', lwd = 2)
lines(test_x, mu - sd, col = 'blue', lty = 'dashed', lwd = 2)
lines(test_x, mu1 + testr_std, col = 'darkgreen', lty = 'dashed', lwd = 2)
lines(test_x, mu1 - testr_std, col = 'darkgreen', lty = 'dashed', lwd = 2)
points(gpv$x, gpv$mean + 2 *gpv$sd, col = 'green', pch = 19, cex = 0.3)
points(gpv$x, gpv$mean - 2 * gpv$sd, col = 'green', pch = 19, cex = 0.3)

## M-spline calculation using R package splines2

# Generate a data set of size n=200 as in the previous example
neg.log.u <- -log(runif(200))
y <- (neg.log.u)^(1/3)

# Select 4 interior knots, and estimate the theta's
library(splines2)
int.kn <- quantile(y, seq(0.1, 0.9, length.out=4))
bound.kn <- c(0, max(y) + 1e-4)
psi <- mSpline(y, knots = int.kn, Boundary.knots = bound.kn, degree = 3)
target_func <- function(theta, psi, t){
  diff = sum((3*t^2 - psi %*% (theta))^2)
  return(diff)
}

opts <- list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-4, "maxeval" = 100)

est_theta = nloptr(x0 = matrix(1, nrow = 7), eval_f = target_func,
                   lb = (rep(0, 7)), ub = rep(Inf, 7), 
                   opts = opts, psi = psi, t = y)$solution
est_theta

# Plot the true and estimated baseline hazard and cumulative baseline hazard functions
v <- seq(from = 0, to = 1.5, by = 0.01)
psi_v = mSpline(v, knots = int.kn, Boundary.knots = bound.kn, degree = 3)
plot(3*v^2 ~ v, type = "l")
est_h0t = psi_v %*% est_theta
lines(est_h0t ~ v, type = "l", col = "red")

Psi_v = mSpline(v, knots = int.kn, Boundary.knots = bound.kn, degree = 3, integral = TRUE)
plot(v^3 ~ v, type = "l")
est_H0t = Psi_v %*% est_theta
lines(est_H0t ~ v, type = "l", col = "red")


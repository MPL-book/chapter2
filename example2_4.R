## Comparing computational burden

# Generate data
neg.log.u <- -log(runif(200))
y <- (neg.log.u)^(1/3)

int.kn <- quantile(y, seq(0.1, 0.9, length.out=4))
bound.kn <- c(0, max(y) + 1e-4)
psi <- mSpline(y, knots = int.kn, Boundary.knots = bound.kn, degree = 3)

# Estimate theta to use for matrix multiplication

target_func <- function(theta, psi, t){
  diff = sum((3*t^2 - psi %*% (theta))^2)
  return(diff)
}

opts <- list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-4, "maxeval" = 100)

est_theta = nloptr(x0 = matrix(1, nrow = 7), eval_f = target_func,
                   lb = (rep(0, 7)), ub = rep(Inf, 7), 
                   opts = opts, psi = psi, t = y)$solution


# Create a grid of times to estimate the cumulative baseline hazard for
v <- seq(from = 0, to = 1.5, by = 0.01)
Psi_v = mSpline(v, knots = int.kn, Boundary.knots = bound.kn, degree = 3, integral = TRUE)

# Record the time taken to estimate the cumulative baseline hazard
start.time <- Sys.time()
est_H0t = Psi_v %*% est_theta
end.time <- Sys.time()
time.taken <- round(end.time - start.time,5)

# Estimate theta to use in exponential transformation

target_func_exp <- function(theta, psi, t){
  alpha = psi %*% (theta)
  diff = sum((3*t^2 - exp(alpha))^2)
  return(diff)
}

opts <- list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-4, "maxeval" = 100)

est_theta_alpha = nloptr(x0 = matrix(1, nrow = 7), eval_f = target_func,
                         lb = (rep(-Inf, 7)), ub = rep(Inf, 7), 
                         opts = opts, psi = psi, t = y)$solution


# Create functions which allow us to integrate the baseline hazard to get the cumulative baseline hazard

h0_t <- function(t, theta, int.kn, bound.kn){
  psi <- bSpline(t, knots = int.kn, Boundary.knots = bound.kn, degree = 3)
  h0 <- exp(psi %*% est_theta_alpha)
}

H0_t <- function(upper_t, f, i, b){
  integrate(f, 0, upper_t, int.kn = i, bound.kn = b)$value
}

# Record the time taken to estimate the cumulative baseline hazard
start.time <- Sys.time()
est_theta_H <- sapply(seq(from = 0, to = 1.5, by = 0.01), H0_t, f = h0_t, i = int.kn, b = bound.kn)
est_H0t_exp = Psi_v %*% est_theta
end.time <- Sys.time()
time.taken <- round(end.time - start.time, 5)


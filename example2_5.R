## Gaussian basis and cumulative basis function

# Generate data as per the previous examples
neg.log.u <- -log(runif(200))
y <- (neg.log.u)^(1/3)

# Write functions to compute Gaussian basis functions
sigma <- 0.5
kn <- quantile(y, seq(0.01, 0.99, length.out=5))
phi.gauss <- NULL

phi_gauss_f <- function(y, kn, a, b){
  phi.gauss <- NULL
  for(u in 1:length(kn)){
    delta_u = pnorm((b - kn[u])/sigma) - pnorm((a - kn[u])/sigma)
    phi.gauss = cbind(phi.gauss, (1/(sigma*delta_u)) * dnorm((y - kn[u])/sigma))
  }
  return(phi.gauss)
}

Phi_gauss_f <- function(y, kn, a, b){
  Phi.gauss <- NULL
  for(u in 1:length(kn)){
    delta_u = pnorm((b - kn[u])/sigma) - pnorm((a - kn[u])/sigma)
    Phi.gauss = cbind(Phi.gauss, (1/delta_u) * (pnorm((y - kn[u])/sigma) - pnorm((a - kn[u])/sigma)))
  }
  return(Phi.gauss)
}

phi.gauss <- phi_gauss_f(y, kn, a = 0, b = max(y))
Phi.gauss <- Phi_gauss_f(y, kn, a = 0, b = max(y))

# Estimate theta

target_func_gauss <- function(theta, phi.gauss, t){
  diff = sum((3*t^2 - phi.gauss %*% (theta))^2)
  return(diff)
}

library(nloptr)
opts <- list("algorithm"="NLOPT_LN_COBYLA",
             "xtol_rel"=1.0e-4, "maxeval" = 1000)

est_theta = nloptr(x0 = matrix(1, nrow = 5), eval_f = target_func_gauss,
                   lb = (rep(0, 5)), ub = rep(Inf, 5), 
                   opts = opts, phi.gauss = phi.gauss, t = y)$solution

# Plot baseline hazard function
v <- seq(from = 0, to = 1.5, by = 0.01)
plot(3*v^2 ~ v, type = "l",                         
     main = "Hazard and Gaussian basis approx.")
psi_v <- phi_gauss_f(v, kn, a = 0, b = max(y))
est_h0t = psi_v %*% est_theta
lines(est_h0t ~ v, type = "l", lty=2)
legend(0, 6.8, legend=c("true", "approx"), lty=1:2, cex=0.9)

# Plot cumulative baseline hazard function
Psi_v = Phi_gauss_f(v, kn, a = 0, b = max(y))
plot(v^3 ~ v, type = "l",                   
     main = "Cumulative hazard and Gaussian basis approx.")
est_H0t = Psi_v %*% est_theta
lines(est_H0t ~ v, type = "l", lty=2)
legend(0, 3.4, legend=c("true", "approx"), lty=1:2, cex=0.9)

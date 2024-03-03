## A simulation study

# Write simulate() function to repeatedly generate a sample
simulate <- function(n, pi_E, g1, g2){
  
  # Simulate covariates
  x1 <- runif(n, 0, 1)
  x2 <- 5 * runif(n, 0, 1)
  X <- cbind(x1, x2)
  beta_true <- c(0.75, -0.5)
  eXtB <- exp(X %*% beta_true)
  
  # Simulate true event times
  neg.log.u <- - log(runif(n, 0, 1))
  y_i <- (neg.log.u/(eXtB))^(1/3)
  
  # Simulate censoring types & times
  u_E <- runif(n, 0, 1)
  u_L <- runif(n, 0, 1)
  u_R <- runif(n, u_L, 1)
  
  t_L <- (y_i^as.numeric(u_E < pi_E)) * (g1 * u_L)^as.numeric(pi_E < u_E & 
                            g1*u_L <= y_i & y_i <= g2*u_R) * (g2 * u_R)^as.numeric(pi_E < u_E & g2*u_R < y_i) * 0^as.numeric(pi_E < u_E & y_i < g1*u_L) 
  
  t_R <- (y_i^as.numeric(u_E < pi_E)) * (g1 * u_L)^as.numeric(pi_E < u_E & y_i < g1*u_L) * (g2 * u_R)^as.numeric(pi_E < u_E 
                          & g1*u_L <= y_i & y_i <= g2*u_R) * Inf^(pi_E < u_E & g2*u_R < y_i)
  
  # Compute midpoint and censoring type for midpoint imputation
  
  midpoint <- t_L + (t_R - t_L)/2
  event <- rep(1, n)
  event[which(is.infinite(midpoint))] <- 0
  midpoint[which(is.infinite(midpoint))] <- t_L[which(is.infinite(midpoint))]
  
  df = data.frame(t_L, t_R, x1, x2, midpoint, event)
  return(df)
}

# Run a simulation study

ch2_ex26_save = matrix(0, nrow = 100, ncol = 8)
baseline_h0 = NULL
library(survivalMPL)

for(s in 1:100){
  
  # Generate data set
  dat <- simulate(n=500, pi_E = 0.5, g1 = 0.9, g2 = 1.3)
  
  # Fit MPL model
  fit.mpl <- coxph_mpl(Surv(t_L, t_R, type = "interval2") ~ 
                         x1 + x2, data = dat, 
                       basis = "m", n.knots = c(5,0))
  
  
  # Save coefficient estimates & SE
  ch2_ex26_save[s,1] <- fit.mpl$coef$Beta[1]
  ch2_ex26_save[s,2] <- fit.mpl$coef$Beta[2]
  
  ch2_ex26_save[s,3] <- fit.mpl$se$Beta$M2HM2[1]
  ch2_ex26_save[s,4] <- fit.mpl$se$Beta$M2HM2[2]
  
  # Fit PH model
  fit.ph <- coxph(Surv(midpoint, event) ~ x1 + x2, data = dat)
  
  # Save coefficient estimates & SE
  ch2_ex26_save[s,5] <- fit.ph$coefficients[1]
  ch2_ex26_save[s,6] <- fit.ph$coefficients[2]
  
  ch2_ex26_save[s,7] <- sqrt(diag(fit.ph$var))[1]
  ch2_ex26_save[s,8] <- sqrt(diag(fit.ph$var))[2]
  
  # Save baseline hazard times and estimates
  psi <- mSpline(survfit(fit.ph)$time, 
                 knots = fit.mpl$knots$Alpha[2:6],
                 Boundary.knots = c(fit.mpl$knots$Alpha[1], 
                                    fit.mpl$knots$Alpha[7]))
  h0t_MPL = psi %*% fit.mpl$coef$Theta
  
  baseline_h0 <- cbind(baseline_h0, 
                       basehaz(fit.ph)$time, basehaz(fit.ph)$hazard, h0t_MPL)
  
}




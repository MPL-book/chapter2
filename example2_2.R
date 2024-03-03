## Piece-wise constant approximation to hazard

# Generate a random sample of n=200 from this
# distribution
neg.log.u <- -log(runif(200))
y <- (neg.log.u)^(1/3)
# Select bins using equal bin counts where
# m=10, 30, 100
bins_f <- function(m, y) {
  bins <- quantile(y, probs = seq(0, 1, length.out = (m+1)))
  
  theta <- cum_theta <- rep(0, m)
  
  for(u in 1:m){
    bin_u_obs <- y[which(bins[u] <= y & y < bins[u+1])]
    #mean hazards in the bin
    theta[u] <- mean(3*bin_u_obs^2)
    cum_theta[u] = mean(bin_u_obs^3)
  }
  
  out <- list(theta = theta, cum_theta = cum_theta, 
              bins = bins)
  return(out)
}

h0t_m10 <- bins_f(10, y)
h0t_m30 <- bins_f(30, y)
h0t_m100 <- bins_f(100, y)

# Plot the true and approximated baseline
# hazard function for the different values of m
v <- seq(from = min(y), to = max(y), by = 0.01)
plot(3 * v^2 ~ v, type = "l", col = "black", 
     main = "hazard and approxmations")
lines(c(h0t_m10$theta, h0t_m10$theta[10]) ~ 
        h0t_m10$bins, type = "s", lty=2)
lines(c(h0t_m30$theta, h0t_m30$theta[30]) ~ 
        h0t_m30$bins, type = "s", lty=3)
lines(c(h0t_m100$theta, h0t_m100$theta[100]) ~ 
        h0t_m100$bins, type = "s", lty=4)
legend("topleft", legend=c("true", "m=10", "m=30", "m=100"), 
       lty=1:4, cex=0.7)

# Plot the true and approximated cumulative
# baseline hazard functions

plot(v^3 ~ v, type = "l", 
     col = "black", main = "cumulative hazard and approxmations")
lines(c(h0t_m10$cum_theta, h0t_m10$cum_theta[10]) ~ 
        h0t_m10$bins, type = "s", lty=2)
lines(c(h0t_m30$cum_theta, h0t_m30$cum_theta[30]) ~ 
        h0t_m30$bins, type = "s", lty=3)
lines(c(h0t_m100$cum_theta, h0t_m100$cum_theta[100]) ~ 
        h0t_m100$bins, type = "s", lty=4)
legend("topleft", legend=c("true", "m=10", "m=30", "m=100"), 
       lty=1:4, cex=0.7)


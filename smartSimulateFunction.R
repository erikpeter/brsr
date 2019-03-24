smart.simulate <- function(n.sim = 1e4, assumptions, alpha, truth) {
  target.power <- assumptions$target.power
  delta.est <- assumptions$delta.est
  sd.initial.est <- assumptions$sd.initial.est
  
  n.per.group <- ceiling(power.t.test(delta = delta.est, sd = sd.initial.est, power = target.power, sig.level = alpha)$n)
  interim.n.per.group <- ceiling(0.5*n.per.group)
  n <- dim(truth)[1]
  blinded.sd.est <- matrix(0, nrow = n.sim, ncol = n)
  new.n.per.group <- matrix(0, nrow = n.sim, ncol = n)
  additional.n.per.group <- matrix(0, nrow = n.sim, ncol = n)
  pval <- matrix(0, nrow = n.sim, ncol = n)
  tval <- matrix(0, nrow = n.sim, ncol = n)
  
  for (i in 1:n.sim) {
    
    
    x.raw <- rnorm(interim.n.per.group)
    y.raw <- rnorm(interim.n.per.group)
    
    # reestimate sample size
    for (j in 1:n) {
      x <- truth$sd[j]*x.raw
      y <- truth$sd[j]*y.raw + truth$delta[j]
      blinded.sd.est[i, j] <- sd(c(x, y))
      nnpg <- ceiling(power.t.test(delta = delta.est, sd = blinded.sd.est[i, j], power = target.power, sig.level = alpha)$n)
      new.n.per.group[i, j] <- max(interim.n.per.group, nnpg)
      additional.n.per.group[i, j] <- new.n.per.group[i, j] - interim.n.per.group
    }
    
    # reestimate sample size
    max.additional <- max(additional.n.per.group[i, ])
    x.raw2 <- rnorm(max.additional)
    y.raw2 <- rnorm(max.additional)
    
    for (j in 1:n) {
      final.x <- c(x.raw, x.raw2) * truth$sd[j]
      final.y <- c(y.raw, y.raw2) * truth$sd[j] + truth$delta[j]
      final.x <- final.x[1:new.n.per.group[i, j]]
      final.y <- final.y[1:new.n.per.group[i, j]]
      test.result <- t.test(final.x, final.y, var.equal = TRUE)
      pval[i, j] <- test.result$p.value
      tval[i, j] <- test.result$statistic
    }
    
    
  }
  
  # summarize results

  power <- mean.sample <- rep(0, n)

  for (j in 1:n) {
    power[j] <- sum(pval[, j] < alpha) / dim(pval)[1]
    mean.sample[j] <- mean(new.n.per.group[, j])
  }
  return(cbind(truth, data.frame(power, mean.sample)))
}

assumptions <- data.frame(target.power = 0.9, delta.est = 0.5, sd.initial.est = 1)
truth <- expand.grid(sd = c(0.5, 0.8, 0.9, 1, 1.2, 1.5, 2.0), delta = c(0, 0.25, 0.5, 0.6, 0.8, 1))
  
  
test.result <- smart.simulate(1e6, assumptions, 0.05, truth)
print(test.result)





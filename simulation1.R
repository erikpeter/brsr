# simple t.test blinded reestimate

target.power <- 0.9
delta.est <- 0.5
sd.initial.est <- 1
alpha <- 0.05

n.per.group <- ceiling(power.t.test(delta = delta.est, sd = sd.initial.est, power = target.power, sig.level = alpha)$n)
interim.n.per.group <- ceiling(0.5*n.per.group)

simulate.power <- function(n.sim = 1e5, true.sd = 1, true.delta = 0.5) {
  blinded.sd.est <- new.n.per.group <- additional.n.per.group <- rep(0, n.sim) 
  tval <- pval <- rep(0, n.sim)
  
  for (i in 1:n.sim) {
    
    # simulate data
    x <- rnorm(interim.n.per.group, 0, true.sd)
    y <- rnorm(interim.n.per.group, true.delta, true.sd)
    
    # reestimate sample size
    blinded.sd.est[i] <- sd(c(x, y))
    new.n.per.group[i] <- ceiling(power.t.test(delta = delta.est, sd = blinded.sd.est[i], power = target.power, sig.level = alpha)$n)
    additional.n.per.group[i] <- max(0, new.n.per.group[i] - interim.n.per.group)
    
    # add new samples & run final analysis
    final.x <- c(x, rnorm(additional.n.per.group[i], 0, true.sd))
    final.y <- c(y, rnorm(additional.n.per.group[i], true.delta, true.sd))
    test.result <- t.test(final.x, final.y, var.equal = TRUE)
    pval[i] <- test.result$p.value
    tval[i] <- test.result$statistic
  }
  
  power <- sum(pval < alpha)/length(pval)
  return(power)
}

print(simulate.power(n.sim = 100))

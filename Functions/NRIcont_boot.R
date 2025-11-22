##### NET RECLASSIFICATIONS IMROVEMENT INDEX BOOT Function ####
# Returns vector with bootstrapped continuous NRI, event-NRI and non-event-NRI from sample data

NRI.cont.boot <- function (NRIdata, runs) {
  boot <- boot(data = NRIdata, statistic = NRIcont, R = runs)
  mean <- boot$t0[1]
  event <- boot$t0[2]
  non.event <- boot$t0[3]
  sd <- sd(boot$t[,1])
  margin <- qt(0.975, df = runs-1)*sd/sqrt(runs)
  ci.lower <- mean - margin
  ci.upper <- mean + margin
  z <- sqrt(runs) * boot$t0[1] / sd
  pval <- 2*pnorm(z, mean, sd, lower.tail = F)
  
  return(c(mean, event, non.event, ci.lower, ci.upper, pval))
}
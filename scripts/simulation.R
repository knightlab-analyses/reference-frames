library(tidyverse)
library(driver) # devtools::install_github("jsilve24/driver")
library(broom)

setwd("~/Research/scale_invariant_analysis/results/2018-07-27_simulation_2/")

set.seed(5)


# helper functions --------------------------------------------------------

# multivariate normal sampler
# 
# n is number of samples
# m is mean - p-vector
# C is covariance
# output: p x n matrix 
rmvnorm <- function(n, m, C){
  # basic input validation
  stopifnot(length(m) == nrow(C))
  stopifnot(nrow(C)==ncol(C))
  
  # main function
  p <- length(m)
  Z <- matrix(rnorm(n*p), p, n)
  L <- t(chol(C))
  Z <- L %*% Z
  Z <- sweep(Z, 1, m, FUN=`+`)
  return(Z)
}

# Start Simulation --------------------------------------------------------

# Basic parameters Lots of samples so no identifiability issues (for simplicity)
# before an after treatment
Nbefore <- 1000 # number of samples in 
Nafter <- 1000
D <- 100



sim <- function(mu){
  # Simulate Truth (log normal with no correlation for simplicity)
  # B represents the true differential abundance with an intercept
  B <- matrix(rnorm(2*D, mean = mu, sd=5), D, 2) # True differential abundance given by B[,2]
  X <- rbind(1, c(rep(0, Nbefore), rep(1, Nafter)))
  Sigma <- diag(rnorm(D, mean=0, sd=4)^2)
  Y <- rmvnorm(Nbefore+Nafter, rep(0, D), Sigma)
  Y <- exp(Y + B%*%X)  # True abundance data Y ~ lognormal(BX, Sigma)
  
  # Convert to Proportions
  Yprop <- t(miniclo(t(Y)))
  
  
  
  # Compare 
  fit <- lm(t(log(Yprop)) ~ t(X)-1)
  # fit <- lm(alr(t(Yprop)) ~ t(X)-1)
  de.fit <- coef(fit)[2,]
  # de.fit <- c(log(alrInv(de.fit)))
  
  return(cbind(B[,2], de.fit))
}

mu <- seq(-8, 8, by=4)
fit <- lapply(as.list(mu), sim)

lse <- map(fit, ~log(sum(exp(.x[,1]))))

map(fit, ~log(sum(exp(.x[,2]))))

dat <- map(fit, as.data.frame) %>%
  map(~mutate(.x, rankV1 = rank(V1), rankde.fit=rank(de.fit))) %>% 
  bind_rows(.id="Mean DE") 


colfunc <- colorRampPalette(c("blue", "red"))
cols <- colfunc(length(unique(dat[,1])))


# Plot results
pdf(file="plot.pdf", width=8, height=6)

plot(0, 0, xlim = c(-20, 20), ylim=c(-20, 20), col="white", asp=1,
     ylab="Inferred Differential Expression", 
     xlab="True Differential Expression")
lines(-100:100, rep(0, 201), col="grey")
lines(rep(0, 201), -100:100, col="grey")
lines(-100:100, -100:100, col="black", lwd=2)
points(dat[,"V1"], dat[,"de.fit"], col=cols[as.integer(dat[,"Mean DE"])], pch=20)
text(-23, 10, "Falsely Positive DE", col="darkgrey")
text(-23, -10, "Truely Negative DE", col="darkgrey")
text(23, 10, "Truely Positive DE", col="darkgrey")
text(23, -10, "Falsely Negative DE", col="darkgrey")
text(20, 18, "Optimal", srt=45, col="black")
legend("right", legend=signif(unlist(lse), 2), 
       fill=cols, title="LSE(True DE)", xpd=TRUE)

dev.off()


# Ranks
pdf(file="plot.ranks.pdf", width=8, height=6)

plot(0, 0, xlim = c(0, D+20), ylim=c(0, D+20), col="white", asp=1,
     ylab="Rank Inferred Differential Expression", 
     xlab="Rank True Differential Expression")
lines(0:(D+20), 0:(D+20), col="black", lwd=2)
points(dat[,"rankV1"], dat[,"rankde.fit"], col=cols[as.integer(dat[,"Mean DE"])], pch=20)
text(D+10, D+6, "Optimal", srt=45, col="black")
legend("right", legend=signif(unlist(lse), 2), 
       fill=cols, title="LSE(True DE)", xpd=TRUE)
dev.off()






# # all in one --------------------------------------------------------------
# 
# # Plot results
# pdf(file="plot.pdf", width=8, height=6)
# par(mfrow=c(2,1))
# plot(0, 0, xlim = c(-20, 20), ylim=c(-20, 20), col="white", asp=1,
#      ylab="Inferred DE", 
#      xlab="True DE")
# lines(-100:100, rep(0, 201), col="grey")
# lines(rep(0, 201), -100:100, col="grey")
# lines(-100:100, -100:100, col="black", lwd=2)
# points(dat[,"V1"], dat[,"de.fit"], col=cols[as.integer(dat[,"Mean DE"])], pch=20)
# text(-23, 10, "Falsely Positive DE", col="darkgrey")
# text(-23, -10, "Truely Negative DE", col="darkgrey")
# text(23, 10, "Truely Positive DE", col="darkgrey")
# text(23, -10, "Falsely Negative DE", col="darkgrey")
# text(20, 18, "Optimal", srt=45, col="black")
# legend("right", legend=signif(unlist(lse), 2), 
#        fill=cols, title="LSE(True DE)", xpd=TRUE)
# 
# 
# # Ranks
# plot(0, 0, xlim = c(0, D+20), ylim=c(0, D+20), col="white", asp=1,
#      ylab="Rank Inferred DE", 
#      xlab="Rank True DE")
# lines(0:(D+20), 0:(D+20), col="black", lwd=2)
# points(dat[,"rankV1"], dat[,"rankde.fit"], col=cols[as.integer(dat[,"Mean DE"])], pch=20)
# text(D+10, D+6, "Optimal", srt=45, col="black")
# legend("right", legend=signif(unlist(lse), 2), 
#        fill=cols, title="LSE(True DE)", xpd=TRUE)
# dev.off()






# Note:  I could have used ILR instead of log proportions - get same result here
# tYilr <- ilr(t(Yprop))
# fit <- lm(tYilr ~ t(X) -1)
# de.fit <- ilrInv(coef(fit)[2,])
# plot(log(de.fit), B[,2], asp=1)
# lines(-100:100, -100:100, col="red", lwd=2)
# lines(-100:100, rep(0, 201), col="grey")
# lines(rep(0, 201), -100:100,, col="grey")

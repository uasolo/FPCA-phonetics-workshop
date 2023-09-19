# Exercise on Functional PCA 
# Author: Michele Gubian
# Last revision: September 2020

# Landmark reg is about time warping and it's SOLELY based on the landmark locations.
# It is NOT based on the curve shapes.

# Data: landmarks extracted from 3 curves - but there ain't no curves here!

# Ex 1: only duration varies (linear time normalisation)
# relative position of mid landmark remains the same.
# (jitter is there because otherwise some matrix becomes singular.. )
land <- data.frame(l1 = rep(0,3), l2 = jitter(1:3, amount = .01), l3 = 2*(1:3))
reg <- landmarkreg.nocurve(land,nhknots = 8, hlambda=1e-8, wlambda =1e-8)
op <- par(mfrow=c(2,1))
plot(reg$warpfd, xlab = "reg. time", ylab = "time", main = "h(t)", las = 1)
plot(reg$logvelfd, xlab = "reg. time", ylab = "log rate", main = "log rate", las = 1)
par(op)
# check results by manual calculation

# Reconstruct durations based on logvel curves
tx <- seq(from = reg$land[1], to = reg$land[length(reg$land)], length.out = 30) # sampling the reg time axis
fdParObj <- fdPar(fdobj = reg$logvelfd$basis, Lfdobj = 2, lambda = 1e-6) # a typical fdPar object based on logvelfd basis
i <- 2 # pick a curve
Interval <- 1 # pick an interval
eval.fd(tx, reg$logvelfd[i]) %>%
  `*`(-1) %>%  exp %>% # integrand in inverse formula
  smooth.basis(tx, ., fdParObj) %>%
  .$fd %>%
  defint.fd(c(reg$land[Interval],reg$land[Interval+1])) # definite integral in inverse formula

# compare with land

# (In this case it might be quicker using exponentiate.fd, but this is more general,
# i.e. when you reconstruct logvelfd curves from emmeans results)


# Ex 2: only position of mid landmark varies
land <- data.frame(l1 = rep(0,3), l2 = c(2,5,8), l3 = jitter(rep(10,3), amount = .01))
# repeat the above (except for manual calculations)

# Ex 3: absolute position of mid landmark constant, total duration varies
land <- data.frame(l1 = rep(0,3), l2 = jitter(rep(2,3), amount = .01), l3 = c(5, 8, 10))
# repeat the above (except for manual calculations)

# Exercise on Functional PCA 
# Author: Michele Gubian
# Last revision: May 18th 2018

# ex9: fundamentals of landmark registration

# Landmark reg is about time warping and it's SOLELY based on the landmark locations.
# It is NOT based on the curve shapes.

# Data: landmarks extracted from 3 curves - but there ain't no curves here!

# only duration varies (linear time normalisation)
# (jitter is there because otherwise some matrix becomes singular.. )
land <- data.frame(l1 = rep(0,3), l2 = jitter(1:3, amount = .01), l3 = 2*(1:3))
reg <- landmarkreg.nocurve(land,nhknots = 8, hlambda=1e-8, wlambda =1e-8)
op <- par(mfrow=c(2,1))
plot(reg$warpfd, xlab = "reg. time", ylab = "time", main = "h(t)", las = 1)
plot(reg$logvelfd, xlab = "reg. time", ylab = "log rate", main = "log rate", las = 1)
par(op)
# check results by manual calculation

# only position of mid landmark varies
land <- data.frame(l1 = rep(0,3), l2 = c(2,5,8), l3 = jitter(rep(10,3), amount = .01))
# repeat the above plot (no manual calculations in this case)
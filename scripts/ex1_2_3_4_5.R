# Exercise on Functional PCA 
# Author: Michele Gubian
# Last revision: May 17th 2018

# ex1: centered hill random height from 0
# ex2: ex1 + random global mean (show global mean subtraction)
# ex3: randomly skewed hill, constant height (show PC1 shape)
# ex4: ex3 + random uncorrelated height (two PCs, show difference w/o y norm)
# ex5: ex4 but the two variations are totally correlated (one PC, show scatter with parabula)

# adjust path or setwd
source('header.R')
ex = 'ex1' # change according to ex number
ex_dir = file.path(data_dir,ex)

# load data

n_curves = 50 # always 50 curves in these exercises
time_list = list()
y_list = list()
duration = c()

mean_y = c() # use in ex2 only

for (i in 1:n_curves) {
	curve = read.table(file.path(ex_dir,paste0('curve',i,'.data')),header=T)
	time_list[[i]] = curve$time
	y_list[[i]] = curve$y
	# from ex2: try with/out the following two lines
	#mean_y = c(mean_y,mean(	curve$y ))
  #y_list[[i]] = y_list[[i]] - mean_y[i]

	duration = c(duration,tail(curve$time,1))
}

# inspect data (optional)
i = 1 # random
plot(time_list[[i]], y_list[[i]], las = 1)

Y_max <- 1
X_max <- max(duration)
plot(c(0,X_max), c(-Y_max,Y_max), type='n')
for (i in sample(n_curves, 10)) { # 1:n_curves
  lines(time_list[[i]], y_list[[i]], col = i)
}

# smoothing

# create common basis

# obligatory linear time normalisation
mean_dur = mean(duration) 
norm_range <- c(0,mean_dur)
n_knots = 8 # try many
lambda = 1e-8 # try many

#knots <- seq(0,mean_dur,length.out = n_knots)
Lfdobj <- 3 # 2 + order of derivative expected to be used. We might need velocity, thus order = 1
norder <- 2 + Lfdobj  # a fixed relation about B-splines
nbasis <- n_knots + norder - 2 # a fixed relation about B-splines
basis <- create.bspline.basis(norm_range, nbasis, norder) #, knots)

fdPar <- fdPar(basis, Lfdobj, lambda)

# visual inspection on random curve(s)

i = 10
t_norm = (time_list[[i]] / duration[i]) * mean_dur # linear time normalisation
curve_fd = smooth.basis(t_norm,y_list[[i]],fdPar)$fd
# plot
plot(t_norm,y_list[[i]],col = 'red', ylim=c(0,Y_max),las=1,xlab='lin. norm. time',ylab='y',main=paste('curve',i,', n_knots =',n_knots,'lambda =',lambda))
lines(curve_fd,lwd=2)

# parameters set to:
# n_knots = 8
# lambda = 1e-8

# smooth all curves
# smooth.basis() does not accept different time samples for different curves.
# Thus we create smooth curves one by one on the same basis, store the spline coefficients and compose an fd object at the end.

y_coefs = matrix(nrow = nbasis, ncol = n_curves)
for (i in 1:n_curves) {
	t_norm = (time_list[[i]] / duration[i]) * mean_dur # linear time normalisation
    	y_coefs[,i] = c(smooth.basis(t_norm,y_list[[i]],fdPar)$fd$coefs)
}
y_fd = fd(coef=y_coefs, basisobj=basis)

plot(y_fd[sample(n_curves, 10)])



  # FPCA

lambda_pca    <- lambda 
pcafdPar  <- fdPar(basis, 2, lambda_pca) # here Lfdobj = 2 since we don't need derivatives of PC curves.
y_pcafd <- pca.fd(y_fd, nharm=2, pcafdPar) # compute first 2 PCs


# plot PC curves
op <- par(mfrow=c(2,1))
plot.pca.fd(y_pcafd, nx=40)
par(op)


# plot PC scores s1 and s2
xyplot	(y_pcafd$scores[,2] ~  y_pcafd$scores[,1]  ,cex=1.5, 
 	xlab = list(label=expression(s[1]),cex=2), ylab= list(label=expression(s[2]),cex=2)
)


# FPCA-based curve reconstruction

# select a curve, use PCs and PC scores 
i = 1
s1 = y_pcafd$scores[i,1]
s2 = y_pcafd$scores[i,2]


plot(norm_range,c(-2*Y_max,2*Y_max),type='n',las=1,xlab='lin. norm. time',ylab='y',main=paste('curve',i,', reconstructionn'))

lines(y_fd[i],col = 'black',lwd=2) # original curve
lines(y_pcafd$harmonics[1],col='red',lwd=2) # PC1
# lines(y_pcafd$harmonics[2],col='orange',lwd=2) # PC2

lines(y_pcafd$meanfd,col='blue',lwd=2) # mean curve
lines(y_pcafd$meanfd+s1*y_pcafd$harmonics[1],col='green',lwd=2) # mean + s1 * PC1
# lines(y_pcafd$meanfd+s1*y_pcafd$harmonics[1]+s2*y_pcafd$harmonics[2],col='violet',lwd=2) # mean + s1 * PC1 + s2 * PC2

		



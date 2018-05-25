# Exercise on Functional PCA 
# Author: Michele Gubian
# Last revision: May 18th 2018

# ex7: 2D curves, two symmetric curves with correlated heights, when y is higher z is lower, y is a hill, z is a valley.
# ex8: 2D curves, y is a symmetric hill with variable height, z is a skewed valley, if y is high z is skewed to the right

source('header.R')
ex = 'ex8' # change according to ex number
ex_dir = file.path(data_dir,ex)

# load data

n_curves = 50
time_list = list()
y_list = list()
duration = c()
z_list = list() 


for (i in 1:n_curves) {
	curve = read.table(file.path(ex_dir,paste0('curve',i,'.data')),header=T)
	time_list[[i]] = curve$time
	y_list[[i]] = curve$y
	duration = c(duration,tail(curve$time,1))
	z_list[[i]] = curve$z
}


# some plotting parameters, adjust to your needs
T_min = 0.2; T_max = 0.3
Y_min = 1; Y_max = 2


# exploratory plots 
op = par(ask=T,mfrow=c(2,1))

for (i in 1:n_curves) {
	plot(time_list[[i]],y_list[[i]],xlim=c(0,T_max), ylim=c(0,Y_max),las=1,xlab='time',ylab='y')
	plot(time_list[[i]],z_list[[i]],xlim=c(0,T_max), ylim=c(-Y_max,0),las=1,xlab='time',ylab='z')
}
par(op)

# smoothing

# create common basis

# obligatory linear time normalisation
mean_dur = mean(duration) 

norm_range <- c(0,mean_dur)
n_knots = 8 # try many
lambda = 1e-8

#knots <- seq(0,mean_dur,length.out = n_knots)
Lfdobj <- 3 # 2 + order of derivative expected to be used. We might need velocity, thus order = 1
norder <- 2 + Lfdobj  # a fixed relation about B-splines
nbasis <- n_knots + norder - 2 # a fixed relation about B-splines
basis <- create.bspline.basis(norm_range, nbasis, norder) #, knots)
fdPar <- fdPar(basis, Lfdobj, lambda)

# visual inspection on random curve(s)
i = 1
t_norm = (time_list[[i]] / duration[i]) * mean_dur # linear time normalisation
curve_fd = smooth.basis(t_norm,y_list[[i]],fdPar)$fd
# plot
plot(t_norm,y_list[[i]],col = 'red', ylim=c(0,Y_max),las=1,xlab='lin. norm. time',ylab='y',main=paste('curve',i,', n_knots =',n_knots,'lambda =',lambda))
lines(curve_fd,lwd=2)
# repeat with z 

# parameters set to:
n_knots = 8
lambda = 1e-8

# smooth all curves
# yz_coefs has a further dimension with respect to y_coefs, because curves are 2-dim trajectories.
yz_coefs = array(dim = c(nbasis,n_curves,2))
for (i in 1:n_curves) {
	t_norm = (time_list[[i]] / duration[i]) * mean_dur # linear time normalisation
    	yz_coefs[,i,1] = c(smooth.basis(t_norm,y_list[[i]],fdPar)$fd$coefs)
    	yz_coefs[,i,2] = c(smooth.basis(t_norm,z_list[[i]],fdPar)$fd$coefs)
}
yz_fd = fd(coef=yz_coefs, basisobj=basis)
op = par(mfrow=c(2,1))
plot(yz_fd) # plot(yz_fd[sample(1:n_curves, 3)])
par(op)


# FPCA

lambda_pca    <- lambda
pcafdPar  <- fdPar(basis, 2, lambda_pca) # here Lfdobj = 2 since we don't need derivatives of PC curves.
yz_pcafd <- pca.fd(yz_fd, nharm=2, pcafdPar) # compute first 2 PCs

# plot PC curves
op <- par(mfrow=c(2,2))
plot.pca.fd(yz_pcafd, nx=40)
par(op)

# cycle plot
op <- par(mfrow=c(2,1))
plot.pca.fd(yz_pcafd, nx=40,cycle=T)
par(op)

plot.pca.fd(yz_pcafd,harm=1,pointplot=T,cycle=T)





# plot PC scores s1 and s2
xyplot	(yz_pcafd$scores[,2] ~  yz_pcafd$scores[,1]  , cex=1.5, 
	xlab = list(label=expression(s[1]),cex=2), ylab= list(label=expression(s[2]),cex=2)
	)

# FPCA-based curve reconstruction

# select a curve, use PCs and PC scores 
i=3
s1 = yz_pcafd$scores[i,1]
t_norm = (time_list[[i]] / duration[i]) * mean_dur
op = par(mfrow=c(2,1))
# y
plot(t_norm,y_list[[i]],col = 'red', ylim=c(-Y_max,Y_max),las=1,xlab='lin. norm. time',ylab='y',main=paste('curve',i,', reconstruction'))
coefs = yz_pcafd$meanfd$coefs[,1,1]
lines(fd(coefs,basis),col='blue',lwd=2)
coefs = yz_pcafd$harmonics$coefs[,1,1]
lines(fd(coefs,basis),col='black',lty=2)
coefs = yz_pcafd$meanfd$coefs[,1,1] + s1*yz_pcafd$harmonics$coefs[,1,1]
lines(fd(coefs,basis),col='green',lwd=2)
# z
plot(t_norm,z_list[[i]],col = 'red', ylim=c(-Y_max,Y_max),las=1,xlab='lin. norm. time',ylab='z',main=paste('curve',i,', reconstruction'))
coefs = yz_pcafd$meanfd$coefs[,1,2]
lines(fd(coefs,basis),col='blue',lwd=2)
coefs = yz_pcafd$harmonics$coefs[,1,2]
lines(fd(coefs,basis),col='black',lty=2)
coefs = yz_pcafd$meanfd$coefs[,1,2] + s1*yz_pcafd$harmonics$coefs[,1,2]
lines(fd(coefs,basis),col='green',lwd=2)

par(op)



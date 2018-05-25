# Exercise on Functional PCA 
# Author: Michele Gubian
# Last revision: May 18th 2018


# ex6: peak position and height vary according to independent binary categories
source('header.R')
ex = 'ex6' 
ex_dir = file.path(data_dir,ex)


# load data

n_curves = 50
time_list = list()
y_list = list()
duration = c()


for (i in 1:n_curves) {
	curve = read.table(file.path(ex_dir,paste0('curve',i,'.data')),header=T)
	time_list[[i]] = curve$time
	duration = c(duration,tail(curve$time,1))	
	y_list[[i]] = curve$y
	#mean_y = c(mean_y,mean(	curve$y ))
	#y_list[[i]] = y_list[[i]] - mean_y[i]
}

metadata = read.table(file.path(ex_dir,'metadata.txt'),header=T)
# the two binary categories (A/B and C/D) are totally independent
chisq.test(table(metadata[,-1]))


# some plotting parameters, adjust to your needs
T_min = 0.2; T_max = 0.3
Y_min = 1; Y_max = 2


# exploratory plots 
op = par(ask=T)

for (i in 1:n_curves) {
	plot(time_list[[i]],y_list[[i]],xlim=c(0,T_max), ylim=c(0,Y_max),las=1,xlab='time',ylab='y')
	}
par(op)


# create common basis

# obligatory linear time normalisation
mean_dur = mean(duration) 

norm_range <- c(0,mean_dur)
n_knots = 8 
lambda = 1e-8

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
n_knots = 8
lambda = 1e-8

# smooth all curves
# smooth.basis() does not accept different time samples for different curves.
# Thus we create smooth curves one by one on the same basis, store the spline coefficients and compose an fd object at the end.

y_coefs = matrix(nrow = nbasis, ncol = n_curves)
for (i in 1:n_curves) {
	t_norm = (time_list[[i]] / duration[i]) * mean_dur # linear time normalisation
    	y_coefs[,i] = c(smooth.basis(t_norm,y_list[[i]],fdPar)$fd$coefs)
}
y_fd = fd(coef=y_coefs, basisobj=basis)

# plot according to categories

col_AB = list(A = 'blue', B = 'orange')
col_CD = list(C = 'black', D = 'green')

plot(norm_range,c(0,Y_max),las=1,xlab='lin. norm. time',ylab='y',type='n')
for (i in 1:n_curves) {
	lines(y_fd[i],col = col_AB[[metadata$category_AB[i]]])
}
legend('topleft',legend=names(col_AB),col=unlist(col_AB),lwd=2)

plot(norm_range,c(0,Y_max),las=1,xlab='lin. norm. time',ylab='y',type='n')
for (i in 1:n_curves) {
	lines(y_fd[i],col = col_CD[[metadata$category_CD[i]]])
}
legend('topleft',legend=names(col_CD),col=unlist(col_CD),lwd=2)



# FPCA

lambda_pca    <- lambda
pcafdPar  <- fdPar(basis, 2, lambda_pca) # here Lfdobj = 2 since we don't need derivatives of PC curves.
y_pcafd <- pca.fd(y_fd, nharm=2, pcafdPar) # compute first 2 PCs


# plot PC curves
op <- par(mfrow=c(2,1))
plot.pca.fd(y_pcafd, nx=40)
par(op)


# plot PC scores s1 and s2
xyplot	(y_pcafd$scores[,2] ~  y_pcafd$scores[,1]  , cex=1.5,  
	xlab = list(label=expression(s[1]),cex=2), ylab= list(label=expression(s[2]),cex=2),
	groups= metadata$category_AB, pch = c('A','B'), col = unlist(col_AB)
	#groups= metadata$category_CD, pch = c('C','D'), col = unlist(col_CD)
	)
# plot PC score s1 or s2 
bwplot( y_pcafd$scores[,2] ~ metadata$category_AB, xlab ='',ylab=list(label=expression(s[1]),cex=2), scales=list(cex=1.5) )


# FPCA-based curve reconstruction

# select a curve, use PCs and PC scores 
i = 1
s1 = y_pcafd$scores[i,1]
s2 = y_pcafd$scores[i,2]


plot(norm_range,c(-1*Y_max,1*Y_max),type='n',las=1,xlab='lin. norm. time',ylab='y',main=paste('curve',i,', reconstructionn'))

lines(y_fd[i],col = 'black',lwd=2)
lines(y_pcafd$harmonics[1],col='red',lwd=2)
lines(y_pcafd$harmonics[2],col='orange',lwd=2)

lines(y_pcafd$meanfd,col='blue',lwd=2)
lines(y_pcafd$meanfd+s1*y_pcafd$harmonics[1],col='green',lwd=2)
lines(y_pcafd$meanfd+s1*y_pcafd$harmonics[1]+s2*y_pcafd$harmonics[2],col='violet',lwd=2)

# average curve per category
# try s1, s2, AB, CD in all combinations
plot(norm_range,c(0,1.5*Y_max),type='n',las=1,xlab='lin. norm. time',ylab='y')
for (AB in c('A','B')) {
	mean_s2 =  mean( y_pcafd$scores[metadata$category_AB == AB,2] )
	lines(y_pcafd$meanfd + mean_s2 * y_pcafd$harmonics[2], lwd=2, col=col_AB[[AB]])
}
legend('topleft',legend=names(col_AB),col=unlist(col_AB),lwd=2)



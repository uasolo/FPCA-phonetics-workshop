# Exercise on Functional PCA 
# Author: Michele Gubian
# Last revision: May 18th 2018

# ex10: landmark reg, symmetric hill, one mid landmark randomly before or after peak, related to category A/B
# ex11: landmark reg, symmetric hill, one mid landmark at same position, before/after peak and duration
# related to category A/B

source('header.R')
ex = 'ex10' # change according to ex number
ex_dir = file.path(data_dir,ex)

# load data

n_curves = 50
time_list = list()
y_list = list()
duration = c()

for (i in 1:n_curves) {
	curve = read.table(file.path(ex_dir,paste0('curve',i,'.data')),header=T)
	time_list[[i]] = curve$time
	y_list[[i]] = curve$y
	duration = c(duration,tail(curve$time,1))
}

metadata = read.table(file.path(ex_dir,'metadata.txt'),header=T)
land = as.matrix(read.table(file.path(ex_dir,'landmarks.txt'),header=F))

# some plotting parameters, adjust to your needs
T_min = 0.2; T_max = 0.3
Y_min = 1; Y_max = 2


# exploratory plots 
op = par(ask=T)
for (i in 1:n_curves) {
	plot(time_list[[i]],y_list[[i]],xlim=c(0,T_max), ylim=c(0,Y_max),las=1,xlab='time',ylab='y', main = metadata$category_AB[i])
	abline(v=land[i,],lty=3,col='black')
}
par(op)


# ex10: check that landmark tends to be before/after half time interval for A/B curves (optional)
library("data.table")
land.table <- as.data.table(land)
land.table[, c("d_ratio", "AB") := list(V2/V3, metadata$category_AB)]
library("ggplot2")
ggplot(land.table, aes(x=AB, y=d_ratio)) + geom_boxplot()


# landmark registration

n_landmarks=dim(land)[2]
reg = landmarkreg.nocurve(land,nhknots = 8, hlambda=1e-8, wlambda =1e-8)

# plot h(t)'s
col_AB = list(A = 'blue', B = 'orange')
set.seed(1234)
plot(c(0,mean(duration)), c(0, max(duration)), type = 'n', xlab = "reg. time", ylab = "time", las = 1)
abline(v=reg$land[2], lty=2)

for (i in sample(n_curves, 20)) {
  points(reg$land[2], land[i,2], pch = 20, col = 'red')
  lines(reg$warpfd[i], col = col_AB[[metadata$category_AB[i]]])
}
legend('topleft',legend=names(col_AB),col=unlist(col_AB),lwd=2)
# plot log rate = - log (dh(t)/dt)
set.seed(1234)
plot(c(0,mean(duration)), c(-1,1), type = 'n', xlab = "reg. time", ylab = "log rate", las = 1)
abline(v=reg$land[2], lty=2)
for (i in sample(1:n_curves, 20)) {
  lines(reg$logvelfd[i], col = col_AB[[metadata$category_AB[i]]])
}
legend('topleft',legend=names(col_AB),col=unlist(col_AB),lwd=2)

# comparere h(t) with log rate
i = 2 # pick one 
op <- par(mfrow=c(2,1))
plot(c(0,mean(duration)), c(0, max(duration)), type = 'n', xlab = "reg. time", ylab = "time", las = 1, main = "h(t)")
abline(v=reg$land[2], lty=2)
lines(c(0,mean(duration)),c(0,mean(duration)), lty = 2, col = "red")
lines(reg$warpfd[i], col = col_AB[[metadata$category_AB[i]]])
plot(c(0,mean(duration)), c(-1,1), type = 'n', xlab = "reg. time", ylab = "log rate", las = 1, main = "log rate = - log (dh(t)/dt)")
abline(v=reg$land[2], lty=2)
lines(c(0,mean(duration)),c(0,0), lty = 2, col = "red")
lines(reg$logvelfd[i], col = col_AB[[metadata$category_AB[i]]])
par(op)



# smoothing

# create common basis

# common duration based on landmark registration
mean_dur = reg$land[n_landmarks]

norm_range <- c(0,mean_dur)
n_knots = 8 # try many
lambda = 1e-8

#knots <- seq(0,mean_dur,length.out = n_knots)
Lfdobj <- 3 # 2 + order of derivative expected to be used. We might need velocity, thus order = 1
norder <- 2 + Lfdobj  # a fixed relation about B-splines
nbasis <- n_knots + norder - 2 # a fixed relation about B-splines
basis <- create.bspline.basis(norm_range, nbasis, norder) #, knots)
fdPar <- fdPar(basis, Lfdobj, lambda)


# parameters set to:
n_knots = 8
lambda = 1e-8

# smooth all curves and apply time warping computed above
y_coefs = matrix(nrow = nbasis, ncol = n_curves)
y_nonreg_fd = list() 
for (i in 1:n_curves) {
	t = time_list[[i]]
	y = y_list[[i]]
	range_i = range(c( reg$hfunmat[,i], t) ) # prevent rounding errors
	basis_i = create.bspline.basis(range_i,nbasis,norder)
	lambda_i = lambda
	fdPar_i <- fdPar(basis_i,Lfdobj,lambda=lambda_i)
	y_nonreg_fd[[i]] = smooth.basis(t,y,fdPar_i)$fd
	# see notes/landmark_reg.pdf
	y_coefs[,i] = c(smooth.basis(reg$x, eval.fd(reg$hfunmat[,i], y_nonreg_fd[[i]]),fdPar)$fd$coefs)
}

y_fd = fd(coef=y_coefs, basisobj=basis)
# curves are linearly time normalized, their duration is mean_dur
plot(y_fd)


### code to produce plots in notes/landmark_reg.pdf
# ex 10
i = 6
t_ind <- c(10,20,30,40,50)
n_comb <- length(t_ind)
t_comb <- reg$x[t_ind]
h_comb <- reg$hfunmat[t_ind, i]
cex_points <- 1.3


png(file.path(plots_dir, "h.png"))
plot(c(0,mean(duration)), c(0, max(duration)), type = 'n', xlab = "reg. time", ylab = "time", las = 1, main = "")
lines(reg$warpfd[i], col = 'black', lwd = 2)
segments(0, h_comb, t_comb, h_comb, lty = 3, col = 'black')
segments(t_comb, 0, t_comb, h_comb, lty = 3, col = 'black')
points(rep(0, n_comb), h_comb, pch = 17, col = 'red', cex = cex_points)
points(t_comb, rep(0, n_comb), pch = 15, col = 'blue', cex = cex_points)
dev.off()

y_samp <- eval.fd(h_comb, y_nonreg_fd[[i]])

png(file.path(plots_dir, "y_orig.png"))
plot(y_nonreg_fd[[i]], xlab = "time", ylab = "y", las = 1, lwd = 2)
segments(0, y_samp, h_comb, y_samp, lty = 3, col = 'black')
segments(h_comb, 0, h_comb, y_samp, lty = 3, col = 'black')
points(rep(0, n_comb), y_samp, pch = 19, col = 'green', cex = cex_points)
points(h_comb, rep(0, n_comb), pch = 17, col = 'red', cex = cex_points)
dev.off()


png(file.path(plots_dir, "y_reg.png"))
plot(y_fd[i], xlab = "reg. time", ylab = "y", las = 1, lwd = 2)
segments(0, y_samp, t_comb, y_samp, lty = 3, col = 'black')
segments(t_comb, 0, t_comb, y_samp, lty = 3, col = 'black')
points(rep(0, n_comb), y_samp, pch = 19, col = 'green', cex = cex_points)
points(t_comb, rep(0, n_comb), pch = 15, col = 'blue', cex = cex_points)
dev.off()

### end of plots for notes/landmark_reg.pdf


# plot according to categories

col_AB = list(A = 'blue', B = 'orange')

plot(norm_range,c(0,Y_max),las=1,xlab='lin. norm. time',ylab='y',type='n')
abline(v=reg$land[2], lty=2)
for (i in 1:n_curves) {
	lines(y_fd[i],col = col_AB[[metadata$category_AB[i]]])
}
legend('topleft',legend=names(col_AB),col=unlist(col_AB),lwd=2)


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
	)
# plot PC score s1
bwplot( y_pcafd$scores[,1] ~ metadata$category_AB, xlab ='',ylab=list(label=expression(s[1]),cex=2), scales=list(cex=1.5) )



# FPCA-based curve reconstruction

# select a curve, use PCs and PC scores 
i = 2
s1 = y_pcafd$scores[i,1]
s2 = y_pcafd$scores[i,2]


plot(norm_range,c(-1*Y_max,1*Y_max),type='n',las=1,xlab='lin. norm. time',ylab='y',main=paste('curve',i,', reconstructionn'))

lines(y_fd[i],col = 'black',lwd=2)
lines(y_pcafd$harmonics[1],col='red',lwd=2)

lines(y_pcafd$meanfd,col='blue',lwd=2)
lines(y_pcafd$meanfd+s1*y_pcafd$harmonics[1],col='green',lwd=2)
lines(y_pcafd$meanfd+s1*y_pcafd$harmonics[1]+s2*y_pcafd$harmonics[2],col='violet',lwd=2)

# average curve per category

plot(norm_range,c(0,1.5*Y_max),type='n',las=1,xlab='lin. norm. time',ylab='y')
for (AB in c('A','B')) {
	mean_s1 =  mean( y_pcafd$scores[metadata$category_AB == AB,1] )
	lines(y_pcafd$meanfd + mean_s1 * y_pcafd$harmonics[1], lwd=2, col=col_AB[[AB]])
}
legend('topleft',legend=names(col_AB),col=unlist(col_AB),lwd=2)


# class analysis



d_temp = t(diff(t(land)))
colnames(d_temp) = c('d1','d2')
metadata = cbind(metadata, d_temp)
metadata$s1 = y_pcafd$scores[,1]

splom(~ cbind(d1 , d2 , s1), groups= category_AB, data=metadata)
metadata$ratio_d1d2 = metadata$d1/metadata$d2
metadata$sum_d1d2 = metadata$d1+metadata$d2
splom(~ cbind(ratio_d1d2, sum_d1d2, s1), groups= category_AB, data=metadata)



######## joint FPCA on registered contours and their relative speech rate (duration)

t_samp = reg$x

# convert logvelfd to basis
logvelfd = smooth.basis(t_samp, eval.fd(t_samp,reg$logvelfd),fdPar)$fd
# plot log velocity according to categories
col_AB = list(A = 'blue', B = 'orange')
plot(norm_range,c(-log(2),log(2)),las=1,xlab='lin. norm. time',ylab='rel. velocity',type='n')
for (i in 1:n_curves) {
	lines(logvelfd[i],col = col_AB[[metadata$category_AB[i]]])
}
abline(v=reg$land)
legend('topleft',legend=names(col_AB),col=unlist(col_AB),lwd=2)



# 2D fd object
Y_coefs = array(dim = c(dim(y_fd$coefs),2))
Y_coefs[,,1] = y_fd$coefs
Y_coefs[,,2] = logvelfd$coefs 
# make the two different units comparable
w2 = mean(eval.fd(t_samp,sd.fd(y_fd)) / eval.fd(t_samp,sd.fd(logvelfd)))
#w2 =  mean(eval.fd(t_samp,sd.fd(y_fd))) / mean( eval.fd(t_samp,sd.fd(logvelfd)))
Y_coefs[,,2] = Y_coefs[,,2] * w2

Y_fd = fd(coef=Y_coefs, basisobj=basis, fdnames=list('norm. time (s)' ,as.character(metadata$curve), c('y', 'rel. velocity') ))


# FPCA

lambda_pca    <- lambda
pcafdPar  <- fdPar(basis, 2, lambda_pca) # here Lfdobj = 2 since we don't need derivatives of PC curves.
Y_pcafd <- pca.fd(Y_fd, nharm=2, pcafdPar) # compute first 2 PCs

# plot PC curves
op <- par(mfrow=c(2,2))
plot.pca.fd.corr(Y_pcafd,xlab = 'norm. time (s)',ylab=c('y','rel. log vel.'),land = reg$land ,nx=20)
par(op)

# plot PC scores s1 and s2
xyplot	(s2 ~ s1, data = data.frame(s1 = Y_pcafd$scores[,1], s2 = Y_pcafd$scores[,2])  , cex=1.5, 
	xlab = list(label=expression(s[1]),cex=2), ylab= list(label=expression(s[2]),cex=2),
	groups= metadata$category_AB, pch = c('A','B'), col = unlist(col_AB)
	)

# ex10:
# PC1 captures variation of landmark position w.r.t. peak
# PC2 captures global duration variation independent of landmark position

# plot representative curves for each category

# plot y

plot(basis$rangeval,c(0,Y_max),type='n',las = 1,xlab = 'normalised time (s)', ylab = 'y',cex.lab = 1.5,cex.axis=1.5)
abline(v=reg$land,lty=2,col='black')
for (AB in c('A','B')) {
	mean_s1 =  mean( Y_pcafd$scores[metadata$category_AB == AB,1] )
    	coefs = Y_pcafd$meanfd$coefs[,1,1] + mean_s1 * Y_pcafd$harmonics$coefs[,1,1] 
	lines(fd(coefs,basis), lwd=2, col=col_AB[[AB]])
    }
legend('topleft',legend=names(col_AB),col=unlist(col_AB),lwd=2,bg='white')


# duration visualised as barplot
# dur contains the durations of each segment in each computed curve
# nrow = number of segments, ncol = number of categories
dur = matrix(nrow = n_landmarks-1, ncol = 2)
colnames(dur) = c('A','B')
for (AB in c('A','B') ) {
	mean_s1 =  mean( Y_pcafd$scores[metadata$category_AB == AB,1] )
    	coefs_logvel = Y_pcafd$meanfd$coefs[,1,2] + mean_s1 * Y_pcafd$harmonics$coefs[,1,2] 
	rel_dur_fd = smooth.basis(t_samp,exp(-(1/w2)*eval.fd(t_samp,fd(coefs_logvel,basis))),fdPar)$fd
    	for (l in 1:(n_landmarks-1)) {
        dur[l, AB] = defint.fd(rel_dur_fd,c(reg$land[l],reg$land[l+1]))
    }
}

# barplot
bp = barplot(dur[,rev(seq_len(ncol(dur)))],names.arg=c('B','A'),horiz =TRUE,las=1,cex.axis=1.5,cex.names=1.5,cex.lab=1.5,col= c('gray85','gray75'),xlab = 'time (s)',ylab = 'category_AB') # the categories are presented to barplot in reverse order because otherwise they would be ordered from bottom to top

# alternative representation: reverse landmark register the representative curves
land_AB = apply(dur,2,cumsum) # cumulative durations of segments
land_AB = rbind(c(0,0),land_AB)
land_AB = as.data.frame(land_AB)
colnames(land_AB) = c('A','B')


plot(c(0,T_max),c(0,Y_max),type='n',las = 1,xlab = 'time (s)', ylab = 'y',cex.lab = 1.5,cex.axis=1.5)
for (AB in c('A','B')) {
	mean_s1 =  mean( Y_pcafd$scores[metadata$category_AB == AB,1] )
	coefs_y = Y_pcafd$meanfd$coefs[,1,1] + mean_s1 * Y_pcafd$harmonics$coefs[,1,1]
	y_reg_fd = fd(coefs_y,basis) # curve in registered time
	
	
	reg_y = landmarkreg.nocurve(matrix(reg$land,ncol=3),land_AB[[AB]],nhknots = 5, hlambda=1e-8, wlambda =1e-8)
	h = reg_y$hfunmat[,1]
	h[1] = 0 ; h[length(h)] = reg$land[3] # prevent rounding errors
	n_knots = 8 
	lambda = 1e-8
	Lfdobj <- 2
	norder <- 2 + Lfdobj
	nbasis <- n_knots + norder - 2
	basis_y <- create.bspline.basis(range(land_AB[[AB]]), nbasis, norder)
	fdPar_y <- fdPar(basis_y, Lfdobj, lambda)
	
	y_unreg_fd = smooth.basis(reg_y$x, eval.fd(h,y_reg_fd),fdPar_y)$fd
	
	lines(y_unreg_fd, lwd=2, col=col_AB[[AB]])
	points(land_AB[[AB]],eval.fd(land_AB[[AB]],y_unreg_fd),pch=18, cex=2,col=col_AB[[AB]])
}
legend('topleft',legend=names(col_AB),col=unlist(col_AB),lwd=2,bg='white')



plot.pca.fd.corr <-  
function (x, nx = 128,height=240, pointplot = TRUE, harm = 0, expand = 0, pcweight = NULL,
    cycle = FALSE, land = NA , xlab = NULL, ylab = NULL, landlab = NULL, png = FALSE, plots_dir = NULL, basename = NULL,ylog = NULL,...)
{
	# removing the long titles on top of the figures
	# enlarging '+' and '-' signs 
	# pcweight is a multiplicative factor in order to represent PC functions in their original scale.
	# Useful if multidim signals were used and they were previously scaled to e.g. similar st dev.
	# ylog = 'n' no log scale. ylog = 'y' used for logvel, meaning that y values are already logs and are transformed back to exp, and a log scale in y is provided.
	# in case of multidim plots, ylog is a vector, e.g. = c('n','y') to have log scale only on the second dim.
	
	log_fun = list (n = identity, y = exp)
	log_str = list (n = '', y = 'y')
	cex = 1.3
	

    pcafd <- x

    if (!(inherits(pcafd, "pca.fd")))
        stop("Argument PCAFD is not a pca.fd object.")
	
    harmfd <- pcafd[[1]]
    basisfd <- harmfd$basis
    rangex <- basisfd$rangeval
    x <- seq(rangex[1], rangex[2], length = nx)
    fdmat <- eval.fd(x, harmfd)
    meanmat <- eval.fd(x, pcafd$meanfd)
    dimfd <- dim(fdmat)
    nharm <- dimfd[2]
    harm <- as.vector(harm)
	if (is.null(xlab))
	xlab = pcafd$meanfd$fdnames[[1]]
	if (is.null(ylab))
	ylab = pcafd$harmonics$fdnames[[3]]
    if (harm[1] == 0)
        harm <- (1:nharm)
    if (length(dimfd) == 2) { # one dim signals
         if (is.null(ylog)) {   ylog = 'n' }
	if (is.null(pcweight)) 
		pcweight = 1
        for (iharm in harm) {
            if (expand == 0)
             #   fac <- sqrt(pcafd$values[iharm])
		fac = sd(pcafd$scores[,iharm])
            else fac <- expand
            vecharm <- fdmat[, iharm]
            pcmat <- cbind(meanmat + fac * vecharm, meanmat -
                fac * vecharm)
            if (pointplot)
                plottype <- "p"
            else plottype <- "l"
            percentvar <- round(100 * pcafd$varprop[iharm], 1)
		if (!is.null(plots_dir))
		if (png) {png(paste(plots_dir,basename,iharm,'.png',sep=''),height=height)}
		   
            plot(x, sapply(pcweight * meanmat,log_fun[[ylog]]), type = "l", 
            ylim = sapply(c(min(pcmat), max(pcmat)),log_fun[[ylog]]), # ylab = paste("PC", iharm),
		xlab = xlab, ylab = ylab,
                main = paste("PC", iharm, " (",
                  percentvar, "% var)"),cex=cex ,las=1,cex.main=1,log = log_str[[ylog]],...)
	# plot landmarks	
		if (!is.null(land)) {
			for (l in land) {
				abline(v = l, lty=2)
			}
		}
		if (!is.null(landlab)) {
			at = c()
			for (i in 1:(length(land)-1)) {
				at = c(at, mean(land[i:(i+1)]))
			}
			axis(3,tick=F,at=at, labels=landlab,cex.axis=0.8)
		}
	     
            if (pointplot) {
		# hack here + and - positions if necessary
                points(x, sapply(pcweight * pcmat[, 1],log_fun[[ylog]]), pch = "+",cex=cex) # original: +
                points(x, sapply(pcweight * pcmat[, 2],log_fun[[ylog]]), pch = "-",cex=cex) # original: -
            }
            else {
                lines(x, sapply(pcweight[jvar] * pcmat[, 1],log_fun[[ylog]]), lty = 2)
                lines(x, sapply(pcweight[jvar] * pcmat[, 2],log_fun[[ylog]]), lty = 3)
            }
		#if (!is.null(plots_dir))
		if (png) {dev.off()}
        }
    }
    else {
        if (cycle && dimfd[3] == 2) { # multi dim, phase plot
            meanmat <- drop(meanmat)
            for (iharm in harm) {
                if (expand == 0)
                  fac <- 2 * sqrt(pcafd$values[iharm])
                else fac <- expand
                matharm <- fdmat[, iharm, ]
                mat1 <- meanmat + fac * matharm
                mat2 <- meanmat - fac * matharm
                if (pointplot)
                  plottype <- "p"
                else plottype <- "l"
                percentvar <- round(100 * pcafd$varprop[iharm],
                  1)
                plot(meanmat[, 1], sapply(meanmat[, 2],log_fun[[ylog]]), type = plottype,
                  xlim = c(min(c(mat1[, 1], mat2[, 1])), max(c(mat1[,
                    1], mat2[, 1]))), ylim = sapply(c(min(c(mat1[, 2],
                    mat2[, 2])), max(c(mat1[, 2], mat2[, 2]))),log_fun[[ylog]]),
                  main = paste("PCA function", iharm, "(Percentage of variability",
                    percentvar, ")"),log=log_str[[ylog]], ...)
                if (pointplot) {
                  points(mat1[, 1], sapply(mat1[, 2],log_fun[[ylog]]), pch = "+",cex=cex)
                  points(mat2[, 1], sapply(mat2[, 2],log_fun[[ylog]]), pch = "-",cex=cex)
                }
                else {
                  lines(mat1[, 1], sapply(mat1[, 2],log_fun[[ylog]]), lty = 2)
                  lines(mat2[, 1], sapply(mat2[, 2],log_fun[[ylog]]), lty = 3)
                }
            }
        }
        else {  # multi dim, regular plot
	        for (iharm in harm) {
                if (expand == 0)
                 # fac <- sqrt(pcafd$values[iharm])
		fac = sd(pcafd$scores[,iharm])
                else fac <- expand
                meanmat <- drop(meanmat)
                matharm <- fdmat[, iharm, ]
                nvar <- dim(matharm)[2]
                if (is.null(ylog)) {ylog = rep('n',nvar)} 
		if	(is.null(pcweight))
		pcweight = rep(1,nvar)
                for (jvar in 1:nvar) {
                  pcmat <- cbind(meanmat[, jvar] + fac * matharm[,
                    jvar], meanmat[, jvar] - fac * matharm[,
                    jvar])
		ylim = sapply(pcweight[jvar] * range( as.vector( pcmat )),log_fun[[ylog[jvar]]])
                 # if (pointplot)
                 #   plottype <- "p"
                #  else
		 plottype <- "l"
                  percentvar <- round(100 * pcafd$varprop[iharm],
                    1)
		if (!is.null(plots_dir))
		if (png) {png(paste(plots_dir,basename,iharm,jvar,'.png',sep=''),height=height)}
                  plot(x, sapply(pcweight[jvar] * meanmat[, jvar],log_fun[[ylog[jvar]]]), type = plottype,  # ylab = paste("Harmonic", iharm),
		main = paste("PC", iharm, " (", percentvar, "% var)"),
		xlab = xlab, ylab = ylab[jvar],
		ylim = ylim,cex = cex, las=1, cex.main=1,log=log_str[[ylog[jvar]]],...)
			# plot landmarks	
		if (!is.null(land)) {
			for (l in land) {
				abline(v = l, lty=2)
			}
		}
		if (!is.null(landlab)) {
			at = c()
			for (i in 1:(length(land)-1)) {
				at = c(at, mean(land[i:(i+1)]))
			}
			axis(3,tick=F,at=at, labels=landlab,cex.axis=0.8)
		}
	

                  if (pointplot) {
                    points(x, sapply(pcweight[jvar] * pcmat[, 1],log_fun[[ylog[jvar]]]), pch = "+",cex=cex)
                    points(x, sapply(pcweight[jvar] * pcmat[, 2],log_fun[[ylog[jvar]]]), pch = "-",cex=cex)
                  }
                  else {
                    lines(x,sapply(pcweight[jvar] * pcmat[, 1],log_fun[[ylog[jvar]]]), lty = 2)
                    lines(x,sapply(pcweight[jvar] * pcmat[, 2],log_fun[[ylog[jvar]]]), lty = 3)
                  }
		#if (!is.null(plots_dir))
		if (png) {dev.off()}
                }
            }
        }
    }
    invisible(NULL)
}


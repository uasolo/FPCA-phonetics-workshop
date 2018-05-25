pca.fd <- function(fdobj, nharm = 2, harmfdPar=fdPar(fdobj),
                   centerfns = TRUE)
{
#  Carry out a functional PCA with regularization
#  Arguments:
#  FDOBJ      ... Functional data object
#  NHARM     ... Number of principal components or harmonics to be kept
#  HARMFDPAR ... Functional parameter object for the harmonics
#  CENTERFNS ... If TRUE, the mean function is first subtracted from each function
#
#  Returns:  An object PCAFD of class "pca.fd" with these named entries:
#  harmonics  ... A functional data object for the harmonics or eigenfunctions
#  values     ... The complete set of eigenvalues
#  scores     ... A matrix of scores on the principal components or harmonics
#  varprop    ... A vector giving the proportion of variance explained
#                 by each eigenfunction
#  meanfd     ... A functional data object giving the mean function
#

#  Last modified:  3 January 2008 by Jim Ramsay
# Previously modified 2007 April 26 by Spencer Graves

  #  Check FDOBJ

  if (!(inherits(fdobj, "fd"))) stop(
    "Argument FD  not a functional data object.")

  #  compute mean function and center if required

  meanfd <- mean.fd(fdobj)
  if (centerfns) fdobj <- center.fd(fdobj)

  #  get coefficient matrix and its dimensions

  coef  <- fdobj$coefs
  coefd <- dim(coef)
  ndim  <- length(coefd)
  nrep  <- coefd[2]
  coefnames <- dimnames(coef)
  if (nrep < 2) stop("PCA not possible without replications.")

  basisobj <- fdobj$basis
  nbasis   <- basisobj$nbasis
  type     <- basisobj$type

  #  set up HARMBASIS
  #  currently this is required to be BASISOBJ

  harmbasis <- basisobj

  #  set up LFDOBJ and LAMBDA

  Lfdobj <- harmfdPar$Lfd
  lambda <- harmfdPar$lambda

  #  compute CTEMP whose cross product is needed

  if (ndim == 3) {
    nvar <- coefd[3]
    ctemp <- matrix(0, nvar * nbasis, nrep)
    for(j in 1:nvar) {
      index <- 1:nbasis + (j - 1) * nbasis
      ctemp[index,  ] <- coef[,  , j]
    }
  } else {
    nvar  <- 1
    ctemp <- coef
  }

  #  set up cross product and penalty matrices

  Cmat <- crossprod(t(ctemp))/nrep
  Jmat <- eval.penalty(basisobj, 0)
  if(lambda > 0) {
    Kmat <- eval.penalty(basisobj, Lfdobj)
    Wmat <- Jmat + lambda * Kmat
  } else {
    Wmat <- Jmat
  }
  Wmat <- (Wmat + t(Wmat))/2

  #  compute the Choleski factor of Wmat

  Lmat    <- chol(Wmat)
  Lmatinv <- solve(Lmat)

  #  set up matrix for eigenanalysis

  if(nvar == 1) {
    if(lambda > 0) {
            Cmat <- t(Lmatinv) %*% Jmat %*% Cmat %*% Jmat %*% Lmatinv
    } else {
            Cmat <- Lmat %*% Cmat %*% t(Lmat)
    }
  } else {
    for (i in 1:nvar) {
      indexi <- 1:nbasis + (i - 1) * nbasis
      for (j in 1:nvar) {
        indexj <- 1:nbasis + (j - 1) * nbasis
        if (lambda > 0) {
          Cmat[indexi, indexj] <- t(Lmatinv) %*% Jmat %*%
          Cmat[indexi, indexj] %*% Jmat %*% Lmatinv
        } else {
          Cmat[indexi, indexj] <- Lmat %*% Cmat[indexi,indexj] %*% t(Lmat)
        }
      }
    }
  }

  #  eigenalysis

  Cmat    <- (Cmat + t(Cmat))/2
  result  <- eigen(Cmat)
  eigvalc <- result$values
  eigvecc <- as.matrix(result$vectors[, 1:nharm])
  sumvecc <- apply(eigvecc, 2, sum)
  eigvecc[,sumvecc < 0] <-  - eigvecc[, sumvecc < 0]

  varprop <- eigvalc[1:nharm]/sum(eigvalc)

  if (nvar == 1) {
    harmcoef <- Lmatinv %*% eigvecc
    harmscr  <- t(ctemp) %*% t(Lmat) %*% eigvecc
  } else {
    harmcoef <- array(0, c(nbasis, nharm, nvar))
    harmscr  <- matrix(0, nrep, nharm)
    for (j in 1:nvar) {
      index <- 1:nbasis + (j - 1) * nbasis
      temp <- eigvecc[index,  ]
      harmcoef[,  , j] <- Lmatinv %*% temp
      harmscr <- harmscr + t(ctemp[index,  ]) %*% t(Lmat) %*% temp
    }
  }
  harmnames <- rep("", nharm)
  for(i in 1:nharm)
    harmnames[i] <- paste("PC", i, sep = "")
  if(length(coefd) == 2)
    harmnames <- list(coefnames[[1]], harmnames,"values")
  if(length(coefd) == 3)
    harmnames <- list(coefnames[[1]], harmnames, coefnames[[3]])
  harmfd   <- fd(harmcoef, basisobj, harmnames)

  pcafd        <- list(harmfd, eigvalc, harmscr, varprop, meanfd)
  class(pcafd) <- "pca.fd"
  names(pcafd) <- c("harmonics", "values", "scores", "varprop", "meanfd")

  return(pcafd)
}

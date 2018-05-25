landmarkreg.nocurve <- function(ximarks, x0marks=xmeanmarks,
                        WfdPar=NULL, nhknots = 30, hlambda=1e-7, wlambda =1e-7 )
{

# monwrd=FALSE , let us suppress the monotone one, it behaves strangely.
# This is a modified version of fda::landmarkreg
# by Michele Gubian, 1 October 2010
# It allows to perform a landmark registration of time axes without 
# any associated curve. It is used as an intermediate step to build 
# relative velocity curves, e.g. once a forced alignment is performed 
# on realizations of the same lexical material.

# Parameters: the fdobj parameter is then missing from the origial fda::landmarkreg paramenter set.


#  Arguments:

#  XIMARKS ... N by (NL + 2) matrix of times of landmarks (NL is the number of interior landmarks) 
#		In contrast with the landmark reg command, also the first and last lanmarks are provided,
#		because in general global durations can differ
#  XOMARKS ... vector of length (NL + 2 ) output landmarks, i.e. where all the input landmarks should be aligned
#                 
#  WFDPAR  ... a functional parameter object defining a warping function
#  MONWRD  ... If TRUE, warping functions are estimated by monotone smoothing,
#                 otherwise by regular smoothing.  The latter is faster, but
#                 not guaranteed to produce a strictly monotone warping
#                 function.  If MONWRD is 0 and an error message results
#                 indicating nonmonotonicity, rerun with MONWRD = 1.
#                 Default:  TRUE
#  YLAMBDA ... smoothing parameter to be used in computing the registered
#                 functions.  For high dimensional bases, local wiggles may be 
#                 found in the registered functions or its derivatives that are
#                 not seen in the unregistered functions.  In this event, this
#                 parameter should be increased.             
#  Returns:
#  WARPFD  ... a functional data object for the warping functions
#  WFD     ... a functional data object for the W functions defining the 
#              warping functions
# LOGSPEED ...  a functional data object for the log of the relative speed referred to the XOMARKS landmarks 

 #  Last modified 15 May 2009 by Jim Ramsay

  #  check FDOBJ

#  if (!(inherits(fdobj,  "fd"))) stop(
#		"Argument fdobj  not a functional data object.")

  #  extract information from curve functional data object and its basis

#  coef   <- fdobj$coefs
#  coefd  <- dim(coef)
#  ndim   <- length(coefd)
#  ncurve <- coefd[2]
#  if (ndim > 2) {
#      nvar <- coefd[3]
#  } else {
#      nvar <- 1
#  }


#  basisobj <- fdobj$basis
#  type     <- basisobj$type
#  nbasis   <- basisobj$nbasis
#  rangeval <- basisobj$rangeval
#  fdParobj <- fdPar(basisobj, 2, ylambda)

  #  check landmarks

#  if (is.vector(ximarks)) ximarks = as.matrix(ximarks)
#  ximarksd <- dim(ximarks)
#  if (ximarksd[1] != ncurve) stop(
#     "Number of rows of XIMARKS is incorrect.")
#    if (any(ximarks <= rangeval[1]) || any(ximarks >= rangeval[2])) stop(
#     "Some landmark values are not within the range.")

  nlandm <- dim(ximarks)[2]
  ncurve = dim(ximarks)[1]
# let all landmark sequences start from 0
ximarks =  t(scale(t(ximarks), center = ximarks[,1] , scale = FALSE))
landerr = matrix(nrow= ncurve, ncol=nlandm) # contains the diff between the average landmarks and the aligned ones

  xmeanmarks <- apply(ximarks,2,mean)
# mean of landmarks incrementally
#  xmeanmarks <- diffinv(apply(diff(t(land)),1,mean))

  if (length(x0marks) != nlandm) stop(
     "Number of target landmarks not equal to number of curve landmarks.")
rangex = range(x0marks)
 rangeval = rangex
  #  set up default WfdPar
  
#  if (is.null(WfdPar)) {
#    basisobj  <- fdobj$basis
#    rangex    <- basisobj$rangeval
    worder = 4 
    wnbasis   <- length(x0marks) + worder - 2
 #    wnbasis   <- 10* length(x0marks) # test
    wbasis    <- create.bspline.basis(rangex, wnbasis, worder, x0marks) # knots on norm landmarks, de Boor theorem
    WfdPar <- fdPar(wbasis,lambda=wlambda)
#  }
   nbasis = wnbasis
  #  check WFDPAR
   
  #WfdPar <- fdParcheck(WfdPar)
		
  #  set up WFD0 and WBASIS

  Wfd0   <- WfdPar$fd
  wLfd   <- WfdPar$Lfd
  wbasis <- Wfd0$basis

  #  set up WLAMBDA


  wlambda <- WfdPar$lambda

  #  check landmark target values

#  wrange <- wbasis$rangeval
#  if (any(rangeval != wrange)) stop(
#		"Ranges for FD and WFDPAR do not match.")

  #  set up analysis

  n   <- min(c(101,10*nbasis))
  x   <- seq(rangeval[1],rangeval[2],length=n)

#  y       <- eval.fd(x, fdobj)
#  yregmat <- y
  hfunmat <- matrix(0,n,ncurve)
  wlambda  <- max(wlambda,1e-10)

  xval    <- x0marks
  nwbasis <- wbasis$nbasis
  Wcoef   <- matrix(0,nwbasis,ncurve)
  nval    <- length(xval)
  wval    <- rep(1,nval)

  #  --------------------------------------------------------------------
  #                  Iterate through curves to register
  #  --------------------------------------------------------------------

  cat("Progress:  Each dot is a curve\n")

  for (icurve in 1:ncurve) {
    cat(".")
    #  set up landmark times for this curve
    yval   <- ximarks[icurve,]
    #  smooth relation between this curve"s values and target"s values
 #   if (monwrd) {
       #  use monotone smoother
       Wfd       <- smooth.morph(xval, yval, WfdPar)$Wfdobj
       h         <- monfn(x, Wfd)
 #      b         <- (rangeval[2]-rangeval[1])/(h[n]-h[1])
       b         <- (diff(range(yval)))/(h[n]-h[1])
 #      a         <- rangeval[1] - b*h[1]
	a=0
       h         <- a + b*h
	h[1] = 0
	h[n] = max(yval)
 #      h[c(1,n)] <- rangeval
 #      wcoefi    <- Wfd$coef
 #      Wcoef[,icurve] <- wcoefi
 #   } else {
       #  use unconstrained smoother
 #      warpfd <- smooth.basis(xval, yval, WfdPar, wval)$fd
       #  set up warping function by evaluating at sampling values
 #      h         <- as.vector(eval.fd(x, warpfd))
#       b         <- (rangeval[2]-rangeval[1])/(h[n]-h[1])
#       a         <- rangeval[1] - b*h[1]
#       h         <- a + b*h
#       h[c(1,n)] <- rangeval
       #  check for monotonicity
       deltah <- diff(h)
       if (any(deltah <= 0)) stop(
           paste("Non-increasing warping function estimated for curve",icurve,
                 " Try setting MONWRD to TRUE."))
  #     wcoefi    <- warpfd$coef
  #     Wcoef[,icurve] <- wcoefi
  #  }
    hfunmat[,icurve] <- h

    #  compute h-inverse  in order to register curves

 #   if (monwrd) {
 #      wcoef        <- Wfd$coefs
 #      Wfdinv       <- fd(-wcoef,wbasis)
 #      WfdParinv    <- fdPar(Wfdinv, wLfd, wlambda)
 #      Wfdinv       <- smooth.morph(h, x, WfdParinv)$Wfdobj
 #      hinv         <- monfn(x, Wfdinv)
#       b            <- (rangeval[2]-rangeval[1])/(hinv[n]-hinv[1])
#       a            <- rangeval[1] - b*hinv[1]
#       hinv         <- a + b*hinv
#       hinv[c(1,n)] <- rangeval
 #  } else {
#	hinvbasis =    create.bspline.basis(range(h),wnbasis, worder)
#	WfdinvPar =  fdPar(hinvbasis,lambda=wlambda)
  horder = 5 # first derivatives needed
   hnbasis = nhknots + horder -2
   hbasis =  create.bspline.basis(range(h), hnbasis, horder)
  WfdinvPar = fdPar(hbasis,lambda=hlambda)

       hinvfd       <- smooth.basis(h, x, WfdinvPar)$fd
       hinv         <- as.vector(eval.fd(h, hinvfd))
       b            <- (rangeval[2]-rangeval[1])/(hinv[n]-hinv[1])
       a            <- rangeval[1] - b*hinv[1]
       hinv         <- a + b*hinv
       hinv[c(1,n)] <- rangeval
#      deltahinv <- diff(hinv)
 #      if (any(deltahinv <= 0)) stop(
 #          paste("Non-increasing warping function estimated for curve",icurve))
	for (iland in 1:nlandm) {
		l = ximarks[icurve,iland]
	if (l > max(h) || l < min(h) ) {
		landerr[icurve,iland] = NA
	}
	else {
		landerr[icurve,iland] = eval.fd(l,hinvfd ) - x0marks[iland]
	}
	}
 #   }

    #  compute registered curves

#    if (length(dim(coef)) == 2) {
#        #  single variable case
#        yregfd <- smooth.basis(hinv, y[,icurve], fdParobj)$fd
#        yregmat[,icurve] <- eval.fd(x, yregfd)
#    }
#    if (length(dim(coef)) == 3) {
#        #  multiple variable case
#        for (ivar in 1:nvar) {
#            # evaluate curve as a function of h at sampling points
#            yregfd <- smooth.basis(hinv, y[,icurve,ivar], fdParobj)$fd
#            yregmat[,icurve,ivar] <- eval.fd(x, yregfd)
#        }
#     }
  }

  cat("\n")

  #  create functional data objects for the registered curves

#  fdParobj    <- fdPar(basisobj, 2, ylambda)
#  regfdobj    <- smooth.basis(x, yregmat, fdParobj)$fd
#  regnames    <- fdobj$fdnames
#  names(regnames)[3] <- paste("Registered",names(regnames)[3])
#  regfdobj$fdnames <- regnames

  #  create functional data objects for the warping functions
   horder = 5 # first derivatives needed
   hnbasis = nhknots + horder -2
   hbasis =  create.bspline.basis(rangex, hnbasis, horder)
  fdParobj = fdPar(hbasis,lambda=hlambda)
  warpfdobj             <- smooth.basis(x, hfunmat, fdParobj)$fd
 # warpfdnames           <- fdobj$fdnames
 # names(warpfdnames)[3] <- paste("Warped",names(regnames)[1])
  #warpfdobj$fdnames     <- warpfdnames      
  
#  Wfd <- fd(Wcoef, wbasis)
  logvelfd =  smooth.basis(x, - log(eval.fd(x,warpfdobj ,1)) ,fdParobj)$fd

#  return( list( "warpfd" = warpfdobj, "Wfd" = Wfd) )
  return( list( "warpfd" = warpfdobj,  "logvelfd" = logvelfd, "land" = x0marks, "landerr" = landerr,
"hfunmat" = hfunmat, "x"= x) )
}

# to register functions on warpfdobj just do the following. 
# First smooth curves separately on their own time interval, say y_orig_fd
# then smooth.basis( x, eval.fd(h, y_orig_fd), fdParobj)
# but paying attention to rounding errors that may bring range(h) outside the def interval of y_orig_fd
# (not h but hfunmat)

# - log D h = lh
#  log D h = - lh
#  D h = exp ( - hl )
# h = int ()

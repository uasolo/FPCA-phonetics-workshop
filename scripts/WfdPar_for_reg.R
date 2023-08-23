WfdPar_for_reg <- function(marks, lambda = 1e-10) {
  Wbasis <- create.bspline.basis(breaks = marks)
  # Lfdobj <- 2
  # nOrder <- 2 + Lfdobj
  # nBasis <- length(marks) + nOrder - 2
  # Wbasis <- create.bspline.basis(rangeval = range(marks), nbasis = nBasis, norder = nOrder)
  Wfd <- fd(matrix(0,Wbasis$nbasis,1), Wbasis)
  return(fdPar(Wfd, 2, lambda))
}
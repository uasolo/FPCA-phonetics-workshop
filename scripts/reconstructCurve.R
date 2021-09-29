# Reconstruct a curve using the FPCA reconstruction formula

# aux function: reconstruct spline coef given fpcaObj, PC, dimension (where applicable) and scores.
# scores is a vector whose indices are PC indices
# for 1-D curves `dimension` is ignored
# Used internally in reconstructCurve() and reconstructDurations()
reconstructCoef <- function(fpcaObj, scores, dimension = 1) {
  # cases for 1-D or multi-D
  d <- fpcaObj$meanfd$coefs %>% dim() %>% length()
  if (! d %in% c(2,3)) {stop(paste("fpcaObj coefficients have illegal dimension:", d))}
  if (d  == 2) { # 1-D
  fpcaObj$meanfd$coefs[, 1] + 
    sapply(seq_along(scores), function(PCidx) {
      scores[PCidx] * fpcaObj$harmonics$coefs[, PCidx]
    }) %>% apply(1, sum)
  } else if (d  == 3) { # multi-D
    fpcaObj$meanfd$coefs[, 1, dimension] + 
      sapply(seq_along(scores), function(PCidx) {
        scores[PCidx] * fpcaObj$harmonics$coefs[, PCidx, dimension]
      })
  }
}


# helper function: reconstruct a curve given fpcaObj, PC, dimension (where applicable), scores and time samples
# scores is a vector whose indices are PC indices
# for 1-D curves `dimension` is ignored
reconstructCurve <- function(fpcaObj, scores, tx, dimension = 1) {
  reconstructCoef(fpcaObj, scores, dimension) %>% 
    fd(coef = ., basisobj = fpcaObj$meanfd$basis) %>%
    eval.fd(tx, .) %>% 
    as.numeric
}


# helper function: reconstruct durations given fpcaObj, PC, dimension, scores and time samples,
# landmarks and lograte_scale
# lograte_dimension should point to the dim of fpcaObj corresponding to log rate - ignored if lograte is the only dimension (1-D curves)
# scores is a vector whose indices are PC indices
# lograte_scale is the scale factor that was applied to the lograte curves
# prior to constructing the milti-dim fd object (if applied) (in ex2D.R it is called w2)
# land is the vector of registered landmarks, typically from reg$land
# Returns a vector of durations
reconstructDuration <- function(fpcaObj, scores, lograte_dimension, tx, land, lograte_scale = 1) {
  # reconstruct h(t)
  h_fd <- reconstructCoef(fpcaObj, scores, lograte_dimension) %>% 
    `*`(-1/lograte_scale) %>%
    fd(coef = ., basisobj = fpcaObj$meanfd$basis) %>%
    eval.fd(tx, .) %>%
    as.numeric %>%
    exp %>%
    smooth.basis(tx, ., fdParObj) %>%
    .$fd
  # compute durations integrating h(t) between pairs of adjacent landmarks
  sapply(1:(length(land)-1), function (Interval) {
    defint.fd(h_fd, c(land[Interval],land[Interval+1]))
  }, simplify = TRUE)
}


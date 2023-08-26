landmarkreg_timeSamples <- function(timeSamples, inputMarks, targetMarks, WfdPar=NULL, wlambda=1e-14) {
  # Landmark registration applied directly to the input time axis samples timeSamples.
  #
  # ARGUMENTS
  # timeSamples: vector (not matrix) of time input samples. 
  # They need to be included in range(inputMarks).
  # They don't need to be ordered (though it is strange if they are not).
  # Other arguments: see one_landmarkreg_nocurves()
  # 
  # VALUE
  # h(timeSamples), i.e. a vector of the same length as timeSamples with time-registered values
  # corresponding to timeSamples
  
  h <- one_landmarkreg_nocurves(inputMarks, targetMarks, WfdPar, wlambda)
  hinv_fd <- Data2fd(h, seq(0, max(targetMarks), length.out = length(h)))
  # hinv_fd <- smooth.basis(h,
  #                         seq(0, max(targetMarks), length.out = length(h)),
  #                         create.bspline.basis(breaks = h)
  # )
  return(eval.fd(timeSamples, hinv_fd))
}
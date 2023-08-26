landmarkreg_nocurves <- function(inputMarks, targetMarks=NULL, hinv=FALSE, njobs=1,
                                 WfdPar=NULL, wlambda=1e-14) {
  # Applies one_landmarkreg_nocurves to each row of inputMarks.
  # Returns all h(t) and log rates in fd form. 
  # The basis for these fd objects is set roughly = number of targetMarks.
  # If hinv=TRUE, returns also all h_inv(t), as well as a matrix of errors
  # as difference between targetMarks and the actual registered positions.
  
  # ARGUMENTS
  # inputMarks: matrix or dataframe, each row obeys the same conventions as the namesake
  # argument of one_landmarkreg_nocurves
  # targetMarks: see one_landmarkreg_nocurves.
  # If NULL, it is set as the mean by column of inputMarks.
  # WfdPar, wlambda: see one_landmarkreg_nocurves
  # hinv: if TRUE return h_inv and landmark errors
  # njobs: if > 1, use a parallel procedure 
  
  # VALUE
  # A list with the following keys:
  # h: a fd object with as many h(t) as the number of rows of inputMarks
  # lograte: - log (dh/dt) 
  # hinv: (optional) a list of fd objects, one for each rows of inputMarks
  # landerr: (optional) a matrix of the same dimensions as inputMarks reporting
  # the difference between targetMarks and the actual position of hinv(inputMarks)
  
  inputMarks <- as.matrix(inputMarks)
  # checks
  l1_not_zeros <- which(inputMarks[,1] != 0)
  if (length(l1_not_zeros) > 0) {
    stop(paste("row(s)",
               paste(l1_not_zeros, collapse = ", "),
               "of inputMarks have nonzero first landmarks"))
  }
  non_increasing <- which(t(apply(inputMarks, 1, diff)) <= 0, arr.ind = TRUE)
  if (nrow(non_increasing) > 0) {
    print(non_increasing)
    stop("The inputMarks positions above are not strictly increasing")
  }
  na <- which(is.na(inputMarks), arr.ind = TRUE)
  if (nrow(na) > 0) {
    print(na)
    stop("The inputMarks positions above contain NAs")
  }
  
  if (!is.null(targetMarks)) {
    if (targetMarks[1] != 0) {
      stop("The first targetMarks should be zero.")
    }
    non_increasing <- which(diff(targetMarks) <= 0)
    if (length(non_increasing) > 0) {
      stop(paste(
        "The targetMarks position(s)",
        paste(non_increasing, collapse = ", "),
        "are not strictly increasing"
      ))
    }
    na <- which(is.na(targetMarks))
    if (length(na) > 0) {
      stop(paste(
        "The targetMarks position(s)",
        paste(na, collapse = ", "),
        "are NAs"
      ))
    }
  }
  
  
}
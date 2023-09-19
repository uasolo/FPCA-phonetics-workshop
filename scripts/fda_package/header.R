library(fda)
if (packageVersion('fda') < 5.1) {
  stop("Please install a recent version of 'fda', 5.1 or newer.")
}
library(tidyverse)
library(emmeans)
library(abind)


dir.create('plots', showWarnings = FALSE)


getPCscores <- function(fpcaObj) {
  if (is.matrix(fpcaObj$scores)) {
    fpcaObj$scores
  } else if (is.array(fpcaObj$scores) & length(dim(fpcaObj$scores)) == 3) {
    apply(fpcaObj$scores, 1:2, sum)
  } else {
    stop("getPCscores: invalid input fpcaObj")
  }
}

source(file.path(scripts_dir,'fdPar.R'))
# use pca.fd version from package fda_2.2.5.tar.gz or earlier (you find a copy in the scripts/ dir)
# source(file.path(scripts_dir,'pca.fd.R'))
# this is a modified version of the landmarkreg() command 
source(file.path(scripts_dir,'landmarkreg.nocurve.R'))
# this is a slightly modified version of the plot.pca.fd() command,
# but you can also use the standard one.
# source(file.path(scripts_dir,'plot.pca.fd.corr.R'))
# function that computes definite integrals
source(file.path(scripts_dir,'defint.fd.R'))
source(file.path(scripts_dir,'reconstructCurve.R'))




library(fda)
packageVersion('fda')
# ‘5.1.5.1’
library(tidyverse)
library(emmeans)
library(abind)


# set your own home directory here:
home_dir = "/homes/m.gubian/Dokumente/Courses/Graz2020/FPCA-phonetics-workshop-master"
data_dir = file.path(home_dir,'data')
scripts_dir =  file.path(home_dir,'scripts')
plots_dir =  file.path(home_dir,'plots')
dir.create(plots_dir, showWarnings = FALSE)


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




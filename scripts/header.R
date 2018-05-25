library(fda)
library(lattice)

# set your own home directory here:
home_dir = "C:/Users/Michele/Dropbox/FDA/Workshop Kent 2017"
data_dir = file.path(home_dir,'data')
scripts_dir =  file.path(home_dir,'scripts')
plots_dir =  file.path(home_dir,'plots')
dir.create(plots_dir, showWarnings = FALSE)

# use pca.fd version from package fda_2.2.5.tar.gz or earlier (you find a copy in the scripts/ dir)
source(file.path(scripts_dir,'pca.fd.R'))
# this is a modified version of the landmarkreg() command 
source(file.path(scripts_dir,'landmarkreg.nocurve.R'))
# this is a slightly modified version of the plot.pca.fd() command,
# but you can also use the standard one.
source(file.path(scripts_dir,'plot.pca.fd.corr.R'))
# function that computes definite integrals
source(file.path(scripts_dir,'defint.fd.R'))




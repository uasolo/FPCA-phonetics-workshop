# testing landmark reg
library(fda)
library(tidyverse)
source('WfdPar_for_reg.R')
source('one_landmarkreg_nocurves.R')
source('landmarkreg_timeSamples.R')

inputMarks <- c(0, 1,  2) * 100
targetMarks <-  c(0, 1.9,  2) * 10
h <- one_landmarkreg_nocurves(inputMarks, targetMarks)

plot(seq(targetMarks[1], targetMarks[length(targetMarks)], length.out = length(h)),
     h,
     type = 'l')
abline(v = targetMarks, col ='red')

t_reg <- landmarkreg_timeSamples(10*(0:20), inputMarks, targetMarks)
points(t_reg, 10*(0:20), col = 'blue')



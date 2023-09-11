# testing landmark reg
library(fda)
library(tidyverse)
library(parallel)
library(doParallel)
source('WfdPar_for_reg.R')
source('one_landmarkreg_nocurves.R')
source('landmarkreg_timeSamples.R')
source('landmarkreg_nocurves.R')

inputMarks <- c(0, 1,  2) * 100
targetMarks <-  c(0, 1.9,  2) * 10
h <- one_landmarkreg_nocurves(inputMarks, targetMarks)

plot(seq(targetMarks[1], targetMarks[length(targetMarks)], length.out = length(h)),
     h,
     type = 'l')
abline(v = targetMarks, col ='red')

t_reg <- landmarkreg_timeSamples(10*(0:20), inputMarks, targetMarks)
points(t_reg, 10*(0:20), col = 'blue')

inputMarks <- matrix(c(0, 0.5, 1, 2, 0, 0.7, 1.2, 1.9, 0, 0.4, 1.1, 2.2), nrow=3, byrow = T)

for (i in 1:5) inputMarks <- rbind(inputMarks, inputMarks)

system.time({
reg <- landmarkreg_nocurves(inputMarks, njobs = 1)
})

ex <- 1 # change according to ex number
curves <- read_csv(file.path("../data/", paste("ex1D", ex, "csv", sep = '.')))
curves3 <- curves %>% 
  filter(curveId %in% 1:3)

curves3 <- curves3 %>% filter(
  curveId == 1 & time < 1.5 | curveId == 2 & time < 1.7 | curveId == 3
)

curves3 <- curves3 %>% mutate(curveId = as.factor(curveId))


land3 <- tribble(
~curveId,  ~l1, ~l2, ~l3,
  1, 0, 1, 1.5,
  2, 0, 1.7/2, 1.7,
  3, 0, 4/3, 2
) %>% 
  mutate(curveId = as.factor(curveId))

ggplot(curves3) +
  aes(time, y, color = curveId) +
  facet_grid(curveId ~ .) +
  geom_line() +
  geom_vline(aes(xintercept = value),
             data = land3 %>% pivot_longer(cols = starts_with("l"))) +
  mytheme



reg <- landmarkreg_nocurves(land3 %>% select(!curveId), c(0,1,2), compute_hinv = TRUE)
reg$h %>% plot()

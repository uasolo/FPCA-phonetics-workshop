# Exercise on GAMMs
# Author: Michele Gubian
# Last revision: March 2022

library(tidyverse)
library(magrittr)
library(mgcv)
library(itsadug)
library(tictoc)

Category.colors <- c("slategray4", "orangered")


# data_dir <- "C:/Users/Michele/Dropbox/scambio_temp/work/FDA/FPCA-phonetics-workshop/data"
data_dir <- "/vdata/ERC2/FPCA/FPCA-phonetics-workshop/data/"
ex <- 2
curves <- read_csv(file.path(data_dir, paste("ex1D", ex, "csv", sep = '.'))) %>%
  mutate(across(c(curveId,  Category), ~ factor(.x))) 

nCurves <- curves %>% select(curveId) %>% n_distinct()

# plot a few curves
ggplot(curves %>% filter(curveId %in% sample(nCurves, 20))) +
  aes(x = time, y = y, group = curveId, color = Category) +
  geom_line() +
  scale_color_manual(values = Category.colors) +
  theme_light() +
  theme(text = element_text(size = 15),
        legend.position = "bottom")

# GAM
mod <- bam(y ~ Category + s(time, by = Category, k= 10)
           # + s(time, curveId, bs = "fs", m=1, k = 15)
           ,rho = rho, AR.start = curves$time == 0
           # ,discrete=TRUE, family="scat"
           # , nthreads = 4
           , data = curves)
summary(mod)
# plot(mod, pages = 1)
plot_smooth(mod, view = "time", plot_all = "Category", col = Category.colors, rug = FALSE)
op <- par(mfrow=c(2,2))
gam.check(mod)
par(op)

rho <- acf_resid(mod)[2]

# binary smooth version
curves$IsPEAK <- (curves$Category == "PEAK") * 1
curves$IsWIDE_PEAK <- (curves$Category == "WIDE_PEAK") * 1
curves$IsLATE_PEAK <- (curves$Category == "LATE_PEAK") * 1
curves$IsUP_PEAK <- (curves$Category == "UP_PEAK") * 1

tic()
mod.bin <- bam(y ~ s(time, k = 15) + s(time, by = IsPEAK, k = 15) 
               + s(time, subjId, bs = "fs", m=1, k = 15)
               + s(time, subjId, by = IsPEAK, bs = "fs", m=1, k = 15)
               , nthreads = 4
               , data = curves)
toc()
summary(mod.bin)

get_modelterm(mod.bin, select = 2) %>%
  ggplot() +
  aes(x = time, y = fit) +
  geom_line() +
  geom_ribbon(mapping = aes(ymin = fit - se.fit, ymax = fit + se.fit), fill = "red", alpha = 0.5) +
  theme_light() +
  theme(text = element_text(size = 15))

inspect_random(mod.bin, select = 4, cond=list(subjId = 1:5 %>% factor()),
               # col = "red",
               col = colorblind_pal()(5),
               lty = 1, lwd = 2, print.summary = T)



### 2D

ex <- 2
curves <- read_csv(file.path(data_dir, paste("ex2D", ex, "csv", sep = '.'))) %>%
  mutate(across(c(curveId,  Category), ~ factor(.x)))

nCurves <- curves %>% select(curveId) %>% n_distinct()

# plot a few curves
ggplot(curves %>% filter(curveId %in% sample(nCurves, 20)) %>%
         pivot_longer(starts_with("y"), names_to = "DIM", values_to = "y")) +
  aes(x = time, y = y, group = curveId, color = Category) +
  facet_grid(DIM ~ .) +
  geom_line() +
  scale_color_manual(values = Category.colors) +
  theme_light() +
  theme(text = element_text(size = 15),
        legend.position = "bottom")




curves %<>% mutate(IsSHALLOW = 1 * (Category == "SHALLOW_y1_TROUGH"))

mod2D <- gam(
  # list(y1 ~ s(time) + s(time, by = IsSHALLOW, k = 20),
  #      y2 ~ s(time) + s(time, by = IsSHALLOW, k = 20))
  list(y1 ~ Category + s(time, by=Category, k = 20) + s(curveId, bs = "re"),
       y2 ~ Category + s(time, by=Category, k = 20) + s(curveId, bs = "re"))
  , family = mvn(d=2)
  # ,rho = AR1, AR.start = curves$time == 0
  , data = curves
)
summary(mod2D)

# plot_smooth(mod2D, view = "time", plot_all = "Category", col = Category.colors, rug = FALSE)
plot(mod2D, pages = 1)

# Dim factor trick
# Wieling, Martijn, et al.
# "Investigating dialectal differences using articulography." 
# Journal of Phonetics 59 (2016): 122-143.



curves %<>%
  pivot_longer(y1:y2, names_to = "Dim", values_to = "y") %>%
  mutate(Dim = factor(Dim)) %>%
  arrange(curveId, Dim, time) 
curves$CategoryDim <- interaction(curves$Category, curves$Dim)

mod <- bam(y ~ CategoryDim + s(time, by=CategoryDim, k = 20) 
           # + s(time, by=Dim, k = 20)
           # + s(time, curveId, bs = "fs", m=1, k = 20)
           # + s(curveId, bs = "re")
           # ,rho = AR1, AR.start = curves$time == 0
           , data = curves
           # ,discrete = TRUE, nthreads = 8
)
summary(mod)


# Exercise on GAMMs
# Author: Michele Gubian
# Last revision: March 2022

library(tidyverse)
library(magrittr)
library(mgcv)
library(itsadug)

Category.colors <- c("slategray4", "orangered")


# data_dir <- "C:/Users/Michele/Dropbox/scambio_temp/work/FDA/FPCA-phonetics-workshop/data"
data_dir <- "/vdata/ERC2/FPCA/FPCA-phonetics-workshop/data/"
ex <- 5
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
mod <- bam(y ~ Category + s(time, by = Category, k= 20)
           # + s(time, curveId, bs = "fs", m=1, k = 15)
           # ,rho = AR1, AR.start = curves$time == 0
           # ,discrete=TRUE, family="scat"
           # , nthreads = 4
           , data = curves)
summary(mod)
# plot(mod, pages = 1)
plot_smooth(mod, view = "time", plot_all = "Category", col = Category.colors)
op <- par(mfrow=c(2,2))
gam.check(mod)
par(op)
# binary smooth version
curves$IsPEAK <- (curves$Category == "PEAK") * 1
curves$IsWIDE_PEAK <- (curves$Category == "WIDE_PEAK") * 1
curves$IsLATE_PEAK <- (curves$Category == "LATE_PEAK") * 1

mod.bin <- bam(y ~ s(time) + s(time, by = IsLATE_PEAK, k = 20), data = curves)
summary(mod.bin)

get_modelterm(mod.bin, select = 2) %>%
  ggplot() +
  aes(x = time, y = fit) +
  geom_line() +
  geom_ribbon(mapping = aes(ymin = fit - se.fit, ymax = fit + se.fit), fill = "red", alpha = 0.5) +
  theme_light() +
  theme(text = element_text(size = 15))

AR1 <- acf_resid(mod)[2]

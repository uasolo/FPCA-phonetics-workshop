library(fda)
library(funData)
library(fpca)
library(tidyverse)
library(landmarkregUtils)
library(emmeans)
library(mgcv)
library(itsadug)

mytheme <- theme_light() +
  theme(text = element_text(size = 16))

Category.colors <- c("slategray4", "forestgreen")


plots_dir <- "presentations/plots/"
data_dir <- "data/"

ex <- 1 # change according to ex number
curves <- read_csv(file.path(data_dir, paste("ex1D", ex, "csv", sep = '.'))) %>% 
  mutate(across(c(curveId, Category), ~ factor(.x)))
nCurves <- curves %>% distinct(curveId) %>% nrow()




# plot a few curves
pl <- ggplot(curves %>% filter(curveId %in% sample(nCurves, 20))) +
  aes(x = time, y = y, group = curveId, color = Category) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("ex1D", ex, "curves", 'png', sep = '.')), pl,
       width = 1800, height = 1600, units = "px"
)

# build a funData object
curvesFun <- long2irregFunData(curves, id = "curveId", time = "time", value = "y") %>% 
    as.funData()

# Compute FPCA
fpca <- PACE(curvesFun, npc = 2)

# Prop of explained var
fpca$values  / sum( fpca$values)

# scores st. dev.
sdFun <- fpca$scores %>% apply(2, sd) 
# PC curves to be plotted
PCcurves <- expand_grid(PC = 1:2,
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, fractionOfStDev) %>%
  reframe(
    funData2long1(fpca$mu + fractionOfStDev * sdFun[PC] * fpca$functions[PC],
                  time = "time", value = "y")
  )
# Plot
pl <- ggplot(PCcurves) +
  aes(x = time, y = y, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  facet_wrap(~ PC, nrow = 1,
             labeller = labeller(PC = ~ str_glue("PC{.x}"))) +
  labs(color = expression(frac(s[k], sigma[k]))) +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0), color = 'black', linewidth = 1.2) +
  mytheme +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("ex1D", ex, "FPCA_curves", 'png', sep = '.')), pl,
       width = 1800, height = 1200, units = "px"
)

# collect PC scores
PCscores <- fpca$scores %>%
  `colnames<-`( paste0("s", 1:2)) %>%
  as_tibble() %>%
  bind_cols(curves %>% distinct(curveId, Category), .)

# scatterplot PC scores s1 and s2 by Category
pl <- ggplot(PCscores) +
  aes(x = s1, y = s2, color = Category) +
  geom_point() +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "right")

ggsave(file.path(plots_dir, str_c("ex1D", ex, "PCscores_scatter", 'png', sep = '.')), pl,
       width = 1800, height = 1200, units = "px"
)


# model s1
mod <- lm(s1 ~ Category, data = PCscores)
emmeans(mod, pairwise ~ Category)

# reconstruct predicted curves
emm <- emmeans(mod, pairwise ~ Category)$emmeans %>%
  as_tibble() %>% 
  select(Category, emmean) %>% 
  rename(s1= emmean)

predCurves <- emm %>% 
  group_by(Category) %>% 
  reframe(funData2long(fpca$mu + s1 * fpca$functions[1])) %>%
  select(!id) %>% 
  rename(y = value)

pl <- ggplot(predCurves) +
  aes(time, y, color = Category) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("ex1D", ex, "pred_curves", 'png', sep = '.')), pl,
       width = 1500, height = 1200, units = "px"
)

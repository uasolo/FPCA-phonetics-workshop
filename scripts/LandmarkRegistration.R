library(fda)
library(funData)
library(MFPCA)
library(tidyverse)
library(emmeans)
library(landmarkregUtils)

mytheme <- theme_light() +
  theme(text = element_text(size = 16))
        # legend.position = "bottom")

Category.colors <- c("slategray4", "orangered")

plots_dir <- "presentations/plots/"

# build a dataset with curves with 2 peaks, different durations and time distorted.
# Distort placing landmark on the two peaks, on the trough
# (and at the end) and moving them at random.

data_dir <- "data/"
ex <- 1
curves <- read_csv(file.path(data_dir, paste("ex1D", ex, "csv", sep = '.')))
curves <- curves %>% 
  filter(Category == "PEAK") %>% 
  select(!Category) %>% 
  mutate(curveId = factor(curveId))

nCurves <- curves %>% distinct(curveId) %>% nrow()

refLand <- c(0, 0.5, 1.25, 1.5, 2)
land <- curves %>% 
  distinct(curveId) %>% 
  mutate(l1 = 0,
         l2 = refLand[2] - 0.1 + abs(rnorm(nCurves, 0, 0.2)),
         l3 = l2 + refLand[3] - refLand[2] -0.1 + abs(rnorm(nCurves, 0, 0.2)),
         l4 = l3 + refLand[4] - refLand[3] -0.1 + abs(rnorm(nCurves, 0, 0.2)),
         l5 = l4 + refLand[5] - refLand[4] -0.1 + abs(rnorm(nCurves, 0, 0.2)))

# reverse landmark registration
curves <-  curves %>% 
  inner_join(land, by = "curveId") %>% 
  group_by(across(c(curveId, starts_with("l")))) %>% 
  mutate(time = landmarkreg_timeSamples(time,
                                        refLand,
                                        cur_group() %>% 
                                          select(starts_with("l")) %>% 
                                          as.numeric())) %>% 
    ungroup() %>% 
    select(!starts_with("l"))

# test
i <- 6
land_i <- land %>%
  filter(curveId == i) %>% 
  select(!curveId) %>% 
  as.numeric()
pl <- ggplot(curves %>% filter(curveId == i)) +
  aes(time, y) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = land_i, color = 'red') +
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = land_i,
                                         labels = str_c("l", seq_along(land_i)))) +
  mytheme

ggsave(file.path(plots_dir, str_c("one_curve_lands", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)


pl <- ggplot(curves %>% filter(curveId %in% 1:4)) +
  aes(time, y, group = curveId, color = curveId) +
  geom_line(linewidth = 1) +
  mytheme +
  theme(legend.position = 'none')

ggsave(file.path(plots_dir, str_c("few_unreg_curves", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

reg <- landmarkreg_nocurves(inputMarks = land %>% select(!curveId),
                            njobs = 3)
curvesReg <- applyReg(dat = curves, reg = reg,
                      id = "curveId", time = "time", value = "y")

pl <- ggplot(curvesReg %>% filter(curveId %in% 1:4)) +
  aes(time, y, group = curveId, color = curveId) +
  geom_line(linewidth = 1) +
  mytheme +
  theme(legend.position = 'none')

ggsave(file.path(plots_dir, str_c("few_reg_curves", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

### ExLand1
ex <- 1
curves <- read_csv(file.path(data_dir, paste("exLand", ex, "csv", sep = '.')))
land <- read_csv(file.path(data_dir,  paste("exLand", ex, "land", "csv", sep = '.')))

curves <- curves %>% mutate(curveId = factor(curveId))
land <- land %>% mutate(curveId = factor(curveId))

curveSample <- c(1:3, 51:53)
pl <- ggplot(curves %>% filter(curveId %in% curveSample)) +
  aes(time, y, group = curveId, color = Category) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  geom_vline(data = land %>%
               filter(curveId %in% curveSample) %>%
               pivot_longer(cols = starts_with("l")),
             mapping = aes(xintercept = value, color = Category)) +
  mytheme +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("early_late_unreg_curves", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

reg <- landmarkreg_nocurves(inputMarks = land %>% select(starts_with("l")),
                            njobs = 2)
curvesReg <- applyReg(dat = curves, reg = reg, grid = seq(0, 2, by = 0.01),
                      id = "curveId", time = "time", value = "y") %>% 
  left_join(curves %>% distinct(curveId, Category), by = "curveId")

pl <- ggplot(curvesReg %>% filter(curveId %in% curveSample)) +
  aes(time, y, group = curveId, color = Category) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  geom_vline(xintercept = reg$landmarks) + 
  xlab("registered time") +
  mytheme +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("early_late_reg_curves", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

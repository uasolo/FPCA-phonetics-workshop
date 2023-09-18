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
curves <- read_csv(file.path(data_dir, paste("exLand", ex, "csv", sep = '.'))) %>% 
  mutate(curveId = factor(curveId))

land <- read_csv(file.path(data_dir, paste("exLand", ex, "land", "csv", sep = '.'))) %>% 
  mutate(curveId = factor(curveId))



# curves <- curves %>% 
#   filter(Category == "PEAK") %>% 
#   select(!Category) %>% 
#   mutate(curveId = factor(curveId))
# 
# nCurves <- curves %>% distinct(curveId) %>% nrow()
# 
# refLand <- c(0, 0.5, 1.25, 1.5, 2)
# land <- curves %>% 
#   distinct(curveId) %>% 
#   mutate(l1 = 0,
#          l2 = refLand[2] - 0.1 + abs(rnorm(nCurves, 0, 0.2)),
#          l3 = l2 + refLand[3] - refLand[2] -0.1 + abs(rnorm(nCurves, 0, 0.2)),
#          l4 = l3 + refLand[4] - refLand[3] -0.1 + abs(rnorm(nCurves, 0, 0.2)),
#          l5 = l4 + refLand[5] - refLand[4] -0.1 + abs(rnorm(nCurves, 0, 0.2)))
# 
# # reverse landmark registration
# curves <-  curves %>% 
#   inner_join(land, by = "curveId") %>% 
#   group_by(across(c(curveId, starts_with("l")))) %>% 
#   mutate(time = landmarkreg_timeSamples(time,
#                                         refLand,
#                                         cur_group() %>% 
#                                           select(starts_with("l")) %>% 
#                                           as.numeric())) %>% 
#     ungroup() %>% 
#     select(!starts_with("l"))




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
  geom_vline(xintercept = reg$landmarks) +
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = str_c("l", seq_along(reg$landmarks)))) +
  mytheme +
  theme(legend.position = 'none')

ggsave(file.path(plots_dir, str_c("few_reg_curves", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

durationReg <- diff(reg$landmarks)

durations <- land %>% 
  pivot_longer(cols = starts_with("l"), names_to = "rightBoundary", values_to = "time") %>% 
  group_by(curveId) %>%
  mutate(before = time - lag(time)) %>% 
  filter(!is.na(before)) %>% 
  mutate(after = durationReg) %>% 
  select(!time) %>% 
  ungroup() %>% 
  pivot_longer(cols = c("before", "after"), names_to = "registration", values_to = "duration") %>% 
  mutate(registration = factor(registration, levels = c("before", "after")),
         rightBoundary = factor(rightBoundary))
  



pl <- ggplot(durations %>% filter(curveId %in% 1:4)) +
  aes(x = curveId, y = duration, color = curveId) + 
  geom_bar(stat="identity", fill = 'white') +
  geom_text(aes(label = duration %>% round(digits = 2)), size = 3.5,
            position = position_stack(vjust = 0.5), show.legend = FALSE) +
  facet_grid(~ registration) +
  ylab("time") +
  coord_flip() + 
  mytheme +
  theme(legend.position = "none")

ggsave(file.path(plots_dir, str_c("few_reg_dur", '.png')), pl,
       width = 2000, height = 1000, units = "px"
)




### landmark and fpca
ex <- 2
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
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = str_c("l", seq_along(reg$landmarks)))) +
  mytheme +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("early_late_reg_curves", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

# regular FPCA on registered y(t)

curvesFun <- long2irregFunData(curvesReg, id = "curveId", time = "time", value = "y") %>% 
  as.funData()

autoplot(curvesFun[curveSample])

fpca <- PACE(curvesFun, npc = 2)
# Prop of explained var
fpca$values / sum( fpca$values)
# scores st. dev.
sdFun <- fpca$scores %>% apply(2, sd) 
# PC curves to be plotted
PCcurves <- expand_grid(PC = 1:2,
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, fractionOfStDev) %>%
  reframe(time = fpca$mu@argvals[[1]],
          y = (fpca$mu + fractionOfStDev * sdFun[PC] * fpca$functions[PC])@X %>% as.numeric()
  )
# Plot
pl <- ggplot(PCcurves) +
  aes(x = time, y = y, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  facet_wrap(~ PC,
             ncol = 1,
             scales = "free_y",
             labeller = labeller(PC = ~ str_glue("PC{.x}"))) +
  labs(color = expression(frac(s[k], sigma[k]))) +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0), color = 'black', linewidth = 1.2) +
  geom_vline(xintercept = reg$landmarks) + 
  xlab("registered time") +
  mytheme +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("FPCA_curves_early_late_simple", '.png')), pl,
       width = 1500, height = 2000, units = "px"
)

# collect PC scores
PCscores <- fpca$scores %>% `colnames<-`( paste0("s", 1:2)) %>% as_tibble %>%
  bind_cols(curvesReg %>% distinct(curveId, Category), .)

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
  rename(`registered time` = time, y = value)

pl <- ggplot(predCurves) +
  aes(`registered time`, y, color = Category) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  geom_vline(xintercept = reg$landmarks) + 
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = str_c("l", seq_along(reg$landmarks)))) +
  mytheme +
  theme(legend.position = "bottom")

ggsave(file.path(plots_dir, str_c("early_late_pred_curves_simple", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

# duration plots
durations <- landmarks2durations(land %>% select(!Category),
                                 id = "curveId",
                                 targetMarks = reg$landmarks)



pl <- ggplot(durations %>% filter(curveId %in% c(1,51))) +
  aes(x = curveId, y = duration, color = curveId) + 
  geom_bar(stat="identity", fill = 'white') +
  geom_text(aes(label = duration %>% round(digits = 2)), size = 3.5,
            position = position_stack(vjust = 0.5), show.legend = FALSE) +
  facet_grid(~ registration) +
  scale_color_manual(values=Category.colors) +
  ylab("time") +
  coord_flip() + 
  mytheme +
  theme(legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggsave(file.path(plots_dir, str_c("early_late_reg_dur", '.png')), pl,
       width = 2000, height = 1000, units = "px"
)


# h(t)
pl <- reg$h[c(1, 51)] %>% 
  fd2funData(seq(0, reg$h$basis$rangeval[2], length.out = 100)) %>% 
  funData2long() %>% 
  ggplot(aes(time, value, group = id, color = id)) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  annotate("segment", x = reg$landmarks[2], xend = reg$landmarks[2],
           y = 0, yend = eval.fd(reg$landmarks[2], reg$h[1]),
           linetype = "dotted") +
  annotate("segment", x = reg$landmarks[2], xend = reg$landmarks[2],
           y = 0, yend = eval.fd(reg$landmarks[2], reg$h[51]),
           linetype = "dotted") +
  annotate("segment", x = 0, xend = reg$landmarks[2],
           y = eval.fd(reg$landmarks[2], reg$h[1]), yend = eval.fd(reg$landmarks[2], reg$h[1]),
           linetype = "dotted") +
  annotate("segment", x = 0, xend = reg$landmarks[2],
           y = eval.fd(reg$landmarks[2], reg$h[51]), yend = eval.fd(reg$landmarks[2], reg$h[51]),
           linetype = "dotted") +
  xlab("registered time") +
  ylab("original time") +
  mytheme +
  theme(legend.position = "none")

ggsave(file.path(plots_dir, str_c("early_late_h", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

pl <- ggplot(curves %>% filter(curveId %in% c(1,51))) +
  aes(time, y, group = curveId, color = curveId) +
  geom_line(linewidth = 1) +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = 'none')

ggsave(file.path(plots_dir, str_c("early_late_two_curves", '.png')), pl,
       width = 1500, height = 1200, units = "px"
)

pl <- bind_rows(
  h = reg$logvelfd %>% # reg$h
  fd2funData(seq(0, reg$h$basis$rangeval[2], length.out = 100)) %>% 
  funData2long(id = "curveId", time = "time", value = "y"),
  f = curvesReg %>% 
    select(!Category),
  .id = "func"
) %>% 
  mutate(func = factor(func, levels = c("f", "h"), labels = c("registered curves", "log rate(t)"))) %>% 
  filter(curveId %in% c(1,51)) %>% 
  ggplot() +
  aes(time, y, group = curveId, color = curveId) +
  geom_line(linewidth = 1) +
  facet_wrap(~ func, ncol = 1, scales = "free_y") +
  scale_color_manual(values=Category.colors) +
  xlab("registered time") +
  ylab("") +
  mytheme +
  theme(legend.position = 'none')

ggsave(file.path(plots_dir, str_c("early_late_reg_and_lograte", '.png')), pl,
       width = 1500, height = 2000, units = "px"
)

# put together a 2D fd object
yRegMult <- multiFunData(list(
  curvesFun,
  fd2funData(reg$logvelfd, argvals = curvesFun@argvals[[1]])
))

# multidim FPCA
mfpca <- MFPCA(yRegMult,
               M = 2,
               uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA"))
)
# Prop of explained var
mfpca$values  / sum( mfpca$values)




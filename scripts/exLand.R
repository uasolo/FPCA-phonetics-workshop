library(fda)
library(funData)
library(MFPCA)
library(tidyverse)
library(emmeans)

# install.packages("devtools")
# devtools::install_github("uasolo/landmarkregUtils")
library(landmarkregUtils)

mytheme <- theme_light() +
  theme(text = element_text(size = 16))

Category.colors <- c("slategray4", "orangered")


plots_dir <- "presentations/plots/"
data_dir <- "data/"

ex <- 4 # change according to ex number
curves <- read_csv(file.path(data_dir, paste("exLand", ex, "csv", sep = '.')))
land <- read_csv(file.path(data_dir,  paste("exLand", ex, "land", "csv", sep = '.')))

curves <- curves %>% mutate(across(c(curveId, Category), ~ factor(.x)))
land <- land %>% mutate(across(c(curveId, Category), ~ factor(.x)))

nCurves <- curves %>% distinct(curveId) %>% nrow()
curveSample <- sample(nCurves, 10)


# plot a few curves
ggplot(curves %>% filter(curveId %in% curveSample)) +
  aes(time, y, group = curveId, color = Category) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  # geom_vline(data = land %>%
  #              filter(curveId %in% curveSample) %>%
  #              pivot_longer(cols = starts_with("l")),
  #            mapping = aes(xintercept = value, color = Category)) +
  mytheme +
  theme(legend.position = "bottom")

# plot one curve with labelled landmarks
id <- 3
land_i <- land %>%
  filter(curveId == id) %>% 
  select(l1:last_col()) %>% 
  as.numeric()
ggplot(curves %>% filter(curveId == id)) +
  aes(time, y) +
  geom_line() +
  geom_vline(xintercept = land_i, color = 'red') +
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = land_i,
                                         labels = str_c("l", seq_along(land_i)))) +
  mytheme

# Landmark reg
reg <- landmarkreg_nocurves(inputMarks = land %>% select(starts_with("l")),
                            njobs = 2)
curvesReg <- applyReg(dat = curves, reg = reg,
                      grid = seq(0, last(reg$landmarks), length.out = 100),
                      id = "curveId", time = "time", value = "y") %>% 
  left_join(curves %>% distinct(curveId, Category), by = "curveId")

# Explore reg
reg$h %>% plot()
reg$lograte %>% plot()
reg$landmarks


# plot a few reg curves
ggplot(curvesReg %>% filter(curveId %in% curveSample)) +
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

# curve dimension
curvesFun <- long2irregFunData(curvesReg, id = "curveId", time = "time", value = "y") %>% 
  as.funData()

# put together a 2D fd object: curve and lograte
yRegMult <- multiFunData(list(
  curvesFun,
  fd2funData(reg$lograte, argvals = curvesFun@argvals[[1]])
))

yRegMult[1] %>% plot()

# multidim FPCA
mfpca <- MFPCA(yRegMult,
               M = 2,
               uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA"))
)
# Prop of explained var
mfpca$values  / sum( mfpca$values)

# scores st. dev.
sdFun <- mfpca$scores %>% apply(2, sd) 
# PC curves to be plotted
PCcurves <- expand_grid(PC = 1:2,
                        DIM = 1, # 1:2 or 1
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, DIM, fractionOfStDev) %>%
  reframe(time = mfpca$meanFunction[[DIM]]@argvals[[1]],
          y = (mfpca$meanFunction[[DIM]] +
                 fractionOfStDev * sdFun[PC] * mfpca$functions[[DIM]][PC])@X %>% as.numeric()
  )
# Plot
DIMlabels <- c(`1` = "y", `2` = "log rate")
ggplot(PCcurves) +
  aes(x = time, y = y, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  facet_grid(DIM ~ PC,
             scales = "free_y",
             labeller = labeller(PC = ~ str_glue("PC{.x}"),
                                 DIM = DIMlabels)) +
  labs(color = expression(frac(s[k], sigma[k]))) +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0), color = 'black', linewidth = 1.2) +
  # geom_vline(xintercept = reg$landmarks) + 
  xlab("registered time") +
  mytheme +
  theme(legend.position = "bottom")

# PC durations to be plotted
DIMlograte <- 2 # lograte is the second dimension in mfpca

PCdur <- expand_grid(PC = 1:2,
                     fractionOfStDev = seq(-1, 1, by=.5),
                     landmarks2long(reg$landmarks)
) %>%
  group_by(PC, fractionOfStDev, leftBoundary) %>%
  mutate(duration = lograte2duration(
    lograte = mfpca$meanFunction[[DIMlograte]] +
      fractionOfStDev * sdFun[PC] * mfpca$functions[[DIMlograte]][PC],
    from = from, to = to)) %>% 
  ungroup() %>% 
  select(!c(from, to, id))

ggplot(PCdur) +
  aes(x = fractionOfStDev %>% factor(labels = ""),
      y = duration, color = fractionOfStDev) + 
  geom_bar(stat="identity", fill = 'white') +
  geom_text(aes(label = duration %>% round(digits = 2)), size = 3.5,
            position = position_stack(vjust = 0.5), show.legend = FALSE) +
  facet_grid(~ PC, labeller = labeller(PC = ~ str_glue("PC{.x}"))) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  labs(color = expression(frac(s[k], sigma[s[k]]))) +
  ylab("time") +
  coord_flip() + 
  mytheme +
  theme(legend.position = "bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


# collect PC scores
PCscores <- mfpca$scores %>% `colnames<-`( paste0("s", 1:2)) %>% as_tibble %>%
  bind_cols(curvesReg %>% distinct(curveId, Category), .)

# scatterplot PC scores s1 and s2 by Category
ggplot(PCscores) +
  aes(x = s1, y = s2, color = Category) +
  geom_point() +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "right")


# model 
mod <- lm(s2 ~ Category, data = PCscores)
emmeans(mod, pairwise ~ Category)

# reconstruct predicted curves
emm <- emmeans(mod, pairwise ~ Category)$emmeans %>%
  as_tibble() %>% 
  select(Category, emmean) %>% 
  rename(s2= emmean)

DIMy <- 1
predCurves <- emm %>% 
  group_by(Category) %>% 
  reframe(funData2long(mfpca$meanFunction[[DIMy]] + s2 * mfpca$functions[[DIMy]][2])) %>%
  select(!id) %>% 
  rename(`registered time` = time, y = value)


# plot pred. curves
ggplot(predCurves) +
  aes(`registered time`, y, color = Category) +
  geom_line() +
  scale_color_manual(values=Category.colors) +
  geom_vline(xintercept = reg$landmarks) + 
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = str_c("l", seq_along(reg$landmarks)))) +
  mytheme +
  theme(legend.position = "bottom")

# pred. durations
predDur <- expand_grid(emm, landmarks2long(reg$landmarks)) %>% 
  group_by(Category, leftBoundary) %>% 
  mutate(duration = lograte2duration(
    lograte = mfpca$meanFunction[[DIMlograte]] +
      s2 * mfpca$functions[[DIMlograte]][2],
    from = from, to = to)) %>% 
  ungroup() %>% 
  select(!c(from, to, id))

ggplot(predDur) +
  aes(x = Category, y = duration, color = Category) + 
  geom_bar(stat="identity", fill = 'white') +
  geom_text(aes(label = duration %>% round(digits = 2)), size = 3.5,
            position = position_stack(vjust = 0.5), show.legend = FALSE) +
  ylab("time") +
  scale_color_manual(values=Category.colors) +
  coord_flip() + 
  mytheme +
  theme(legend.position = "bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

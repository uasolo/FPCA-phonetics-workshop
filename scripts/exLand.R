library(fda)
library(funData)
library(MFPCA)
library(tidyverse)
library(emmeans)
library(gridExtra)
library(RColorBrewer)

# install.packages("devtools")
# devtools::install_github("uasolo/landmarkregUtils")
library(landmarkregUtils)

mytheme <- theme_light() +
  theme(text = element_text(size = 16))

Category.colors <- c("darkslategray", "orangered")


plots_dir <- "presentations/plots"
data_dir <- "data"

ex <- 5 # change according to ex number
raw_curves <- readRDS(file.path(data_dir, str_c("ex1D", ex, "rds", sep = '.'))) %>% ungroup() %>% 
  mutate(across(c(curveId, Category), ~ factor(.x)))


sp <- 0.01 # unified sampling period
maxT <- max(raw_curves$time)
grid <- seq(0, maxT, by = sp) # unified sampling grid 
curves <- raw_curves %>% 
  group_by(curveId, Category) %>% # all the factors at the level of curveId or higher (e.g. speaker)
  reframe(approx(time, y, grid) %>% as_tibble()) %>% # linear interpolation on grid 
  ungroup() %>% 
  rename(time = x, y = y)

land <- readRDS(file.path(data_dir, str_c("land1D", ex, "rds", sep = '.'))) %>% ungroup() %>% 
  mutate(across(c(curveId, Category), ~ factor(.x)))


# plot a few curves
set.seed(123)
subset_curveId <- raw_curves %>%
  ungroup() %>% 
  distinct(curveId) %>%
  slice_sample(n = 20)

ggplot(curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = y, group = curveId, color = Category) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  mytheme  +
  theme(legend.position = "bottom")



# with landmark position
land_y <- land %>%
  pivot_longer(l1:last_col(), names_to = "landmark", values_to = "time") %>% 
  inner_join(subset_curveId, by = "curveId") %>% 
  group_by(curveId) %>% 
  mutate(y = {
    y <- curves %>%
      inner_join(cur_group(), by = 'curveId') %>% 
      pull(y)
    approx(grid, y, time, rule = 2)$y
  }) 
           


ggplot(curves %>% inner_join(subset_curveId, by = "curveId")) +
  aes(x = time, y = y, group = curveId) +
  geom_line(linewidth = 0.6, color = 'slategray4') +
  geom_point(data = land_y,
             mapping = aes(time, y, color = landmark, group = curveId),
             inherit.aes = FALSE,
             size = 2) +
  scale_color_brewer(palette = "Dark2") + 
  mytheme  +
  theme(legend.position = "bottom")





# plot one curve with labelled landmarks
id <- 73
land_id <- land %>% filter(curveId == id)
landmarks <- land_id %>% select(l1:last_col()) %>% as.numeric()
landlabels <- land_id %>% select(l1:last_col()) %>% colnames()
Cat_id <- land_id %>% pull(Category)

ggplot(curves %>% filter(curveId == id)) +
  aes(time, y) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = landmarks 
             , color = 'red') +
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = landmarks,
                                         labels = landlabels)) +
  ggtitle(str_glue("CurveId: {id}; Category: {Cat_id}")) +
  mytheme





# Landmark reg
# 1. compute time warping curves h(t), based on landmarks (not on curves!)
reg <- landmarkreg_nocurves(inputMarks = land %>% select(starts_with("l")),
                            njobs = 2)
# 2. apply time warping to curves
curvesReg <- applyReg(dat = curves %>% filter(!is.na(y)), # remove NAs otherwise applyReg fails
                      reg = reg,
                      grid = seq(0, last(reg$landmarks), length.out = 100),
                      id = "curveId", time = "time", value = "y") %>% 
  left_join(curves %>% distinct(curveId, Category), by = "curveId")

# Explore reg
reg$h %>% plot()
reg$lograte %>% plot()
reg$landmarks


# plot a few reg curves
ggplot(curvesReg %>% inner_join(subset_curveId, by = "curveId")) +
  aes(time, y, group = curveId, color = Category) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values=Category.colors) +
  geom_vline(xintercept = reg$landmarks) + 
  xlab("registered time") +
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = reg$landmarks %>% names())) +
  mytheme +
  theme(legend.position = "bottom")

# curve dimension
curvesFun <- long2irregFunData(curvesReg, id = "curveId", time = "time", value = "y") %>% 
  as.funData()


# Here you can jump to ex1D.R and proceed with FPCA (and/or GAM)...
# .. or continue below for joint curve/duration FPCA

# put together a 2D fd object: curve and lograte
yRegMult <- multiFunData(list(
  curvesFun,
  fd2funData(reg$lograte, argvals = curvesFun@argvals[[1]])
))

yRegMult[61] %>% plot()

# multidim FPCA
nPC <- 2
mfpca <- MFPCA(yRegMult,
               M = nPC,
               uniExpansions = list(list(type = "uFPCA"),list(type = "uFPCA"))
)
# Prop of explained var
mfpca$values  / sum( mfpca$values)

# scores st. dev.
sdFun <- mfpca$values %>% sqrt()
# PC curves to be plotted
DimCurves <- c(1)
DimLograte <- 2 # lograte is the second dimension in mfpca
PCcurves <- expand_grid(PC = 1:nPC,
                        Dim = DimCurves , #c(DimCurves), # DimLograte),
                        fractionOfStDev = seq(-1, 1, by=.25)) %>%
  group_by(PC, Dim, fractionOfStDev) %>%
  reframe(time = mfpca$meanFunction[[Dim]]@argvals[[1]],
          y = (mfpca$meanFunction[[Dim]] +
                 fractionOfStDev * sdFun[PC] * mfpca$functions[[Dim]][PC])@X %>% as.numeric()
  )
# Plot
DimLabels <- c(`1` = "y", `2` = "log rate")
PCcurves_plot <- ggplot(PCcurves) +
  aes(x = time, y = y, group = fractionOfStDev, color = fractionOfStDev) +
  geom_line() +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  facet_grid(Dim ~ PC,
             scales = "free_y",
             labeller = labeller(PC = ~ str_glue("PC{.x}"),
                                 Dim = DimLabels)) +
  labs(color = expression(frac(s[k], sigma[k]))) +
  geom_line(data = PCcurves %>% filter(fractionOfStDev == 0), color = 'black', linewidth = 1.2) +
  geom_vline(xintercept = reg$landmarks) +
  xlab("registered time") + 
  mytheme +
  theme(legend.position = "right",
        axis.title.y = element_blank())

PCcurves_plot

# PC durations to be plotted

PCdur <- expand_grid(PC = 1:nPC,
                     fractionOfStDev = seq(-1, 1, by=.5),
                     landmarks2long(reg$landmarks)
) %>%
  group_by(PC, fractionOfStDev, leftBoundary) %>%
  mutate(duration = lograte2duration(
    lograte = mfpca$meanFunction[[DimLograte]] +
      fractionOfStDev * sdFun[PC] * mfpca$functions[[DimLograte]][PC],
    from = from, to = to)) %>% 
  ungroup() %>% 
  select(!c(from, to, id))

PCdur_plot <- ggplot(PCdur) +
  aes(x = fractionOfStDev %>% factor(labels = ""),
      y = duration, color = fractionOfStDev) + 
  geom_bar(stat="identity", fill = 'white') +
  geom_text(aes(label = duration %>% round(digits = 2)), size = 5,
            position = position_stack(vjust = 0.5), show.legend = FALSE) +
  facet_grid(~ PC, labeller = labeller(PC = ~ str_glue("PC{.x}"))) +
  scale_color_gradient2(low = "blue", mid = "grey", high = "orangered") +
  labs(color = expression(frac(s[k], sigma[s[k]]))) +
  ylab("time") +
  coord_flip() + 
  mytheme +
  theme(legend.position = "right",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

PCdur_plot

grid.arrange(PCcurves_plot, PCdur_plot, nrow=2)

# Plot PC curves against original time axis + landmark positions
# (this is an example for 1-dim curves)

# add landmark position (except for the first one at time = 0)
PCdur <- PCdur %>% 
  group_by(PC, fractionOfStDev) %>% 
  mutate(land = cumsum(duration))

# compute original time samples from registered time samples
PCcurves <- PCcurves %>% 
  inner_join(PCdur %>% distinct(PC, fractionOfStDev), by = c("PC", "fractionOfStDev")) %>% 
  group_by(PC, Dim, fractionOfStDev) %>% 
  mutate(origTime = {
    landmarks <- c(0, # first landmark at time = 0
                   PCdur %>% 
                     inner_join(cur_group(), by = c("PC", "fractionOfStDev")) %>% 
                     pull(land))
    landmarkreg_timeSamples(timeSamples = time,
                            inputMarks = reg$landmarks,
                            targetMarks = landmarks)
  })



# collect PC scores
PCscores <- mfpca$scores %>%
  `colnames<-`( paste0("s", 1:nPC)) %>%
  as_tibble() %>%
  bind_cols(curves %>% distinct(curveId, Category), .)

# scatterplot PC scores s1 and s2 by Category
ggplot(PCscores) +
  aes(x = s1, y = s2, color = Category) +
  geom_point() +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "right")

# boxplots PC scores by Category
PCscores %>% 
  pivot_longer(cols = s1:all_of(str_glue("s{nPC}")), 
               names_to = "PC", values_to = "score") %>% 
  ggplot() +
  aes(x = Category, y = score, color = Category) +
  geom_boxplot() +
  facet_wrap(~ PC) +
  scale_color_manual(values=Category.colors) +
  mytheme +
  theme(legend.position = "bottom")

# Model
s <- 1 # score index
model_eq <- str_glue("s{s} ~ Category") %>% as.formula()
mod <- lm(model_eq, data = PCscores)
mod %>% summary()
emmeans(mod, pairwise ~ Category)

# reconstruct predicted curves

predCurves <- emmeans(mod, pairwise ~ Category)$emmeans %>%
  as_tibble() %>%
  expand_grid(Dim = DimCurves) %>% 
  group_by(Category, Dim) %>% 
  reframe(bind_cols(
    funData2long1(mfpca$meanFunction[[Dim]] +
                    emmean * mfpca$functions[[Dim]][s], value = "y"),
    funData2long1(mfpca$meanFunction[[Dim]] +
                    lower.CL * mfpca$functions[[Dim]][s], value = "yl") %>% 
      select(yl),
    funData2long1(mfpca$meanFunction[[Dim]] +
                    upper.CL * mfpca$functions[[Dim]][s], value = "yu") %>% 
      select(yu)
  ))

predCurves_plot <- predCurves %>% 
  ggplot() +
  aes(time, y, color = Category) +
  geom_line(linewidth = 0.8) +
  geom_ribbon(aes(x = time, ymin = yl, ymax = yu, fill = Category),
              alpha = 0.3, inherit.aes = FALSE) +
  scale_color_manual(values=Category.colors) +
  scale_fill_manual(values=Category.colors) +
  xlab("registered time") +
  geom_vline(xintercept = reg$landmarks) + 
  scale_x_continuous(sec.axis = dup_axis(name = "landmarks",
                                         breaks = reg$landmarks,
                                         labels = reg$landmarks %>% names())) +
  # ggtitle(str_glue("Reconstructed curves according to regr model: s{s} ~ Category")) +
  mytheme +
  theme(legend.position = "bottom")

predCurves_plot

# pred. durations
predDur <- expand_grid(emmeans(mod, pairwise ~ Category)$emmeans %>% 
                         as_tibble()
                       , landmarks2long(reg$landmarks)) %>% 
  group_by(Category, leftBoundary) %>% 
  mutate(duration = lograte2duration(
    lograte = mfpca$meanFunction[[DimLograte]] +
      emmean * mfpca$functions[[DimLograte]][s],
    from = from, to = to)) %>% 
  ungroup() %>% 
  select(!c(from, to, id))

predDur_plot <- ggplot(predDur) +
  aes(x = Category, y = duration, color = Category) + 
  geom_bar(stat="identity", fill = 'white') +
  geom_text(aes(label = duration %>% round(digits = 2)), size = 5,
            position = position_stack(vjust = 0.5), show.legend = FALSE) +
  ylab("time") +
  scale_color_manual(values=Category.colors) +
  coord_flip() + 
  mytheme +
  theme(legend.position = "bottom",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

predDur_plot

grid.arrange(predCurves_plot, predDur_plot, nrow=2)

# Predicted curves on the original time axis

# add landmark position (except for the first one at time = 0)
predDur <- predDur %>% 
  group_by(Category) %>% 
  mutate(land = cumsum(duration))

# compute original time samples from registered time samples
predCurves <- predCurves %>% 
  group_by(Category) %>% 
  mutate(origTime = {
    landmarks <- c(0, # first landmark at time = 0
                   predDur %>% 
                     inner_join(cur_group(), by = "Category") %>% 
                     pull(land))
    landmarkreg_timeSamples(timeSamples = time,
                            inputMarks = reg$landmarks,
                            targetMarks = landmarks)
  })

# compute (time, y) coordinates of landmark points on predicted curves
predLand <- predDur %>% 
  group_by(Category) %>% 
  mutate(land = cumsum(duration)) %>% 
  select(Category, land) %>% 
  reframe(land = list(c(0, land))) %>% 
  inner_join(predCurves %>% select(Category, origTime, y),
             by = 'Category',
             multiple = "all") %>% 
  nest(.by = c(Category, land), .key = "curves") %>% 
  group_by(Category) %>%
  mutate(y = list(approx(curves[[1]]$origTime, curves[[1]]$y, land[[1]], rule = 2)$y)) %>% 
  select(!curves) %>% 
  unnest(cols = c(land, y))
  


ggplot(predCurves) +
  aes(origTime, y, color = Category) +
  geom_line(linewidth = 0.8) +
  geom_point(data = predLand,
             mapping = aes(land, y, color = Category),
             inherit.aes = FALSE,
             size = 3
            ) +
  scale_color_manual(values=Category.colors) +
  xlab("time") +
  # ggtitle(str_glue("Reconstructed curves according to regr model: s{s} ~ Category")) +
  mytheme +
  theme(legend.position = "bottom")
